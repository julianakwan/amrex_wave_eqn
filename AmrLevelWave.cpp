#include "AmrLevelWave.H"
#include <AMReX_ParmParse.H>
#include <Derive.H>
#include <numeric>

#ifdef USE_CATALYST
#include "CatalystAdaptor.h"
#endif

using namespace amrex;

// constexpr int AmrLevelWave::ncomp;

constexpr int AmrLevelWave::nghost;
int AmrLevelWave::verbose = 0;
int AmrLevelWave::rk_order = 4;
Real AmrLevelWave::cfl = 0.2;
Vector<float> AmrLevelWave::ampl;
Vector<float> AmrLevelWave::width;
Real AmrLevelWave::scalar_mass = 1.0;
int AmrLevelWave::ncomp = NFIELDS * 2;
Real AmrLevelWave::tagging_criterion = 1e3;
Real AmrLevelWave::k_r = 1.0;
Real AmrLevelWave::alpha = 1.0;
Vector<std::string> AmrLevelWave::diagnostics; // this is for error checking

namespace {
struct WaveBCFill {
  AMREX_GPU_DEVICE
  void operator()(const IntVect &iv, Array4<Real> const &dest,
                  const int /*dcomp*/, const int /*numcomp*/,
                  GeometryData const &geom, const Real /*time*/,
                  const BCRec * /*bcr*/, const int /*bcomp*/,
                  const int /*orig_comp*/) const {
    // removed because periodic

    // In this test, we only need to fill the x-direction bounary,
    // because it's periodic in other directions.  We also could
    // have used BCType::reflect_odd, and then we would not need to
    // do anything here.  However, this is an example of how to fill
    // external Dirichlet BC.
    // const int ilo = geom.Domain().smallEnd(0);
    // const int ihi = geom.Domain().bigEnd(0);
    // const auto [i,j,k] = iv.dim3();
    // if (i < ilo) {
    //     dest(i,j,k,0) = -dest(2*ilo-i-1,j,k,0);
    //     dest(i,j,k,1) = -dest(2*ilo-i-1,j,k,1);
    // }
    // if (i > ihi) {
    //     dest(i,j,k,0) = -dest(2*ihi-i+1,j,k,0);
    //     dest(i,j,k,1) = -dest(2*ihi-i+1,j,k,1);
    // }
  }
};

void wave_bcfill(Box const &bx, FArrayBox &data, const int dcomp,
                 const int numcomp, Geometry const &geom, const Real time,
                 const Vector<BCRec> &bcr, const int bcomp, const int scomp) {
  GpuBndryFuncFab<WaveBCFill> gpu_bndry_func(WaveBCFill{});
  gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}
} // namespace

AmrLevelWave::AmrLevelWave(Amr &amr, int lev, const Geometry &gm,
                           const BoxArray &ba, const DistributionMapping &dm,
                           Real time)
    : AmrLevel(amr, lev, gm, ba, dm, time) {}

AmrLevelWave::~AmrLevelWave() {}

void AmrLevelWave::variableSetUp() {
  read_params();

  desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                         StateDescriptor::Point, nghost, ncomp,
                         &cell_quartic_interp);

  // int lo_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,    // external
  // Dirichlet
  //                                        BCType::int_dir,    // periodic
  //                                        BCType::int_dir) }; // periodic
  // int hi_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,
  //                                        BCType::int_dir,
  //                                        BCType::int_dir) };
  int lo_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::int_dir,   // periodic
                                         BCType::int_dir,   // periodic
                                         BCType::int_dir)}; // periodic
  int hi_bc[BL_SPACEDIM] = {
      AMREX_D_DECL(BCType::int_dir, BCType::int_dir, BCType::int_dir)};

  Vector<BCRec> bcs(ncomp, BCRec(lo_bc, hi_bc));

  StateDescriptor::BndryFunc bndryfunc(wave_bcfill);
  bndryfunc.setRunOnGPU(true);

  Vector<std::string> param_names(ncomp); // = {"u1", "v1", "u2", "v2"};

  // amrex::Print() << ncomp <<  "\n";
  // amrex::Print() << nfields <<  "\n";

  for (int n = 0; n < NFIELDS; n++) {
    char name[6];
    sprintf(name, "phi%d", n);
    param_names[2 * n] = name;
    sprintf(name, "dphi%d", n);
    param_names[2 * n + 1] = name;
  }

  // for (int n = 0; n < ncomp; n++)
  //   amrex::Print() << param_names[n] << " " << n <<  "\n";

  desc_lst.setComponent(State_Type, 0, param_names, bcs, bndryfunc);

  // New diagnostic variable for testing interpolation between levels (against
  // analytic solution)

  derive_lst.add(
      "frac_error", amrex::IndexType::TheCellType(), 1, diagnostics,
      derive_func_fab,
      [=](const amrex::Box &box) { return amrex::grow(box, nghost); },
      &amrex::cell_quartic_interp);

  derive_lst.addComponent("frac_error", desc_lst, State_Type, 0, 1);
}

void AmrLevelWave::variableCleanUp() {
  desc_lst.clear();
  derive_lst.clear();
}

void AmrLevelWave::init(AmrLevel &old) {
  Real dt_new = parent->dtLevel(Level());
  Real cur_time = old.get_state_data(State_Type).curTime();
  Real prev_time = old.get_state_data(State_Type).prevTime();
  Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  for (int k = 0; k < NUM_STATE_TYPE; ++k) {
    MultiFab &S_new = get_new_data(k);
    FillPatch(old, S_new, 0, cur_time, k, 0, ncomp);
  }
}

void AmrLevelWave::init() {
  Real dt = parent->dtLevel(Level());
  Real cur_time = getLevel(Level() - 1).state[State_Type].curTime();
  Real prev_time = getLevel(Level() - 1).state[State_Type].prevTime();
  Real dt_old = (cur_time - prev_time) / (Real)parent->MaxRefRatio(Level() - 1);
  setTimeLevel(cur_time, dt_old, dt);

  for (int k = 0; k < NUM_STATE_TYPE; ++k) {
    MultiFab &S_new = get_new_data(k);
    FillCoarsePatch(S_new, 0, cur_time, k, 0, ncomp);
  }
}

void AmrLevelWave::computeInitialDt(int finest_level, int /*sub_cycle*/,
                                    Vector<int> &n_cycle,
                                    const Vector<IntVect> & /*ref_ratio*/,
                                    Vector<Real> &dt_level, Real stop_time) {
  if (Level() > 0) {
    return;
  } // Level 0 does this for every level.

  Vector<int> nsteps(
      n_cycle.size()); // Total number of steps in one level 0 step
  std::partial_sum(n_cycle.begin(), n_cycle.end(), nsteps.begin(),
                   std::multiplies<int>());

  Real dt_0 = std::numeric_limits<Real>::max();
  for (int ilev = 0; ilev <= finest_level; ++ilev) {
    const auto dx = parent->Geom(ilev).CellSizeArray();
    Real dtlev = cfl * std::min({AMREX_D_DECL(dx[0], dx[1], dx[2])});
    dt_0 = std::min(dt_0, nsteps[ilev] * dtlev);
  }
  // dt_0 will be the time step on level 0 (unless limited by stop_time).

  if (stop_time > 0) {
    // Limit dt's by the value of stop_time.
    const Real eps = 0.001 * dt_0;
    const Real cur_time = get_state_data(State_Type).curTime();
    if ((cur_time + dt_0) > (stop_time - eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  for (int ilev = 0; ilev <= finest_level; ++ilev) {
    dt_level[ilev] = dt_0 / Real(nsteps[ilev]);
  }
}

void AmrLevelWave::computeNewDt(int finest_level, int sub_cycle,
                                Vector<int> &n_cycle,
                                const Vector<IntVect> &ref_ratio,
                                Vector<Real> & /*dt_min*/,
                                Vector<Real> &dt_level, Real stop_time,
                                int /*post_regrid_flag*/) {
  // For this code we can just call computeInitialDt.
  computeInitialDt(finest_level, sub_cycle, n_cycle, ref_ratio, dt_level,
                   stop_time);
}

void AmrLevelWave::post_timestep(int iteration) {
  if (Level() < parent->finestLevel()) {
    auto &fine_level = getLevel(Level() + 1);
    MultiFab &S_fine = fine_level.get_new_data(State_Type);
    MultiFab &S_crse = this->get_new_data(State_Type);
    Real t = get_state_data(State_Type).curTime();

    IntVect ratio = parent->refRatio(Level());
    AMREX_ASSERT(ratio == 2 || ratio == 4);
    if (ratio == 2) {
      // Need to fill one ghost cell for the high-order interpolation below
      FillPatch(fine_level, S_fine, 1, t, State_Type, 0, ncomp);
    }

    // Original interpolation:
    FourthOrderInterpFromFineToCoarse(S_crse, 0, 2, S_fine, ratio);
    // Average between cell faces, also removes need for fill patch;
    //  average_down(S_fine, S_crse, 0, S_crse.nComp(), ratio);
  }

#ifdef USE_CATALYST

  // If using insitu, then save the grid geometry, field data and refinement
  // ratios as containers to pass into Catalyst

  const int output_levs =
      parent->finestLevel() + 1; // always have at least one level, i.e coarsest
  amrex::Vector<const amrex::MultiFab *> mfs(output_levs);
  amrex::Vector<amrex::IntVect> ref_ratios(output_levs);
  amrex::Vector<amrex::Geometry> geoms(output_levs);

  ref_ratios = parent->refRatio();

  for (int lev = 0; lev < output_levs; lev++) {
    auto &level = getLevel(lev);
    MultiFab &S = level.get_new_data(State_Type);
    mfs[lev] = &S;
    geoms[lev] = level.Geom(); // I think this is for the current level
  }

  Real time = get_state_data(State_Type).curTime();
  amrex::Print() << "Number of levels: " << output_levs << "\n";

  // names of variables we want to plot
  std::vector<std::string> catalyst_names;
  const std::list<DeriveRec> &dlist = derive_lst.dlist();
  for (auto const &d : dlist) {
    if (amrex::Amr::isDerivePlotVar(d.name())) {
      catalyst_names.push_back(d.name());
    }
  }

  for (int typ = 0; typ < desc_lst.size(); typ++) {
    for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
      if (amrex::Amr::isStatePlotVar(desc_lst[typ].name(comp)) &&
          desc_lst[typ].getType() == IndexType::TheCellType()) {
        catalyst_names.push_back(desc_lst[typ].name(comp));
      }
    }
  }

  for (int i = 0; i < catalyst_names.size(); i++)
    amrex::Print() << "Plot vars = " << catalyst_names[i] << std::endl;

  CatalystAdaptor::Execute(verbose, nStep(), time, iteration, output_levs,
                           geoms, ref_ratios, mfs);

#endif

  AmrLevel::post_timestep(iteration);
}

void AmrLevelWave::errorEst(TagBoxArray &tags, int /*clearval*/, int /*tagval*/,
                            Real /*time*/, int /*n_error_buf*/, int /*ngrow*/) {
  auto const &S_new = get_new_data(State_Type);

  const char tagval = TagBox::SET;
  auto const &a = tags.arrays();
  auto const &s = S_new.const_arrays();
  auto local_tagging_criterion = tagging_criterion;
  amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) {
    // Just an example, not necessarily good choice.
    // TODO: this is only in the first scalar field, have to check for the other
    // fields also
    if (amrex::Math::abs(s[bi](i, j, k, 1)) > local_tagging_criterion) {
      a[bi](i, j, k) = tagval;
    }
  });
}

void AmrLevelWave::read_params() {
  ParmParse pp("wave");
  pp.query("v", verbose); // Could use this to control verbosity during the run
  pp.query("rk_order", rk_order);
  pp.query("cfl", cfl);
  pp.getarr("initial_amplitude", ampl, 0, NFIELDS);
  pp.getarr("initial_width", width, 0, NFIELDS);
  pp.query("scalar_mass", scalar_mass);
  pp.query("tagging_criterion", tagging_criterion);
  pp.query("wave_vector", k_r);
  pp.query("alpha", alpha);

  ncomp = 2 * NFIELDS;

  // // read array of initial amplitudes
  // Vector<float> ampl;
  // int nx;
  // if (nx=pp.countval("initial_amplitude")) {
  //    // get nx values starting at index 0 and store in ampl.
  //    // dx is automatically resized here.
  //    pp.getarr("initial_amplitude",ampl,0,nx);
}
