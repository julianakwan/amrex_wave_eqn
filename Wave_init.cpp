#include "AmrLevelWave.H"
#include "InitialConditions.H"

#include <AMReX_ParmParse.H>
#include <cmath>

using namespace amrex;

void AmrLevelWave::initData() {
  const auto problo = geom.ProbLoArray();
  const auto probhi = geom.ProbHiArray();
  const auto dx = geom.CellSizeArray();

  Real midpts[3];
  midpts[0] = 0.5 * (probhi[0] - problo[0]);
  midpts[1] = 0.5 * (probhi[1] - problo[1]);
  midpts[2] = 0.5 * (probhi[2] - problo[2]);

  MultiFab &S_new = get_new_data(State_Type);
  auto const &snew = S_new.arrays();

  amrex::Real t0 = -5.4;
  amrex::ParmParse pp("wave");
  pp.query("initial_time", t0);

  InitialConditions SineGordon(alpha, k_r);

  amrex::ParallelFor(S_new, [=] AMREX_GPU_DEVICE(int bi, int i, int j,
                                                 int k) noexcept {
    Real x = problo[0] + (i + 0.5) * dx[0];
    Real y = problo[1] + (j + 0.5) * dx[1];
    Real z = problo[2] + (k + 0.5) * dx[2];


    for (int n = 0; n < NFIELDS; n++) {
      // snew[bi](i,j,k,2*n) = 0.0;
      // snew[bi](i,j,k,2*n+1) = std::exp(-16.*rr2) *
      //std::pow(std::cos(Pi*rr2),6);

      // snew[bi](i,j,k,2*n) = 1+ampl[n]*std::exp(-width[n]*rr2);
      // snew[bi](i,j,k,2*n+1) = 0;

      // snew[bi](i, j, k, 2 * n) = SineGordon.breather_solution(x - midpts[0], 0);
      // snew[bi](i, j, k, 2 * n + 1) =
      //     SineGordon.breather_solution_deriv(x - midpts[0], 0);

      snew[bi](i,j,k,2*n) =
      SineGordon.breather_solution(x-midpts[0], y-midpts[1],
				   z-midpts[2], t0);
      snew[bi](i,j,k,2*n+1) = 0;
    }
  });
}
