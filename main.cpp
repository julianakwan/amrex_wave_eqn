#include "AmrWave.H"
#include <AMReX.H>
#include <AMReX_ParmParse.H>

#define USE_CATALYST 0

#if USE_CATALYST
#include "CatalystAdaptor.h"
#endif

using namespace amrex;

amrex::LevelBld *getLevelBld();

int main(int argc, char *argv[]) {
  amrex::Initialize(argc, argv);

  int max_step = -1;
  Real strt_time = 0.0;
  Real stop_time = -1.0;

  {
    ParmParse pp;
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);
  }

#if USE_CATALYST
  ParmParse pp("paraview");
  std::string catalyst_filename;
  std::string paraview_impl_dir;
  std::vector<std::string> catalyst_options;

  pp.query("input_script", catalyst_filename);
  pp.query("path_to_catalyst_lib", paraview_impl_dir);

  std::string nm;
  int n_opts = pp.countval("options");
  for (int i = 0; i < n_opts; i++) {
    pp.get("options", nm, i);
    catalyst_options.push_back(nm);
  }

  //  pp.query("options", catalyst_options);
  CatalystAdaptor::Initialize(catalyst_filename, catalyst_options,
                              paraview_impl_dir);
#endif

  if (max_step < 0 && stop_time < 0.0) {
    amrex::Abort(
        "Exiting because neither max_step nor stop_time is non-negative.");
  }

  {
    auto amr = std::make_unique<AmrWave>(getLevelBld());

    amr->init(strt_time, stop_time);

    while (amr->okToContinue() &&
           (amr->levelSteps(0) < max_step || max_step < 0) &&
           (amr->cumTime() < stop_time || stop_time < 0.0)) {
      amr->coarseTimeStep(stop_time);
      int current_step = amr->levelSteps(0);
    }

    // Write final checkpoint and plotfile

    if (amr->stepOfLastCheckPoint() < amr->levelSteps(0)) {
      amr->checkPoint();
    }
    if (amr->stepOfLastPlotFile() < amr->levelSteps(0)) {
      amr->writePlotFile();
    }
  }

#if USE_CATALYST
  CatalystAdaptor::Finalize();
#endif

  amrex::Finalize();
}
