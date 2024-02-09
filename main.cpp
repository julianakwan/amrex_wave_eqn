#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include "AmrWave.H"

#ifdef USE_CATALYST
#include "CatalystAdaptor.h"
#endif

using namespace amrex;

amrex::LevelBld* getLevelBld ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);



    int  max_step = -1;
    Real strt_time = 0.0;
    Real stop_time = -1.0;

   
    {
        ParmParse pp;
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);

    }

#ifdef USE_CATALYST
    ParmParse pp;
    std::string catalyst_filename; 
    std::string catalyst_options;
    int catalyst_verbosity_level; 

    pp.query("catalyst_input_script", catalyst_filename);
    pp.query("catalyst_options", catalyst_options);

    CatalystAdaptor::Initialize(catalyst_filename, catalyst_options);
#endif


    if (max_step < 0 && stop_time < 0.0) {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        auto amr = std::make_unique<AmrWave>(getLevelBld());

        amr->init(strt_time,stop_time);

	
        while ( amr->okToContinue() &&
                (amr->levelSteps(0) < max_step || max_step < 0) &&
                (amr->cumTime() < stop_time || stop_time < 0.0) )
        {
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

#ifdef USE_CATALYST
    CatalystAdaptor::Finalize();
#endif


    amrex::Finalize();
}
