#include "AmrLevelWave.H"
#include <cmath>

using namespace amrex;

/*
void
AmrLevelWave::initData ()
{
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();

    MultiFab& S_new = get_new_data(State_Type);
    auto const& snew = S_new.arrays();

    amrex::ParallelFor(S_new,
    [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
    {
        Real x = problo[0] + (i+0.5)*dx[0];
        Real r = x - 0.5;
        constexpr Real Pi = 3.1415926535897932384626;
        snew[bi](i,j,k,0) = 0.0;
        snew[bi](i,j,k,1) = std::exp(-16.*r*r) * std::pow(std::cos(Pi*r),6);
    });
    } */

void
AmrLevelWave::initData ()
{
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();

    MultiFab& S_new = get_new_data(State_Type);
    auto const& snew = S_new.arrays();

    amrex::ParallelFor(S_new,
    [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
    {
        Real x = problo[0] + (i+0.5)*dx[0];
	Real y = problo[1] + (j+0.5)*dx[1]; //the 0.5 makes the value cell centered (i think)
	Real z = problo[2] + (k+0.5)*dx[2];
					
        Real rr2 = (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) + (z - 0.5)*(z - 0.5);  // this is the radius 
	Real rr = sqrt(rr2);

	Real amplitude = 100.;

	Real width = 0.01;
	
        constexpr Real Pi = 3.1415926535897932384626;
        snew[bi](i,j,k,0) = 1.0 + amplitude * std::exp(-rr2/width);
        snew[bi](i,j,k,1) = 0.0;
    });
}
