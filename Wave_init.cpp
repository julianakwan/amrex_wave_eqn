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
    const auto probhi = geom.ProbHiArray();
    const auto dx = geom.CellSizeArray();

    Real midpts[3];
    midpts[0] = 0.5*(probhi[0]-problo[0]);
    midpts[1] = 0.5*(probhi[1]-problo[1]);
    midpts[2] = 0.5*(probhi[2]-problo[2]); 
				       
    amrex::Print() << midpts[0] << midpts[1] << midpts[2] << "\n";

    
    MultiFab& S_new = get_new_data(State_Type);
    auto const& snew = S_new.arrays();

    amrex::ParallelFor(S_new,
    [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
    {
        Real x = problo[0] + (i+0.5)*dx[0];
	Real y = problo[1] + (j+0.5)*dx[1];
	Real z = problo[2] + (k+0.5)*dx[2];


					
        Real rr2 = (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) + (z - 0.5)*(z - 0.5);  // this is the radius 
	//the 0.5's place the gaussian at the center because the box goes from 0.0-1.0 
	Real rr = sqrt(rr2);


	
        constexpr Real Pi = 3.1415926535897932384626;
	Real mu = 0.7;
	Real mu_coeff = mu/(std::sqrt(1-mu*mu));
	constexpr Real t0 = -5.4;


		// constexpr Real k_r = 1;
		// constexpr Real omega = 1;
	for (int n = 0; n < nfields; n++)
	  {
	    //	    snew[bi](i,j,k,2*n) = 0.0;
	    //	    snew[bi](i,j,k,2*n+1) = std::exp(-16.*rr2) * std::pow(std::cos(Pi*rr2),6);

	    // snew[bi](i,j,k,0) = std::cos(k_r*rr2);
	    // snew[bi](i,j,k,1) = omega*std::sin(k_r*rr2);


	    // snew[bi](i,j,k,0) = 1+ampl[0]*std::exp(-width[0]*rr2);
	    // snew[bi](i,j,k,1) = 0;
	    
	    snew[bi](i,j,k,0) = 4*4*4*
	      std::atan(mu_coeff*std::sin((1-mu*mu)*t0)/std::cosh(mu*(x-midpts[0])))*
	      std::atan(mu_coeff*std::sin((1-mu*mu)*t0)/std::cosh(mu*(y-midpts[1])))*
	      std::atan(mu_coeff*std::sin((1-mu*mu)*t0)/std::cosh(mu*(z-midpts[2])));
	    snew[bi](i,j,k,1) = 0;


	  }

    });
}
