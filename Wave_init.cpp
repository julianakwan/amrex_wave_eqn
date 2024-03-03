#include "AmrLevelWave.H"
#include <cmath>
#include "InitialConditions.H"

using namespace amrex;


// void
// AmrLevelWave::initData ()
// {
//     const auto problo = geom.ProbLoArray();
//     const auto dx = geom.CellSizeArray();

//     MultiFab& S_new = get_new_data(State_Type);
//     auto const& snew = S_new.arrays();

//     const amrex::Real alpha = 0.7;
    
//     SineGordon SG_breather(alpha);
    
//     SG_breather.init_data_1D(geom, S_new);

// }


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

    
    MultiFab& S_new = get_new_data(State_Type);
    auto const& snew = S_new.arrays();

    constexpr Real t0 = 0;
    InitialConditions SineGordon(alpha, k_r);    

    amrex::ParallelFor(S_new,
    [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
    {
        Real x = problo[0] + (i+0.5)*dx[0];
	Real y = problo[1] + (j+0.5)*dx[1];
	Real z = problo[2] + (k+0.5)*dx[2];


	Real rr2 = (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5);
	
	//        constexpr Real Pi = 3.1415926535897932384626;
	// Real mu = 0.7;
	// Real mu_coeff = mu/(std::sqrt(1-mu*mu));
	// constexpr Real t0 = 0;



		// constexpr Real k_r = 1;
		// constexpr Real omega = 1;
	for (int n = 0; n < nfields; n++)
	  {
	    //	    snew[bi](i,j,k,2*n) = 0.0;
	    //	    snew[bi](i,j,k,2*n+1) = std::exp(-16.*rr2) * std::pow(std::cos(Pi*rr2),6);



	    snew[bi](i,j,k,2*n) = 1+ampl[n]*std::exp(-width[n]*rr2);
	    snew[bi](i,j,k,2*n+1) = 0;
	    
	    // snew[bi](i,j,k,0) = 4*4*4*
	    //   std::atan(mu_coeff*std::sin((1-mu*mu)*t0)/std::cosh(mu*(x-midpts[0])))*
	    //   std::atan(mu_coeff*std::sin((1-mu*mu)*t0)/std::cosh(mu*(y-midpts[1])))*
	    //   std::atan(mu_coeff*std::sin((1-mu*mu)*t0)/std::cosh(mu*(z-midpts[2])));
	    // snew[bi](i,j,k,1) = 0;

	    // snew[bi](i,j,k,2*n) = SineGordon.breather_solution(x-midpts[0], t0);
	    // snew[bi](i,j,k,2*n+1) = SineGordon.breather_solution_deriv(x-midpts[0], t0);


	  }

    });
}
