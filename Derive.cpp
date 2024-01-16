#include <AMReX_REAL.H>
#include <AmrLevelWave.H>



void derive_func_fab(const amrex::Box & bx, amrex::FArrayBox& derfab,int dcomp, int /*numcomp*/,
		     const amrex::FArrayBox& datfab, const amrex::Geometry& geom, 
		     const amrex::Real time, const int* /*bcomp*/, int /*scomp*/)

{
  constexpr amrex::Real k_r = 1;  //TODO: these need to be passed in somehow...
  constexpr amrex::Real omega = 1;

  const auto problo = geom.ProbLoArray();
  const auto dx = geom.CellSizeArray();

  auto const& s = datfab.array();
  auto const& s_out = derfab.array();

  amrex::ParallelFor(bx, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		     {


		       amrex::Real x = problo[0] + (i+0.5)*dx[0];
		       amrex::Real y = problo[1] + (j+0.5)*dx[1]; 
		       amrex::Real z = problo[2] + (k+0.5)*dx[2];
					
		       amrex::Real rr2 = (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) + (z - 0.5)*(z - 0.5);  // this is the radius 

		       amrex::Real exact_soln = std::cos(k_r*rr2-omega*time);
		       if (time > 0)
			 s_out(i,j,k,dcomp) = amrex::Math::abs(s(i,j,k,0)-exact_soln);
		       else
			 s_out(i,j,k,dcomp) = 0;
		     });


}
