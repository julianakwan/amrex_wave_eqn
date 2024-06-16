//#include <AMReX_REAL.H>
//#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
//#include "AmrLevelWave.H"
#include <cmath>

class InitialConditions
{
 private:
  const amrex::Real m_alpha;  
  const amrex::Real m_beta;

  public:
  InitialConditions(const amrex::Real alpha, const amrex::Real k_r):m_alpha(alpha), m_beta(std::sqrt(1.0-alpha*alpha)){};

  amrex::Real breather_solution(const amrex::Real x, const amrex::Real t) const;
  amrex::Real breather_solution(const amrex::Real x, const amrex::Real y, const amrex::Real z, const amrex::Real t) const;
  amrex::Real breather_solution_deriv(const amrex::Real x, const amrex::Real t) const;


  amrex::Real travelling_wave(const amrex::Real x, const amrex::Real y, const amrex::Real z, const amrex::Real t) const;
  amrex::Real travelling_wave_deriv(const amrex::Real x, const amrex::Real y, const amrex::Real z, const amrex::Real t) const;
};
