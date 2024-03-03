#include "Potential.H"

amrex::Real Potential::phi_sq(const amrex::Vector<amrex::Real> phi) const
{
  amrex::Real phi2 = 0;

  for (auto it = phi.begin(); it != phi.end(); ++it)  //for (int i : phi)
    phi2 += (*it)*(*it);   

  
  return 0.5*m_mass*m_mass*phi2;

}

amrex:: Real Potential::sine_gordon(const amrex::Real phi) const
{

  return std::sin(phi); 

}





