#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <cmath>

class Potential {
private:
  amrex::Real m_mass;

public:
  Potential(const amrex::Real mass) : m_mass(mass) {};
  ~Potential() {};

  [[nodiscard]] AMREX_INLINE AMREX_GPU_DEVICE amrex::Real
  phi_sq(const amrex::Vector<amrex::Real> phi) const {
    amrex::Real phi2 = 0;

    for (auto it = phi.begin(); it != phi.end(); ++it) // for (int i : phi)
      phi2 += (*it) * (*it);

    return 0.5 * m_mass * m_mass * phi2;
  }
};
