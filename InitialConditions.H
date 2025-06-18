#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <cmath>

class InitialConditions {
private:
  const amrex::Real m_alpha;
  const amrex::Real m_beta;

public:
  InitialConditions(const amrex::Real alpha, const amrex::Real k_r)
      : m_alpha(alpha), m_beta(std::sqrt(1.0 - alpha * alpha)) {};

  [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
  breather_solution(const amrex::Real x, const amrex::Real t) const {
    // Sine Gordon 1D breather solution
    amrex::Real x1 = 0;
    amrex::Real x2 = 0;
    amrex::Real v = 0;

    amrex::Real y1 = t - v * x + x1;
    amrex::Real y2 = x - v * t + x2;

    return 4 * std::atan(m_beta * std::cos(m_alpha * t) / m_alpha /
                         std::cosh(m_beta * x));
  }

  [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
  breather_solution(const amrex::Real x, const amrex::Real y,
                    const amrex::Real z, const amrex::Real t) const {
    // Sine Gordon 3D psudo-breather solution from arXiv:1212.2716

    return 4 * 4 * 4 *
           std::atan(m_alpha * std::sin(m_beta * t) / m_beta /
                     std::cosh(m_alpha * x)) *
           std::atan(m_alpha * std::sin(m_beta * t) / m_beta /
                     std::cosh(m_alpha * y)) *
           std::atan(m_alpha * std::sin(m_beta * t) / m_beta /
                     std::cosh(m_alpha * z));
  }

  [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
  breather_solution_deriv(const amrex::Real x, const amrex::Real t) const {
    // First derivative of Sine Gordon 1D breather solution
    amrex::Real x1 = 0;
    amrex::Real x2 = 0;
    amrex::Real v = 0;

    amrex::Real y1 = t - v * x + x1;
    amrex::Real y2 = x - v * t + x2;

    amrex::Real numerator =
        m_alpha * std::sin(m_alpha * y1) * std::cosh(m_beta * y2);
    amrex::Real denominator =
        m_alpha * m_alpha * std::cosh(m_beta * y2) * std::cosh(m_beta * y2) +
        m_beta * m_beta * std::cos(m_alpha * y1) * std::cos(m_alpha * y1);

    return -4 * m_alpha * m_beta * numerator / denominator;
  }

  [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
  travelling_wave(const amrex::Real x, const amrex::Real y, const amrex::Real z,
                  const amrex::Real t) const {
    // for the wave to be at the center of the grid, need to pass in
    // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
    amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

    return std::cos(m_alpha * rr2 - m_alpha * t);
  }

  [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
  travelling_wave_deriv(const amrex::Real x, const amrex::Real y,
                        const amrex::Real z, const amrex::Real t) const {
    amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

    return m_alpha * std::sin(m_alpha * rr2 - m_alpha * t);
  }
};
