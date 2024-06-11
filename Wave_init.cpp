#include "AmrLevelWave.H"
#include <cmath>

using namespace amrex;

void AmrLevelWave::initData() {
  const auto problo = geom.ProbLoArray();
  const auto probhi = geom.ProbHiArray();
  const auto dx = geom.CellSizeArray();

  MultiFab &S_new = get_new_data(State_Type);
  auto const &snew = S_new.arrays();

  amrex::Real alpha = 0.7;
  amrex::Real beta = std::sqrt(1 - alpha * alpha);
  amrex::Real t0 = -5.3;

  amrex::Real x0 = 0.5 * (probhi[0] - problo[0]);
  amrex::Real y0 = 0.5 * (probhi[1] - problo[1]);
  amrex::Real z0 = 0.5 * (probhi[2] - problo[2]);

  amrex::ParallelFor(S_new, [=] AMREX_GPU_DEVICE(int bi, int i, int j,
                                                 int k) noexcept {
    Real x = problo[0] + (i + 0.5) * dx[0];
    Real y = problo[1] + (j + 0.5) * dx[1];
    Real z = problo[2] + (k + 0.5) * dx[2];

    Real rr2 = (x - x0) * (x - x0) + (y - y0) * (y - y0) +
               (z - z0) * (z - z0); // this is the radius
    Real rr = sqrt(rr2);

    constexpr Real Pi = 3.1415926535897932384626;

    for (int n = 0; n < nfields; n++) {
      //          snew[bi](i, j, k, 2 * n) = 1.0 + ampl[n] * std::exp(-rr2 /
      //          width[n]);
      //	  snew[bi](i, j, k, 2 * n) =
      // 4*std::atan(beta*std::cos(alpha*t0)/alpha/std::cosh(beta*(x-x0)));
      snew[bi](i, j, k, 2 * n) = 4 *
                                 std::atan(alpha * std::sin(beta * t0) / beta /
                                           std::cosh(alpha * (x - x0))) *
                                 4 *
                                 std::atan(alpha * std::sin(beta * t0) / beta /
                                           std::cosh(alpha * (y - y0))) *
                                 4 *
                                 std::atan(alpha * std::sin(beta * t0) / beta /
                                           std::cosh(alpha * (z - z0)));
      snew[bi](i, j, k, 2 * n + 1) = 0.0;
    }
  });
}
