#include "AmrLevelWave.H"
#include <cmath>

using namespace amrex;


void AmrLevelWave::initData() {
  const auto problo = geom.ProbLoArray();
  const auto dx = geom.CellSizeArray();

  MultiFab &S_new = get_new_data(State_Type);
  auto const &snew = S_new.arrays();

  amrex::ParallelFor(
      S_new, [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) noexcept {
        Real x = problo[0] + (i + 0.5) * dx[0];
        Real y = problo[1] + (j + 0.5) * dx[1];
        Real z = problo[2] + (k + 0.5) * dx[2];

        Real rr2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) +
                   (z - 0.5) * (z - 0.5); // this is the radius
        Real rr = sqrt(rr2);

        constexpr Real Pi = 3.1415926535897932384626;

        for (int n = 0; n < nfields; n++) {
          snew[bi](i, j, k, 2 * n) = 1.0 + ampl[n] * std::exp(-rr2 / width[n]);
          snew[bi](i, j, k, 2 * n + 1) = 0.0;
        }
      });
}
