#include "BEAM/Main/Beam_Rejection.H"

#include "ATOOLS/Org/Message.H"

#include <cstdlib>
#include <vector>

using namespace BEAM;

namespace {

  // The integrals are plain sums over a uniform grid. These converge
  // exponentially in the number of nodes per skin depth for the smooth and
  // decaying density, and with b on the same grid T_AA needs no interpolation
  // of T_A. For large b, T_A(s) T_A(|b - s|) is flat along the line between the
  // nuclei, so the grid has to reach b itself. It ends at b = 2 R + 20 d, where
  // T_AA has fallen by about exp(-20) from the point where the nuclei touch.
  constexpr int nodes_per_d = 4;
  constexpr double edge = 20.;

  int GridSize(const Woods_Saxon& density)
  {
    return static_cast<int>(
        std::ceil((2. * density.R() / density.d() + edge) * nodes_per_d));
  }

} // namespace

Nucleus_Nucleus_Rejection::Nucleus_Nucleus_Rejection(
    const ATOOLS::Flavour& flav1, const ATOOLS::Flavour& flav2,
    const Woods_Saxon& density, double sNN)
    : Beam_Rejection_Base(flav1, flav2),
      // ln^2 fit of Table I in arxiv:1710.07098, in mb with s in GeV^2
      m_sigmaNN((25.0 + 0.146 * ATOOLS::sqr(std::log(sNN))) * 1.e9 /
                ATOOLS::rpa->Picobarn()),
      m_bmax(GridSize(density) * density.d() / nodes_per_d),
      m_lnTAA(ATOOLS::axis(GridSize(density), 0., m_bmax))
{
  const int n = GridSize(density);
  const double h = density.d() / nodes_per_d;

  // T_A / A at s = h sqrt(i^2 + j^2)
  std::vector<double> table((n + 1) * (n + 1));
  const auto T = [&](int i, int j) -> double& {
    return table[std::abs(i) * (n + 1) + std::abs(j)];
  };
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j)
      T(i, j) = density.Thickness(h * std::sqrt(double(i * i + j * j)));

  // T_AA at b = h m, with the second nucleus centred at (h m, 0)
  const double A = flav1.GetMassNumber();
  for (int m = 0; m <= n; ++m) {
    double sum = 0.;
    for (int i = m - n; i <= n; ++i)
      for (int j = -n; j <= n; ++j)
        sum += T(i, j) * T(i - m, j);
    m_lnTAA.Fill(m, std::log(ATOOLS::sqr(A * h) * sum));
  }

  msg_Info() << "Optical Glauber overlap for " << flav1
             << ": Woods-Saxon R = " << density.R() * ATOOLS::rpa->hBar_c()
             << " fm, d = " << density.d() * ATOOLS::rpa->hBar_c()
             << " fm, sigma_NN = " << m_sigmaNN * ATOOLS::rpa->Picobarn() / 1.e9
             << " mb.\n";
}
