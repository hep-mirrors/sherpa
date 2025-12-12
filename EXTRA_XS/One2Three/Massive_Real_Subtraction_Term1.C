#include "EXTRA_XS/One2Three/Massive_Real_Subtraction_Term1.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Vec4.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;

void Massive_Real_Subtraction_Term1::Calculate_real_subtraction(const ATOOLS::Vec4D_Vector& momenta,
                                                            bool anti)
{
  // implementation based on the formulas in the Catani Dittmaier Seymour Trocsanyi paper from 2002

  // first: some variables needed later (nomination also based on paper)
  ATOOLS::Vec4<double> p_g = momenta[1];
  ATOOLS::Vec4<double> p_b = momenta[2];
  ATOOLS::Vec4<double> p_bb = momenta[3];

  double V_gb_bb = V_ijk(p_g, p_b, p_bb, m_prop);

  double m2_ij = p_b * p_b; // because m_i = 0 (gluon) and m_b = m_bb
  double prefactor = V_gb_bb/ ((p_g + p_b)*(p_g + p_b) - m2_ij);

  double D_gb_bb = prefactor * ME2_Born;
  subtraction_term = D_gb_bb;

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = 0;
  }

  // filling this spin-amplitudes object with the correct values:
  (*this)[0] = born_hel["00"];
  (*this)[1] = born_hel["10"];
  (*this)[2] = born_hel["01"];
  (*this)[3] = born_hel["11"];

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] *= std::sqrt(prefactor);
  }
}