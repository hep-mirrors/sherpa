#include "EXTRA_XS/One2Three/Massive_Real_Subtraction_Term2.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Vec4.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;

void Massive_Real_Subtraction_Term2::Calculate_real_subtraction(const ATOOLS::Vec4D_Vector& momenta,
                                                            bool anti)
{
  // implementation based on the formulas in the Catani Dittmaier Seymour Trocsanyi paper from 2002

  // first: some variables needed later (nomination also based on paper)
  ATOOLS::Vec4<double> p_g = momenta[1];
  ATOOLS::Vec4<double> p_b = momenta[2];
  ATOOLS::Vec4<double> p_bb = momenta[3];
  double m2_ij = p_b * p_b; // because m_i = 0 (gluon) and m_b = m_bb

  double V_gbb_b = V_ijk(p_g, p_bb, p_b, m_prop);

  double prefactor = V_gbb_b/ ((p_g + p_bb)*(p_g + p_bb) - m2_ij);
  double D_gbb_b = prefactor * ME2_Born;
  
  subtraction_term = D_gbb_b;
  
  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] = 0;
  }

  // filling this spin-amplitudes object with the correct values:
  (*this)[0] = born_hel["00"]; // first 4 values: gluon_hel = 0
  (*this)[1] = born_hel["10"];
  (*this)[2] = born_hel["01"];
  (*this)[3] = born_hel["11"];
  (*this)[4] = born_hel["00"]; // last 4 values: gluon_hel = 1
  (*this)[5] = born_hel["10"];
  (*this)[6] = born_hel["01"];
  (*this)[7] = born_hel["11"];
  // all values are taken 2 times, because the spin_Amplitude of S needs to have the same size as R. R has additional gluon degrees of freedom, therefore,
  // these additional degrees of freedom are needed here as well. Divide all entries by 2 in order to correct the double counting.
  if(prefactor < 0.0){
    sign = -1.0;
    prefactor *= -1.0;
  } else sign = 1.0;
  double factor(0.5 * std::sqrt(prefactor * 3));

  for (size_t i = 0; i < this->size(); ++i) {
    (*this)[i] *= factor;   // * 3 for born colour factor * 0.5 (0.5 to correct double counting)
  }
}


double Massive_Real_Subtraction_Term2::getSign(){
  return (-1.0) * sign;    // -1.0 because the S term will be subtracted
}