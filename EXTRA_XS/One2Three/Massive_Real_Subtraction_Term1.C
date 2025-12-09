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

  for (size_t i=0; i<size(); ++i) {
    std::cout << (*this)[i] << std::endl;
  }

  double D_gb_bb = V_gb_bb/ ((p_g + p_b)*(p_g + p_b) - m2_ij) * ME2_Born;
  
  subtraction_term = D_gb_bb;
}