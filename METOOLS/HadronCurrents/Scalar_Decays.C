#include "METOOLS/HadronCurrents/Scalar_Decays.H"
#include "METOOLS/HadronCurrents/Line_Shapes.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <set>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

////////////////////////////////////////////////////////////////////
//
// Pseudo(Scalar) -> (Pseudo)Scalar + (Pseudo)Scalar
//
////////////////////////////////////////////////////////////////////

S_PP::S_PP(const Flavour & inflav,const vector<Flavour> & outflavs,
	   const double & BR) :
  Partial_Width_Base(inflav,outflavs,BR)
{
  FixPrefactor();
}

const double S_PP::Calculate(const double & s) {
  m_lambda2 = ( sqr(s-m_decmasses2[0]-m_decmasses2[1])-
		4.*m_decmasses2[0]*m_decmasses2[1] )/s;
  return Flux(s) * ME2(s) * PS_2(s);
}

const double S_PP::ME2(const double & s)  const {
  return m_prefactor;
}

const double S_PP::PS_2(const double & s) const {
  return sqrt(m_lambda2)/(8.*M_PI*s);
}

