#include "HADRONS++/Current_Library/VA_F_F.H"
#include "METOOLS/Main/XYZFuncs.H"

using namespace std;
using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;

VA_F_F::FF_Base::~FF_Base()
{
}

void VA_F_F::SetModelParameters( struct GeneralModel _md )
{
  m_cR   = Complex(0.,_md("v",1.)+_md("a",-1.));
  m_cL   = Complex(0.,_md("v",1.)-_md("a",-1.));

  switch( int(_md("V_A_FORM_FACTOR", 1)+0.5) ) {
  case 1:
    p_ff = NULL;
    msg_Tracking()<<"Using no form factor for "<<m_name<<std::endl;
    break;
  default:
      msg_Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    abort();
  }
}


void VA_F_F::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  XYZFunc F(moms, m_flavs, m_anti, p_i);
  double factor = 1.0;
  if(p_ff) {
    double q2 = (moms[p_i[0]] - moms[p_i[1]]).Abs2();
    factor    = p_ff->ff(q2);
  }
  Vec4C amp;
  for(int h0=0; h0<2; h0++) {
    for(int h1=0; h1<2; h1++) {
      // the first current index in the decay channel file has to be the
      // "barred spinor" (for the !m_anti case), the second the non-barred one
      amp=factor*F.L(0,h0, 1,h1, m_cR,m_cL);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert( amp ,spins );
    }
  }
}

DEFINE_CURRENT_GETTER(VA_F_F,"VA_F_F")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_F_F>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ (B) \\rightarrow (D) l \\nu_l $ \n\n"
    <<"Order: 0 = bar'ed spinor, 1 = non-bar'ed spinor \n\n"
    <<"\\[\\bar{u}(p_0) \\gamma_\\mu [ v-a\\gamma_5 ] u(p_1) \\] \n\n"
    <<std::endl;
}
