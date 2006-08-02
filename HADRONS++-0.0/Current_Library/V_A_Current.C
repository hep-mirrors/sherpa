#include "V_A_Current.H"
#include "XYZFuncs.H"

using namespace HADRONS;
using namespace ATOOLS;

void V_A_Current::SetModelParameters( struct GeneralModel _md )
{
  m_cR   = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL   = Complex(0.,_md("v",1.)+_md("a",1.));

  switch( int(_md("V_A_FORM_FACTOR", 1)+0.5) ) {
  case 1:
    p_ff = NULL;
    msg.Tracking()<<"Using no form factor for "<<m_name<<std::endl;
    break;
  default:
      msg.Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    abort();
  }
}


void V_A_Current::Calc()
{
  XYZFunc F(m_n-1, p_moms, p_flavs, m_k0n);
  double factor = 1.0;
  if(p_ff) {
    double q2 = (p_moms[0] - p_moms[1]).Abs2();
    factor    = p_ff->ff(q2);
  }
  for(int i=0; i<m_spin_combinations; i++) {
    // 0 is "the barred spinor" in the current, 1 is the not-barred one
    p_results[i] = factor*F.L(0,1,i,m_cR,m_cL);
  }
  F.Delete();
}

DECLARE_GETTER(V_A_Current_Getter, "V_A_Current",
               Current_Base,Flavour_Info);

Current_Base* V_A_Current_Getter::operator()(const Flavour_Info &parameters) const
{
  return new V_A_Current(parameters.flavs, parameters.nout, parameters.indices, "V_A_Current");
}

void V_A_Current_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}
