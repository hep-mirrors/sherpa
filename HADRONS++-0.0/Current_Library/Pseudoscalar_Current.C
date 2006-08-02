#include "Pseudoscalar_Current.H"
#include "Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;

void Pseudoscalar_Current::SetModelParameters( struct GeneralModel _md )
{
  m_pionmode = (p_flavs[0].Kfcode() == kf::pi_plus) ? 1 : 0;
  m_Vxx  = m_pionmode ? 
      _md("Vud", rpa.gen.ComplexMatrixElement(std::string("CKM"), 0, 0).real()) : 
      _md("Vus", rpa.gen.ComplexMatrixElement(std::string("CKM"), 0, 1).real());
  m_fxx  = m_pionmode ? _md("fpi", 0.0924) : _md("fK", 0.113);
}

void Pseudoscalar_Current::Calc()
{
  // 0 is Pseudoscalar
  double factor=-1.0*m_fxx*sqrt(2.0)*m_Vxx;
  p_results[0] = ComplexVec4D( Vec4D(0.0,0.0,0.0,0.0), factor*p_moms[0]);
}



DECLARE_GETTER(Pseudoscalar_Current_Getter, "Pseudoscalar_Current",
               Current_Base,Flavour_Info);

Current_Base* Pseudoscalar_Current_Getter::operator()(const Flavour_Info &parameters) const
{
  return new Pseudoscalar_Current(parameters.flavs, parameters.nout, parameters.indices, "Pseudoscalar_Current");
}

void Pseudoscalar_Current_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}
