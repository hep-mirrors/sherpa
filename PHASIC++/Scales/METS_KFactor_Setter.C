#include "PHASIC++/Scales/KFactor_Setter_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace PHASIC {

  class METS_KFactor_Setter: public KFactor_Setter_Base {
  private:

    size_t m_oqcdlo;

    double m_asref;

  public:

    METS_KFactor_Setter(Process_Base *const proc,
			const size_t &oqcdlo,const size_t &oewlo);

    double KFactor();

  };// end of class KFactor_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(METS_KFactor_Setter_Getter,"METS",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *METS_KFactor_Setter_Getter::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new METS_KFactor_Setter
    (args.p_proc,args.m_oqcdlo,args.m_oewlo);
}

void METS_KFactor_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"METS kfactor scheme\n";
}

METS_KFactor_Setter::METS_KFactor_Setter
(Process_Base *const proc,const size_t &oqcdlo,const size_t &oewlo):
  KFactor_Setter_Base(proc), m_oqcdlo(oqcdlo)
{
}

double METS_KFactor_Setter::KFactor() 
{
  if (!m_on) return 1.0;
  if (!m_kfkey.Assigned()) {
     std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),3,0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
    m_asref=MODEL::as->AlphaS(rpa.gen.CplScale());
  }
  if (m_kfkey.Weight()!=ATOOLS::UNDEFINED_WEIGHT) return m_kfkey.Weight();
  if (p_proc->OrderQCD()<0 || p_proc->OrderEW()<0) {
    THROW(fatal_error,"Couplings not set for process '"+p_proc->Name()+"'");
  }
  m_kfkey<<1.0;
  if (p_proc->OrderQCD()>0) {
    m_kfkey<<pow(MODEL::as->AlphaS(m_kfkey[0])/m_asref,p_proc->OrderQCD());
    msg_Debugging()<<METHOD<<"(): "<<p_proc->Name()<<" ("<<p_proc->NQCD()<<","
		   <<p_proc->OrderQCD()<<") {\n"
		   <<"  \\mu_{fac}   = "<<sqrt(m_kfkey[1])<<"\n"
		   <<"  \\mu_{ren}   = "<<sqrt(m_kfkey[0])<<"\n"
		   <<"} -> as = "<<MODEL::as->AlphaS(m_kfkey[0])
		   <<" => K = "<<m_kfkey.Weight()<<"\n";
  }
  return m_kfkey.Weight();
}

