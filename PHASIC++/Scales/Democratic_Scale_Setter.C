#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"

namespace PHASIC {

  class Democratic_Scale_Setter: public Scale_Setter_Base {
  private:
    std::string                 m_muf2tag, m_mur2tag;
    ATOOLS::Algebra_Interpreter m_muf2calc, m_mur2calc;
    Tag_Setter                  m_tagset;
    ATOOLS::Flavour_Vector      m_f;
    double FindKT2Max();
  public:
    Democratic_Scale_Setter(const Scale_Setter_Arguments &args);

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode);
    void SetScale(const std::string &mu2tag,Tag_Setter &mu2tagset,
		  ATOOLS::Algebra_Interpreter &mu2calc);
  };// end of class Scale_Setter_Base
}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Democratic_Scale_Setter,"Democratic",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,Democratic_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new Democratic_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    Democratic_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"democratic scale scheme";
}

Democratic_Scale_Setter::Democratic_Scale_Setter(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args), m_tagset(this)
{
  size_t pos(args.m_scale.find('{'));
  std::string mur2tag("MU_R2"), muf2tag("MU_F2");
  if (pos!=std::string::npos) {
    muf2tag=args.m_scale.substr(pos+1);
    pos=muf2tag.rfind('}');
    if (pos==std::string::npos)
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    muf2tag=muf2tag.substr(0,pos);
    pos=muf2tag.find("}{");
    if (pos==std::string::npos) {
      mur2tag=muf2tag;
    }
    else {
      mur2tag=muf2tag.substr(pos+2);
      muf2tag=muf2tag.substr(0,pos);
    }
  }
  SetScale(muf2tag,m_tagset,m_muf2calc);
  SetScale(mur2tag,m_tagset,m_mur2calc);
  SetCouplings();
}

double Democratic_Scale_Setter::
Calculate(const std::vector<ATOOLS::Vec4D> &momenta,const size_t &mode) 
{
  while (!m_ampls.empty()) { m_ampls.back()->Delete(); m_ampls.pop_back(); }
  m_ampls.clear();
  m_p  = p_proc->Momenta();
  m_f  = p_proc->Flavours();
  std::vector<std::vector<int> > * cols = p_proc->Colours();
  if (cols==NULL) {
    msg_Error()<<"Error in "<<METHOD<<": "
	       <<"didn't find process-defined colours, will exit the run.\n";
    exit(0);
  }
  for (size_t i(0);i<p_proc->NIn();++i) {
    m_p[i]=-m_p[i];
    m_f[i]=m_f[i].Bar();
    int help = (*cols)[i][1];
    (*cols)[i][1] = (*cols)[i][0];
    (*cols)[i][0] = help;
  }
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  m_ampls.push_back(ampl);
  ampl->SetNIn(p_proc->NIn());
  for (size_t i=0;i<m_p.size();++i) {
    ampl->CreateLeg(m_p[i],m_f[i],ColorID((*cols)[i][0],(*cols)[i][1]));
  }  
  double kt2max = (p_proc->HasInternalScale()?
		   sqr(p_proc->InternalScale()):
		   FindKT2Max());
  ampl->SetMuQ2(Max(1.,kt2max));
  ampl->SetMuF2(Max(1.,kt2max));
  ampl->SetMuR2(Max(1.,kt2max));
  ampl->SetKT2(Max(1.,kt2max));
  ampl->SetMu2(Max(1.,kt2max));
  ampl->SetProc(p_proc);
  ampl->SetMS(p_proc->Generator());
  m_scale[stp::ren]=m_scale[stp::fac]=Max(1.,kt2max);
  msg_Debugging()<<METHOD<<" ("<<p_proc->NIn()<<" --> "<<(m_f.size()-2)<<" process,\n"
		 <<"   sqrt{kt2max} = "<<sqrt(kt2max)<<" from direct = "
		 <<sqr(p_proc->InternalScale())<<" ("<<p_proc->HasInternalScale()<<")\n";
  return m_scale[stp::fac];
}

void Democratic_Scale_Setter::SetScale
(const std::string &mu2tag,Tag_Setter &mu2tagset,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Name()<<"' {\n";
  msg_Indent();
  mu2tagset.SetCalculator(&mu2calc);
  mu2calc.SetTagReplacer(&mu2tagset);
  mu2tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  msg_Debugging()<<"}\n";
}

double Democratic_Scale_Setter::FindKT2Max() {
  std::vector<std::vector<int> > * cols = p_proc->Colours();
  double kt2max = 0., shat = (m_p[0]+m_p[1]).Abs2();
  double kt2test, kt2, DeltaRij2; 
  size_t trip_i, anti_i;
  for (size_t i=0;i<m_p.size();++i) {
    trip_i = (*cols)[i][0]; anti_i = (*cols)[i][1];
    for (size_t j=i+1;j<m_p.size();++j) {
      if (trip_i!=(*cols)[j][1] && anti_i!=(*cols)[j][0]) continue;
      kt2test = kt2 = shat;
      DeltaRij2 = -1.;
      if (i<2) {
	if (j>1) kt2 = m_p[j].MPerp2();
      }
      else {
	DeltaRij2 = cosh(m_p[i].Eta()-m_p[j].Eta())-cos(m_p[i].Phi()-m_p[j].Phi());
	kt2 = Min(m_p[i].MPerp2(),m_p[j].MPerp2())*DeltaRij2;
      }
      if (kt2<kt2test) kt2test = kt2;
    }
    if (kt2test>kt2max) kt2max = kt2test;
  }
  return kt2max;
}
