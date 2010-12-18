#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

namespace PHASIC {
  class Fastjet_Finder : public Selector_Base {
    double m_ptmin,m_etmin,m_delta_r,m_f;
    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;

  public:
    Fastjet_Finder(int nin, int nout,ATOOLS::Flavour * fl,std::string algo,
		   double ptmin, double etmin, double dr, double f, int nn);

    ~Fastjet_Finder();


    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    bool   Trigger(const ATOOLS::Vec4D_Vector &);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,
		      ATOOLS::NLO_subevtlist *const subs);

    void   BuildCuts(Cut_Data *) {}
  };
}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"


using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Fastjet_Finder::Fastjet_Finder(int nin, int nout,ATOOLS::Flavour * fl, std::string algo,
			       double ptmin, double etmin, double dr, double f, int nn) : 
  Selector_Base("Fastjetfinder"), m_ptmin(ptmin), m_etmin(etmin), 
  m_delta_r(dr), m_f(f), p_jdef(0), p_siscplug(0)
{
  fastjet::JetAlgorithm ja(fastjet::kt_algorithm);

  if (algo=="cambridge") ja=fastjet::cambridge_algorithm;
  if (algo=="antikt")    ja=fastjet::antikt_algorithm;
  if (algo=="siscone") p_siscplug=new fastjet::SISConePlugin(m_delta_r,m_f);

  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else p_jdef=new fastjet::JetDefinition(ja,m_delta_r);

  m_fl         = fl;
  m_smin       = Max(sqr(m_ptmin),sqr(m_etmin));
  m_smax       = sqr(rpa.gen.Ecms());

  m_nin        = nin;
  m_nout       = nout;
  m_n          = nn;

  m_sel_log    = new Selector_Log(m_name);
}


Fastjet_Finder::~Fastjet_Finder() {
  delete p_jdef;
  if (p_siscplug) delete p_siscplug;
}


bool Fastjet_Finder::NoJetTrigger(const Vec4D_Vector &p)
{
  if (m_n<1) return true;

  double s=(p[0]+p[1]).Abs2();
  return (s>m_smin*4.);
}

bool Fastjet_Finder::Trigger(const Vec4D_Vector &p)
{
  if (m_n<1) return true;

  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<p.size();++i) {
    if (m_fl[i].Strong())
      input.push_back(fastjet::PseudoJet(p[i][1],p[i][2],p[i][3],p[i][0])); 
  }
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=cs.inclusive_jets();

  int n=0;
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin) n++;
  }

  bool trigger(true);
  if (n<m_n) trigger=false;

  return (1-m_sel_log->Hit(1-trigger));
}

bool Fastjet_Finder::JetTrigger(const Vec4D_Vector &p,
				ATOOLS::NLO_subevtlist *const subs)
{
  if (m_n<1) return true;

  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<subs->back()->m_n;++i) {
    if (subs->back()->p_fl[i].Strong())
      input.push_back(fastjet::PseudoJet(p[i][1],p[i][2],p[i][3],p[i][0]));      
  }
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=cs.inclusive_jets();

  int n=0;
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin) n++;
  }

  bool trigger(true);
  if (n<m_n) trigger=false;
  
  return (1-m_sel_log->Hit(1-trigger));
}


namespace PHASIC{

DECLARE_ND_GETTER(Fastjet_Finder_Getter,"FastjetFinder",Selector_Base,Selector_Key,true);

Selector_Base *Fastjet_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<5) THROW(critical_error,"Invalid syntax");
 
  double f(.75);
  if (key.front().size()>=6) f=ToType<double>(key[0][5]);

  Fastjet_Finder *jf(new Fastjet_Finder(key.p_proc->NIn(),key.p_proc->NOut(),
					(Flavour*)&key.p_proc->Process()->Flavours().front(),
					key[0][0],
					ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2])),
					ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3])),
					ToType<double>(key[0][4]),f,
					ToType<int>(key[0][1])));
  jf->SetProcess(key.p_proc);
  return jf;
}

void Fastjet_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"FastjetFinder algorithm n ptmin etmin dr [f(siscone)=0.75]\n" 
     <<"              algorithm: kt,antikt,cambridge,siscone";
}

}

#endif
