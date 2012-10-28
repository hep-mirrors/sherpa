#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/JadePlugin.hh"

namespace PHASIC {
  class Fastjet_Finder : public Selector_Base {
    double m_ptmin,m_etmin,m_delta_r,m_f,m_eta,m_y;
    int m_nb, m_eekt;
    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;
    fastjet::EECambridgePlugin * p_eecamplug;
    fastjet::JadePlugin * p_jadeplug;

    bool BTag(const fastjet::PseudoJet& jet);

  public:
    Fastjet_Finder(int nin, int nout,ATOOLS::Flavour * fl,std::string algo,
		   double ptmin, double etmin, double dr, double f, double eta, double y, int nn, int nb);

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
			       double ptmin, double etmin, double dr, double f, double eta, double y, int nn, int nb) :
  Selector_Base("Fastjetfinder"), m_ptmin(ptmin), m_etmin(etmin), 
  m_delta_r(dr), m_f(f), m_eta(eta), m_y(y), m_nb(nb), m_eekt(0), p_jdef(0),
  p_siscplug(NULL), p_eecamplug(NULL), p_jadeplug(NULL)
{
  bool ee(rpa->gen.Beam1().IsLepton() && rpa->gen.Beam2().IsLepton());

  fastjet::JetAlgorithm ja(fastjet::kt_algorithm);

  if (algo=="cambridge") ja=fastjet::cambridge_algorithm;
  if (algo=="antikt")    ja=fastjet::antikt_algorithm;
  if (algo=="siscone") p_siscplug=new fastjet::SISConePlugin(m_delta_r,m_f);
  if (ee) {
    if (algo=="eecambridge") p_eecamplug=new fastjet::EECambridgePlugin(dr);
    if (algo=="jade") p_jadeplug=new fastjet::JadePlugin();
  }

  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else if (p_eecamplug) p_jdef=new fastjet::JetDefinition(p_eecamplug);
  else if (p_jadeplug) p_jdef=new fastjet::JetDefinition(p_jadeplug);
  else if (ee) {
    p_jdef=new fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    m_eekt=1;
  }
  else p_jdef=new fastjet::JetDefinition(ja,m_delta_r);

  m_fl         = fl;
  m_smin       = Max(sqr(m_ptmin),sqr(m_etmin));
  m_smax       = sqr(rpa->gen.Ecms());

  m_nin        = nin;
  m_nout       = nout;
  m_n          = nn;

  m_sel_log    = new Selector_Log(m_name);
}


Fastjet_Finder::~Fastjet_Finder() {
  delete p_jdef;
  if (p_siscplug) delete p_siscplug;
  if (p_eecamplug) delete p_eecamplug;
  if (p_jadeplug) delete p_jadeplug;
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
    if (m_fl[i].Strong()) {
      fastjet::PseudoJet tmp(p[i][1],p[i][2],p[i][3],p[i][0]);
      tmp.set_user_index(m_fl[i].HepEvt());
      input.push_back(tmp);
    }
  }
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=cs.inclusive_jets();

  if (m_eekt) {
    int n(0);
    for (size_t i(0);i<input.size();++i)
      if (cs.exclusive_dmerge_max(i)>sqr(m_ptmin)) ++n;
    return (1-m_sel_log->Hit(1-(n>=m_n)));
  }

  int n(0), nb(0);
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin &&
	(m_eta==100 || dabs(pj.Eta())<m_eta) &&
	(m_y==100 || dabs(pj.Y())<m_y)) {
      n++;
      if (BTag(jets[i])) nb++;
    }
  }

  bool trigger(true);
  if (n<m_n) trigger=false;
  if (nb<m_nb) trigger=false;

  return (1-m_sel_log->Hit(1-trigger));
}

bool Fastjet_Finder::JetTrigger(const Vec4D_Vector &p,
				ATOOLS::NLO_subevtlist *const subs)
{
  if (m_n<1) return true;

  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<subs->back()->m_n;++i) {
    if (subs->back()->p_fl[i].Strong()) {
        fastjet::PseudoJet tmp(p[i][1],p[i][2],p[i][3],p[i][0]);
        tmp.set_user_index(m_fl[i].HepEvt());
        input.push_back(tmp);
    }

  }
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=cs.inclusive_jets();

  if (m_eekt) {
    int n(0);
    for (size_t i(0);i<input.size();++i)
      if (cs.exclusive_dmerge_max(i)>sqr(m_ptmin)) ++n;
    return (1-m_sel_log->Hit(1-(n>=m_n)));
  }

  int n(0), nb(0);
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin &&
	(m_eta==100 || dabs(pj.Eta())<m_eta) &&
	(m_y==100 || dabs(pj.Y())<m_y)) {
      n++;
      if (BTag(jets[i])) nb++;
    }
  }

  bool trigger(true);
  if (n<m_n) trigger=false;
  if (nb<m_nb) trigger=false;

  return (1-m_sel_log->Hit(1-trigger));
}

bool Fastjet_Finder::BTag(const fastjet::PseudoJet& jet)
{
  if (m_nb<0) return false; // for performance reasons

#ifdef USING__FASTJET__3
  int nb=0;
  std::vector<fastjet::PseudoJet> cons = jet.constituents();
  for (size_t i=0; i<cons.size(); ++i) {
    if (cons[i].user_index()==5) ++nb;
    if (cons[i].user_index()==-5) --nb;
  }
  return (nb!=0);
#else
  return false;
#endif
}


namespace PHASIC{

DECLARE_ND_GETTER(Fastjet_Finder_Getter,"FastjetFinder",Selector_Base,Selector_Key,true);

Selector_Base *Fastjet_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<5) THROW(critical_error,"Invalid syntax");
 
  double f(.75);
  if (key.front().size()>=6) f=ToType<double>(key[0][5]);
  double eta(100.), y(100.);
  if (key.front().size()>=7) eta=ToType<double>(key[0][6]);
  if (key.front().size()>=8) y=ToType<double>(key[0][7]);
  int nb(-1);
  if (key.front().size()>=9) nb=ToType<int>(key[0][8]);
#ifndef USING__FASTJET__3
  if (nb>0) THROW(fatal_error, "b-tagging in FastjetFinder needs FastJet >= 3.0.");
#endif

  Fastjet_Finder *jf(new Fastjet_Finder(key.p_proc->NIn(),key.p_proc->NOut(),
					(Flavour*)&key.p_proc->Process()->Flavours().front(),
					key[0][0],
					ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2])),
					ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3])),
					ToType<double>(key[0][4]),f,eta,y,
					ToType<int>(key[0][1]),nb));
  jf->SetProcess(key.p_proc);
  return jf;
}

void Fastjet_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"FastjetFinder algorithm n ptmin etmin dr [f(siscone)=0.75 [eta=100 [y=100 [nb=-1]]]]\n"
     <<"              algorithm: kt,antikt,cambridge,siscone";
}

}

#endif
