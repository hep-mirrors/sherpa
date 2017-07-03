#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"
#include "PHASIC++/Selectors/Selector.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/JadePlugin.hh"

namespace PHASIC {
  class Fastjet_Veto : public Selector_Base {
    double m_ptmin,m_etmin,m_delta_r,m_f,m_eta,m_y;
    size_t m_nj, m_nb, m_nb2, m_eekt;
    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;
    fastjet::EECambridgePlugin * p_eecamplug;
    fastjet::JadePlugin * p_jadeplug;

  public:
    Fastjet_Veto(Process_Base *const proc,std::string algo,
                 double ptmin, double etmin, double dr, double f,
                 double eta, double y, int nj, int nb, int nb2);

    ~Fastjet_Veto();


    bool   Trigger(const ATOOLS::Vec4D_Vector &p,
                   ATOOLS::NLO_subevt *const sub=NULL);

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

Fastjet_Veto::Fastjet_Veto(Process_Base *const proc,
			   std::string algo, double ptmin, double etmin,
			   double dr, double f, double eta, double y,
			   int nj, int nb, int nb2) :
  Selector_Base("Fastjetfinder",proc), m_ptmin(ptmin), m_etmin(etmin),
  m_delta_r(dr), m_f(f), m_eta(eta), m_y(y),
  m_nj(nj), m_nb(nb), m_nb2(nb2), m_eekt(0), p_jdef(0),
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

  m_smin       = Max(sqr(m_ptmin),sqr(m_etmin));
}


Fastjet_Veto::~Fastjet_Veto() {
  delete p_jdef;
  if (p_siscplug) delete p_siscplug;
  if (p_eecamplug) delete p_eecamplug;
  if (p_jadeplug) delete p_jadeplug;
}


bool Fastjet_Veto::Trigger(const Vec4D_Vector &p,
                           ATOOLS::NLO_subevt *const sub)
{
  if (m_nj<0) return false;
  size_t n(sub?sub->m_n:m_n);
  const Flavour *const fl (sub?sub->p_fl:p_fl);

  DEBUG_FUNC((p_proc?p_proc->Flavours():Flavour_Vector()));
  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<n;++i) {
    if (Flavour(kf_jet).Includes(fl[i]) ||
        ((m_nb>0 || m_nb2>0) && fl[i].Kfcode()==kf_b)) {
      input.push_back(MakePseudoJet(fl[i], p[i]));
    }
  }
  if (msg_LevelIsDebugging()) {
    for (size_t i(0);i<input.size();++i)
      msg_Out()<<input[i].user_index()<<": "
               <<"("<<input[i].E()<<","<<input[i].px()
               <<","<<input[i].py()<<","<<input[i].pz()<<")"<<std::endl;
  }

  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=fastjet::sorted_by_pt(cs.inclusive_jets());
  msg_Debugging()<<"njets(ini)="<<jets.size()<<std::endl;

  if (m_eekt) {
    size_t nj(0);
    for (size_t i(0);i<input.size();++i)
      if (cs.exclusive_dmerge_max(i)>sqr(m_ptmin)) ++nj;
    return (1-m_sel_log->Hit(1-(nj>=m_nj)));
  }

  size_t nj(0), nb(0), nb2(0);
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    msg_Debugging()<<"Jet "<<i<<": pT="<<pj.PPerp()<<", |eta|="<<dabs(pj.Eta())
                   <<", |y|="<<dabs(pj.Y())<<std::endl;
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin &&
	(m_eta==100 || dabs(pj.Eta())<m_eta) &&
	(m_y==100 || dabs(pj.Y())<m_y)) {
      nj++;
      if (BTag(jets[i], 1)) nb++;
      if (BTag(jets[i], 2)) nb2++;
    }
  }
  msg_Debugging()<<"njets(fin)="<<n<<std::endl;

  bool trigger(true);
  if (nj>m_nj)   trigger=false;
  if (nb>m_nb)   trigger=false;
  if (nb2>m_nb2) trigger=false;

  if (!trigger) msg_Debugging()<<"Point discarded by jet veto"<<std::endl;
  else          msg_Debugging()<<"Point passed"<<std::endl;
  return (1-m_sel_log->Hit(1-trigger));
}


DECLARE_ND_GETTER(Fastjet_Veto,"FastjetVeto",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Veto>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<5) THROW(critical_error,"Invalid syntax");
 
  double f(.75);
  if (key.front().size()>=6) f=ToType<double>(key[0][5]);
  double eta(100.), y(100.);
  if (key.front().size()>=7) eta=ToType<double>(key[0][6]);
  if (key.front().size()>=8) y=ToType<double>(key[0][7]);
  int nb(-1), nb2(-1);
  if (key.front().size()>=9) nb=ToType<int>(key[0][8]);
  if (key.front().size()>=10) nb2=ToType<int>(key[0][9]);
#ifndef USING__FASTJET__3
  if (nb>0 || nb2>0) THROW(fatal_error, "b-tagging needs FastJet >= 3.0.");
#endif

  Fastjet_Veto *jf(new Fastjet_Veto(key.p_proc,key[0][0],
                                    ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2])),
                                    ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3])),
                                    ToType<double>(key[0][4]),f,eta,y,
                                    ToType<int>(key[0][1]),nb,nb2));
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Veto>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"FastjetVeto algorithm n ptmin etmin dr [f(siscone)=0.75 [eta=100 [y=100 [nb=-1 [nb2=-1]]]]\n"
     <<"            algorithm: kt,antikt,cambridge,siscone";
}

#endif
