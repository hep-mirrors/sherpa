#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Selector.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/JadePlugin.hh"

namespace PHASIC {
  class Jet_Selector : public Selector_Base {
    std::string m_algo;
    double m_ptmin,m_R,m_f,m_etamin,m_etamax,m_ymin,m_ymax;
    size_t m_nmin, m_nmax;
    kf_code m_outjetkf;

    ATOOLS::Jet_Inputs          m_jetinput;
    ATOOLS::Jet_Identifications m_jetids;

    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;

    public:
    Jet_Selector(const Selector_Key &key);
    ~Jet_Selector();

    bool   Trigger(ATOOLS::Selector_List &);

    void   BuildCuts(Cut_Data *);
  };
}



using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Jet_Selector::Jet_Selector(const Selector_Key &key) :
  Selector_Base("Jet_Selector",key.p_proc),
  m_algo(""), m_ptmin(0.), m_R(0.), m_f(0.7),
  m_etamin(-std::numeric_limits<double>::max()),
  m_etamax(std::numeric_limits<double>::max()),
  m_ymin(-std::numeric_limits<double>::max()),
  m_ymax(std::numeric_limits<double>::max()),
  m_nmin(0), m_nmax(std::numeric_limits<size_t>::max()),
  m_outjetkf(kf_none), p_jdef(NULL), p_siscplug(NULL)
{
  DEBUG_FUNC("");
  for (size_t k=0;k<key.size();++k) {
    if      (key[k].size()>1 && key[k][0]=="Input_Particles") {
      if (key[k][1]=="all") {
        THROW(not_implemented,"Keyword 'all' not implemented yet.");
      }
      else {
        int kfc(ToType<long int>(key[k][1]));
        Flavour fl(abs(kfc),kfc<0);
        m_jetinput.push_back(Jet_Input(fl));
        for (size_t i(2);i<key[k].size();++i) {
          std::string ds(key[k][i]);
          if (ds.find("=")!=std::string::npos) ds.replace(ds.find("="),1,":");
          std::string var(ds,0,ds.find(':'));
          std::string val(ds,ds.find(':')+1);
          if (!var.size() || !val.size()) THROW(fatal_error,"Input error.");
          if      (var=="PT")  m_jetinput.back().m_ptmin=ToType<double>(val);
          else if (var=="ETA") m_jetinput.back().m_etamax=ToType<double>(val);
          else if (var=="Y")   m_jetinput.back().m_ymax=ToType<double>(val);
        }
      }
    }
    else if (key[k].size()>1 && key[k][0]=="Output_ID") {
      m_outjetkf = ToType<kf_code>(key[k][1]);
    }
    else if (key[k].size()>3 && key[k][0]=="Jet_Algorithm") {
      m_algo=key[k][1];
      for (size_t i(2);i<key[k].size();++i) {
        std::string ds(key[k][i]);
        if (ds.find("=")!=std::string::npos) ds.replace(ds.find("="),1,":");
        std::string var(ds,0,ds.find(':'));
        std::string val(ds,ds.find(':')+1);
        if (!var.size() || !val.size()) THROW(fatal_error,"Input error.");
        if      (var=="PT")     m_ptmin=ToType<double>(val);
        else if (var=="R")      m_R=ToType<double>(val);
        else if (var=="f")      m_f=ToType<double>(val);
        else if (var=="ETA")    m_etamin=-(m_etamax=ToType<double>(val));
        else if (var=="ETAMIN") m_etamin=ToType<double>(val);
        else if (var=="ETAMAX") m_etamax=ToType<double>(val);
        else if (var=="Y")      m_ymin=-(m_ymax=ToType<double>(val));
        else if (var=="YMIN")   m_ymin=ToType<double>(val);
        else if (var=="YMAX")   m_ymax=ToType<double>(val);
      }
    }
    else if (key[k].size()>2 && key[k][0]=="Identify_As") {
      int kfc(ToType<long int>(key[k][1]));
      Flavour fl(abs(kfc),kfc<0);
      std::string input(key[k][2]), rel("");
      if      (input.find('>')!=std::string::npos) rel=">";
      else if (input.find('<')!=std::string::npos) rel="<";
      else THROW(not_implemented,"Unknown relation.");
      std::string var(input.substr(0,input.find(rel)));
      input=input.substr(input.find(rel)+1);
      std::string val(input.substr(0,input.find('[')));
      input=input.substr(input.find('[')+1);
      std::string mode(input.substr(0,input.find(']')));
      double ptmin(0.),etmin(0.),emin(0.);
      if      (var=="PT") ptmin=ToType<double>(val);
      else if (var=="ET") etmin=ToType<double>(val);
      else if (var=="E")  emin=ToType<double>(val);
      else THROW(not_implemented,"Unknown variable.");
      JetIdMode::code jetidmode(JetIdMode::unknown);
      if (rel==">") jetidmode|=JetIdMode::larger;
      if (mode=="rel") jetidmode|=JetIdMode::relative;
      m_jetids.push_back(new Jet_Identification(fl,ptmin,etmin,emin,jetidmode));
    }
    else if (key[k].size()>1 && key[k][0]=="NMin") {
      m_nmin = ToType<int>(key[k][1]);
    }
    else if (key[k].size()>1 && key[k][0]=="NMax") {
      m_nmax = ToType<int>(key[k][1]);
    }
    else {
      msg_Debugging()<<"Read in subselectors."<<std::endl;
      ReadInSubSelectors(key,k);
      break;
    }
  }

  m_smin       = sqr(m_ptmin);

  if (m_nmin>m_nmax) THROW(fatal_error,"Inconsistent setup.");
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"Jet Algorithm: "<<m_algo
             <<", pT>"<<m_ptmin<<", R="<<m_R<<", f="<<m_f
             <<", "<<m_etamin<<"<eta<"<<m_etamax
             <<", "<<m_ymin<<"<y<"<<m_ymax<<std::endl;
    msg_Out()<<"Jet_Input: "<<std::endl;
    for (size_t i(0);i<m_jetinput.size();++i) {
      msg_Out()<<"  "<<m_jetinput[i]<<std::endl;
    }
    msg_Out()<<"Jet_Identification: "<<m_jetids.size()<<"\n";
    for (size_t i(0);i<m_jetids.size();++i) {
      msg_Out()<<"  "<<*m_jetids[i]<<std::endl;
    }
    msg_Out()<<"NMin: "<<m_nmin<<", NMax: "<<m_nmax<<std::endl;
    msg_Out()<<"Additional Selectors: "<<m_sels.size()<<"\n";
    for (size_t i(0);i<m_sels.size();++i)
      msg_Out()<<"  "<<m_sels[i]->Name()<<std::endl;
  }

  // init jet algo
  fastjet::JetAlgorithm ja;
  if (m_algo=="kt")             ja=fastjet::kt_algorithm;
  else if (m_algo=="cambridge") ja=fastjet::cambridge_algorithm;
  else if (m_algo=="antikt")    ja=fastjet::antikt_algorithm;
  else if (m_algo=="siscone")   p_siscplug=new fastjet::SISConePlugin(m_R,m_f);
  else THROW(not_implemented,"Unknown algorithm.");

  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else            p_jdef=new fastjet::JetDefinition(ja,m_R);
}


Jet_Selector::~Jet_Selector() {
  while (m_jetids.size()>0) {
    delete *m_jetids.begin();
    m_jetids.erase(m_jetids.begin());
  }
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}


bool Jet_Selector::Trigger(Selector_List &sl)
{
  size_t n(m_n);

  DEBUG_FUNC((p_proc?p_proc->Flavours():Flavour_Vector()));
  Vec4D_Vector moms(sl.size(),Vec4D(0.,0.,0.,0.));
  for (size_t i(0);i<m_nin;++i) moms[i]=sl[i].Momentum();
  std::vector<fastjet::PseudoJet> input,jets;
  // book-keep where jet input was taken from
  std::vector<size_t> jetinputidx;
  for (size_t i(m_nin);i<n;++i) if (sl[i].Momentum()!=moms[i]) {
    if (ToBeClustered(sl[i].Flavour(), sl[i].Momentum(), m_jetinput)) {
      input.push_back(MakePseudoJet(sl[i].Flavour(),sl[i].Momentum()));
      jetinputidx.push_back(i);
    }
    else moms[i]=sl[i].Momentum();
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

  size_t njets(0);
  std::vector<std::pair<Flavour, Vec4D> > clusidjetmoms;
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pi(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    double pipt(pi.PPerp()),pieta(pi.Eta()),piy(pi.Y());
    msg_Debugging()<<"Jet "<<i<<": pT="<<pipt<<", eta="<<pieta
                   <<", y="<<piy<<std::endl;
    if (pipt>m_ptmin && pieta>m_etamin && pieta<m_etamax
                     && piy>m_ymin     && piy<m_ymax) {
      Flavour jetidflav(FlavourTag(jets[i], m_jetids, m_outjetkf));
      if (jetidflav.Kfcode()==m_outjetkf) ++njets;
      clusidjetmoms.push_back(std::make_pair(jetidflav,pi));
    }
  }
  msg_Debugging()<<"njets(fin)="<<njets<<std::endl;
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"Clustered objects with assigned flavours:\n";
    for (size_t i(0);i<clusidjetmoms.size();++i)
      msg_Out()<<clusidjetmoms[i].first<<": "<<clusidjetmoms[i].second<<std::endl;
  }

  // fill in identified jets into respective places,
  // unidentified jets have m_outjetkf, if not kf_none, also fill them in
  // remove all filled idx from jetinputidx vector
  for (std::vector<std::pair<Flavour,Vec4D> >::iterator cjit(clusidjetmoms.begin());
       cjit!=clusidjetmoms.end();) {
    bool assigned(false);
    if (cjit->first.Kfcode()!=kf_none) {
      for (size_t j(0);j<moms.size();++j) {
        if (moms[j]==Vec4D(0.,0.,0.,0.) && cjit->first.Includes(sl[j].Flavour())) {
          moms[j]=cjit->second;
          jetinputidx.erase(std::remove(jetinputidx.begin(),jetinputidx.end(),j),
                            jetinputidx.end());
          clusidjetmoms.erase(cjit);
          assigned=true;
          break;
        }
      }
    }
    if (!assigned) ++cjit;
  }

  // if unidentified jets have kf_none, fill into remaining slots in jetinputidx
  if (clusidjetmoms.size()>jetinputidx.size())
    THROW(fatal_error,"Too many jets left.");
  for (size_t i(0);i<clusidjetmoms.size();++i)
    moms[jetinputidx[i]]=clusidjetmoms[i].second;

  // remaining spots remain empty
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"Final momenta list:\n";
    for (size_t i(0);i<moms.size();++i) msg_Out()<<sl[i]<<std::endl;
  }

  bool trigger(njets>=m_nmin && njets<=m_nmax);
  if (!trigger) {
    msg_Debugging()<<"Point discarded by jet finder"<<std::endl;
    m_sel_log->Hit(true);
    return false;
  }
  for (size_t k=0;k<m_sels.size();++k) {
    if (!m_sels[k]->Trigger(sl)) {
      msg_Debugging()<<"Point discarded by subselector"<<std::endl;
      m_sel_log->Hit(true);
      return false;
    }
  }
  msg_Debugging()<<"Point passed"<<std::endl;
  m_sel_log->Hit(false);
  return true;
}

void Jet_Selector::BuildCuts(Cut_Data * cuts)
{
  for (size_t i(0);i<m_sels.size();++i) m_sels[i]->BuildCuts(cuts);
}

DECLARE_ND_GETTER(Jet_Selector,"Jet_Selector",
                  Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,
                              Jet_Selector>::operator()
(const Selector_Key &key) const
{
  return new Jet_Selector(key);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Jet_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  std::string w(width+4,' ');
  str<<"Jet_Selector {\n"
     <<w<<"  Input_Particles <kf1> [PT:<ptmin>] [ETA:<etamax>] [Y:<ymax>]\n"
     <<w<<"  [Input_Particles <kf2> [PT:<ptmin>] [ETA:<etamax>] [Y:<ymax>]]\n"
     <<w<<"  [...]\n"
     <<w<<"  Jet_Algorithm <algorithm> PT:<ptmin> R:<dR> [ETA:<etamax>] [Y:<ymax>]\n"
     <<w<<"  Indentify_As <kf> [E/ET/PT><emin>[rel/abs]]\n"
     <<w<<"  NMin <nmin>\n"
     <<w<<"  [NMax <nmax>]\n"
     <<w<<"  [Output_ID <kf>]\n"
     <<w<<"  Selector 1\n"
     <<w<<"  Selector 2\n"
     <<w<<"}";
}

#else

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace PHASIC {
  class Jet_Selector : public Selector_Base {
  public:
    Jet_Selector(const Selector_Key &key) :
      Selector_Base("Jet_Selector",key.p_proc) {}

    bool   Trigger(const ATOOLS::Vec4D_Vector &) { return true; }
    void   BuildCuts(Cut_Data *) {}
  };
}

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_ND_GETTER(Jet_Selector,"Jet_Selector",
                  Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,
                              Jet_Selector>::operator()
(const Selector_Key &key) const
{
  THROW(fatal_error, "To use the Jet_Selector Sherpa must be linked "+
                     std::string("to FastJet during compilation."));
  return new Jet_Selector(key);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Jet_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  std::string w(width+4,' ');
  str<<"Jet_Selector {\n"
     <<w<<"  To use this selector Sherpa must be compiled with FastJet support.\n"
     <<w<<"}";
}

#endif
