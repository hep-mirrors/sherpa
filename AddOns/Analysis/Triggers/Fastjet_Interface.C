#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/SISConePlugin.hh"

using namespace ANALYSIS;
using namespace ATOOLS;

class Fastjet_Interface: public Trigger_Base {
private:

  fastjet::JetDefinition m_jdef;

  fastjet::JetDefinition::Plugin *p_plug;

  size_t m_njets, m_btag;
  double m_ptmin, m_etamax;

public:

  // constructor
  Fastjet_Interface(const std::string &inlist,
		    const std::string &outlist,
		    const fastjet::JetDefinition &jdef,
 		    fastjet::JetDefinition::Plugin *const plug,
		    const size_t &njets,const double &ptmin,
		    const double &etamax,const size_t btag):
    Trigger_Base(inlist,outlist), m_jdef(jdef), p_plug(plug),
    m_njets(njets), m_btag(btag), m_ptmin(ptmin), m_etamax(etamax) {}

  ~Fastjet_Interface()
  {
    if (p_plug) delete p_plug;
  }

  // member functions
  Analysis_Object *GetCopy() const 
  {
    return new Fastjet_Interface
      (m_inlist,m_outlist,m_jdef,NULL,m_njets,m_ptmin,m_etamax,m_btag);
  }

  int BTag(ATOOLS::Particle *const p)
  {
    msg_Indent();
    if (p->ProductionBlob()==NULL ||
	p->ProductionBlob()->NInP()!=1 ||
	p->ProductionBlob()->Type()==btp::Beam) {
      if (p->Flav().IsB_Hadron() ||
	  p->Flav().Kfcode()==5) return p->Flav().IsAnti()?-5:5;
      return 0;
    }
    return BTag(p->ProductionBlob()->InParticle(0));
  }

  void Evaluate(const ATOOLS::Particle_List &plist,
		ATOOLS::Particle_List &outlist,
		double value,double ncount)
  {
    std::vector<fastjet::PseudoJet> input(plist.size()), jets;
    for (size_t i(0);i<input.size();++i) {
      Vec4D p(plist[i]->Momentum());
      input[i]=fastjet::PseudoJet(p[1],p[2],p[3],p[0]);
      input[i].set_user_index(BTag(plist[i]));
    }
    fastjet::ClusterSequence cs(input,m_jdef);
    if (m_njets>0) {
      jets=fastjet::sorted_by_pt(cs.exclusive_jets((int)m_njets));
    }
    else {
      jets=fastjet::sorted_by_pt(cs.inclusive_jets());
    }
    std::vector<double> *ktdrs(new std::vector<double>());
    for (size_t i(input.size());i>0;--i) {
      if      (m_jdef.jet_algorithm()==fastjet::kt_algorithm)
        ktdrs->push_back(cs.exclusive_dmerge(i-1));
      else if (m_jdef.jet_algorithm()==fastjet::antikt_algorithm)
        ktdrs->insert(ktdrs->begin(),1./cs.exclusive_dmerge(i-1));
    }
    std::string key("KtJetrates(1)"+m_outlist);
    p_ana->AddData(key,new Blob_Data<std::vector<double> *>(ktdrs));
    for (size_t i(0);i<jets.size();++i) {
      kf_code flav(kf_jet);
      if (m_btag) {
#ifdef USING__FASTJET__3
	int nb(0);
	const std::vector<fastjet::PseudoJet>
	  &cons(jets[i].constituents());
	for (size_t j=0;j<cons.size();++j) {
	  if (cons[j].user_index()==5) ++nb;
	  if (cons[j].user_index()==-5) --nb;
	}
	if (nb!=0) flav=kf_bjet;
#else
	THROW(fatal_error,"FastJet >= v3 required for b tags");
#endif
      }
      Vec4D jetmom(jets[i][3],jets[i][0],jets[i][1],jets[i][2]);
      if (jetmom.PPerp()>m_ptmin && abs(jetmom.Eta())<m_etamax)
        outlist.push_back(new Particle (1,Flavour(flav),jetmom));
    }
    std::sort(outlist.begin(),outlist.end(),Order_PT());
  }

};

DECLARE_GETTER(Fastjet_Interface,"FastJets",
	       Analysis_Object,Analysis_Key);

Analysis_Object *ATOOLS::Getter
<Analysis_Object,Analysis_Key,Fastjet_Interface>::
operator()(const Analysis_Key& key) const
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto inlist = s["InList"].SetDefault("FinalState").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("FastJets").Get<std::string>();
  const auto njets = s["NJets"].SetDefault(0).Get<size_t>();
  const auto ptmin = s["PTMin"].SetDefault(0.0).Get<double>();
  const auto etamax = s["EtaMax"].SetDefault(1000.0).Get<double>();

  fastjet::JetAlgorithm algo;
  size_t siscone = 0;
  const auto rawalgorithm = s["Algorithm"].SetDefault("kt").Get<std::string>();
  if (rawalgorithm=="kt") algo=fastjet::kt_algorithm;
  else if (rawalgorithm=="cambridge") algo=fastjet::cambridge_algorithm;
  else if (rawalgorithm=="antikt") algo=fastjet::antikt_algorithm;
  else if (rawalgorithm=="siscone") siscone=1;
  else THROW(fatal_error, "Unknown jet algorithm.");

  fastjet::RecombinationScheme recom;
  const auto rawscheme = s["Scheme"].SetDefault("E").Get<std::string>();
  if (rawscheme=="E") recom=fastjet::E_scheme;
  else if (rawscheme=="pt") recom=fastjet::pt_scheme;
  else if (rawscheme=="pt2") recom=fastjet::pt2_scheme;
  else if (rawscheme=="Et") recom=fastjet::Et_scheme;
  else if (rawscheme=="Et2") recom=fastjet::Et2_scheme;
  else if (rawscheme=="BIpt") recom=fastjet::BIpt_scheme;
  else if (rawscheme=="BIpt2") recom=fastjet::BIpt2_scheme;
  else THROW(fatal_error, "Unknown recombination scheme.");

  const auto R = s["R"].SetDefault(0.4).Get<double>();
  const auto f = s["f"].SetDefault(0.75).Get<double>();

  fastjet::Strategy strategy(fastjet::Best);
  const auto rawstrategy = s["Strategy"].SetDefault("Best").Get<std::string>();
  if (rawstrategy=="Best") strategy=fastjet::Best;
  else if (rawstrategy=="N2Plain") strategy=fastjet::N2Plain;
  else if (rawstrategy=="N2Tiled") strategy=fastjet::N2Tiled;
  else if (rawstrategy=="N2MinHeapTiled") strategy=fastjet::N2MinHeapTiled;
  else if (rawstrategy=="NlnN") strategy=fastjet::NlnN;
  else if (rawstrategy=="NlnNCam") strategy=fastjet::NlnNCam;
  else THROW(fatal_error, "Unknown strategy.");

  const auto btag = s["BTag"].SetDefault(0).Get<size_t>();

  if (siscone) {
    fastjet::JetDefinition::Plugin *plug(new fastjet::SISConePlugin(R,f));
    fastjet::JetDefinition jdef(plug);
    return new Fastjet_Interface(inlist,outlist,jdef,plug,njets,ptmin,etamax,btag);
  }
  fastjet::JetDefinition jdef(algo,R,recom,strategy);
  return new Fastjet_Interface(inlist,outlist,jdef,NULL,njets,ptmin,etamax,btag);
}

void ATOOLS::Getter
<Analysis_Object,Analysis_Key,Fastjet_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList: list,\n"
     <<std::setw(width+7)<<" "<<"OutList: list,\n"
     <<std::setw(width+7)<<" "<<"NJets: <n>  # (default 0 -> inclusive mode)\n"
     <<std::setw(width+7)<<" "<<"PTMin: <ptmin>  # (default 0)\n"
     <<std::setw(width+7)<<" "<<"EtaMax: <etamax>  # (default 1000.)\n"
     <<std::setw(width+7)<<" "<<"Algorithm: <algo>  # [kt|antikt|cambridge|siscone] (default kt)\n"
     <<std::setw(width+7)<<" "<<"Scheme: <scheme>  # [E|pt|pt2|Et|Et2|BIpt|BIpt2] (default E)\n"
     <<std::setw(width+7)<<" "<<"R: <R>  # (default 0.4)\n"
     <<std::setw(width+7)<<" "<<"f: <f>  # (siscone only, default 0.75)\n"
     <<std::setw(width+7)<<" "<<"Strategy: <strategy>  # [N2Plain|N2Tiled|N2MinHeapTiled|NlnN|NlnNCam|Best] (default Best)\n"
     <<std::setw(width+7)<<" "<<"BTag: <tag>  # 0|1 (default 0 -> no b-tag)\n"
     <<std::setw(width+4)<<" "<<"}";
}

#endif
