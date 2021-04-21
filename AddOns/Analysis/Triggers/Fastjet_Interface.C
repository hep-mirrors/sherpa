#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace ANALYSIS;
using namespace ATOOLS;

class Fastjet_Interface: public Trigger_Base {
private:

  fjcore::JetDefinition m_jdef;

  size_t m_njets, m_btag;
  double m_ptmin, m_etamax;

public:

  // constructor
  Fastjet_Interface(const std::string &inlist,
		    const std::string &outlist,
		    const fjcore::JetDefinition &jdef,
		    fjcore::JetDefinition::Plugin *const plug,
		    const size_t &njets,const double &ptmin,
		    const double &etamax,const size_t btag):
    Trigger_Base(inlist,outlist), m_jdef(jdef),
    m_njets(njets), m_btag(btag), m_ptmin(ptmin), m_etamax(etamax) {}

  ~Fastjet_Interface()
  {
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
    std::vector<fjcore::PseudoJet> input(plist.size()), jets;
    for (size_t i(0);i<input.size();++i) {
      Vec4D p(plist[i]->Momentum());
      input[i]=fjcore::PseudoJet(p[1],p[2],p[3],p[0]);
      input[i].set_user_index(BTag(plist[i]));
    }
    fjcore::ClusterSequence cs(input,m_jdef);
    if (m_njets>0) {
      jets=fjcore::sorted_by_pt(cs.exclusive_jets((int)m_njets));
    }
    else {
      jets=fjcore::sorted_by_pt(cs.inclusive_jets());
    }
    std::vector<double> *ktdrs(new std::vector<double>());
    for (size_t i(input.size());i>0;--i) {
      if      (m_jdef.jet_algorithm()==fjcore::kt_algorithm)
        ktdrs->push_back(cs.exclusive_dmerge(i-1));
      else if (m_jdef.jet_algorithm()==fjcore::antikt_algorithm)
        ktdrs->insert(ktdrs->begin(),1./cs.exclusive_dmerge(i-1));
    }
    std::string key("KtJetrates(1)"+m_outlist);
    p_ana->AddData(key,new Blob_Data<std::vector<double> *>(ktdrs));
    for (size_t i(0);i<jets.size();++i) {
      kf_code flav(kf_jet);
      if (m_btag) {
	int nb(0);
	const std::vector<fjcore::PseudoJet>
	  &cons(jets[i].constituents());
	for (size_t j=0;j<cons.size();++j) {
	  if (cons[j].user_index()==5) ++nb;
	  if (cons[j].user_index()==-5) --nb;
	}
	if (nb!=0) flav=kf_bjet;
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

  fjcore::JetAlgorithm algo;
  const auto rawalgorithm = s["Algorithm"].SetDefault("kt").Get<std::string>();
  if (rawalgorithm=="kt") algo=fjcore::kt_algorithm;
  else if (rawalgorithm=="cambridge") algo=fjcore::cambridge_algorithm;
  else if (rawalgorithm=="antikt") algo=fjcore::antikt_algorithm;
  else THROW(fatal_error, "Unknown jet algorithm.");

  fjcore::RecombinationScheme recom;
  const auto rawscheme = s["Scheme"].SetDefault("E").Get<std::string>();
  if (rawscheme=="E") recom=fjcore::E_scheme;
  else if (rawscheme=="pt") recom=fjcore::pt_scheme;
  else if (rawscheme=="pt2") recom=fjcore::pt2_scheme;
  else if (rawscheme=="Et") recom=fjcore::Et_scheme;
  else if (rawscheme=="Et2") recom=fjcore::Et2_scheme;
  else if (rawscheme=="BIpt") recom=fjcore::BIpt_scheme;
  else if (rawscheme=="BIpt2") recom=fjcore::BIpt2_scheme;
  else THROW(fatal_error, "Unknown recombination scheme.");

  const auto R = s["R"].SetDefault(0.4).Get<double>();
  const auto f = s["f"].SetDefault(0.75).Get<double>();

  fjcore::Strategy strategy(fjcore::Best);
  const auto rawstrategy = s["Strategy"].SetDefault("Best").Get<std::string>();
  if (rawstrategy=="Best") strategy=fjcore::Best;
  else if (rawstrategy=="N2Plain") strategy=fjcore::N2Plain;
  else if (rawstrategy=="N2Tiled") strategy=fjcore::N2Tiled;
  else if (rawstrategy=="N2MinHeapTiled") strategy=fjcore::N2MinHeapTiled;
  else if (rawstrategy=="NlnN") strategy=fjcore::NlnN;
  else if (rawstrategy=="NlnNCam") strategy=fjcore::NlnNCam;
  else THROW(fatal_error, "Unknown strategy.");

  const auto btag = s["BTag"].SetDefault(0).Get<size_t>();

  fjcore::JetDefinition jdef(algo,R,recom,strategy);
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
