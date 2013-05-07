#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

using namespace ANALYSIS;
using namespace ATOOLS;

class Fastjet_Interface: public Trigger_Base {
private:

  fastjet::JetDefinition m_jdef;

  fastjet::JetDefinition::Plugin *p_plug;

  size_t m_njets, m_btag;
  double m_ptmin;

public:

  // constructor
  Fastjet_Interface(const std::string &inlist,
		    const std::string &outlist,
		    const fastjet::JetDefinition &jdef,
 		    fastjet::JetDefinition::Plugin *const plug,
		    const size_t &njets,const double &ptmin,
		    const size_t btag):
    Trigger_Base(inlist,outlist), m_jdef(jdef), p_plug(plug),
    m_njets(njets), m_ptmin(ptmin), m_btag(btag) {}

  ~Fastjet_Interface()
  {
    if (p_plug) delete p_plug;
  }

  // member functions
  Analysis_Object *GetCopy() const 
  {
    return new Fastjet_Interface
      (m_inlist,m_outlist,m_jdef,NULL,m_njets,m_ptmin,m_btag);
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
      jets=fastjet::sorted_by_pt(cs.inclusive_jets(m_ptmin));
    }
    std::vector<double> *ktdrs(new std::vector<double>());
    for (size_t i(input.size());i>0;--i)
      ktdrs->push_back(cs.exclusive_dmerge(i-1));
    std::string key("KtJetrates(1)"+m_outlist);
    p_ana->AddData(key,new Blob_Data<std::vector<double> *>(ktdrs));
    outlist.resize(jets.size());
    for (size_t i(0);i<outlist.size();++i) {
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
      outlist[i] = new Particle
	(1,Flavour(flav),Vec4D
	 (jets[i][3],jets[i][0],jets[i][1],jets[i][2]));
    }
    std::sort(outlist.begin(),outlist.end(),Order_PT());
  }

};

DECLARE_GETTER(Fastjet_Interface,"FastJets",
	       Analysis_Object,Argument_Matrix);	

Analysis_Object *ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Fastjet_Interface>::
operator()(const Argument_Matrix &parameters) const	
{
  if (parameters.size()<1) return NULL;
  fastjet::JetAlgorithm algo(fastjet::kt_algorithm);
  fastjet::RecombinationScheme recom(fastjet::E_scheme);
  fastjet::Strategy strategy(fastjet::Best);
  double R=0.4, f=0.75, p=1.0, ptmin=0.0;
  size_t njets=0, siscone=0, btag=0;
  std::string inlist="FinalState", outlist="FastJets";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Algorithm") {
      if (parameters[i][1]=="kt") algo=fastjet::kt_algorithm;
      if (parameters[i][1]=="cambridge") algo=fastjet::cambridge_algorithm;
      if (parameters[i][1]=="antikt") algo=fastjet::antikt_algorithm;
      if (parameters[i][1]=="siscone") siscone=1;
    }
    else if (parameters[i][0]=="Scheme") {
      if (parameters[i][1]=="E") recom=fastjet::E_scheme;
      if (parameters[i][1]=="pt") recom=fastjet::pt_scheme;
      if (parameters[i][1]=="pt2") recom=fastjet::pt2_scheme;
      if (parameters[i][1]=="Et") recom=fastjet::Et_scheme;
      if (parameters[i][1]=="Et2") recom=fastjet::Et2_scheme;
      if (parameters[i][1]=="BIpt") recom=fastjet::BIpt_scheme;
      if (parameters[i][1]=="BIpt2") recom=fastjet::BIpt2_scheme;
    }
    else if (parameters[i][0]=="R") R=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="f") f=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="p") p=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Strategy") {
      if (parameters[i][1]=="N2Plain") strategy=fastjet::N2Plain;
      if (parameters[i][1]=="N2Tiled") strategy=fastjet::N2Tiled;
      if (parameters[i][1]=="N2MinHeapTiled") strategy=fastjet::N2MinHeapTiled;
      if (parameters[i][1]=="NlnN") strategy=fastjet::NlnN;
      if (parameters[i][1]=="NlnNCam") strategy=fastjet::NlnNCam;
      if (parameters[i][1]=="Best") strategy=fastjet::Best;
    }
    else if (parameters[i][0]=="NJets") njets=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="PTMin") ptmin=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="BTag") btag=ATOOLS::ToType<int>(parameters[i][1]);
  }
  if (siscone) {
    fastjet::JetDefinition::Plugin *plug(new fastjet::SISConePlugin(R,f));
    fastjet::JetDefinition jdef(plug);
    return new Fastjet_Interface(inlist,outlist,jdef,plug,njets,ptmin,btag);
  }
  fastjet::JetDefinition jdef(algo,R,recom,strategy);
  return new Fastjet_Interface(inlist,outlist,jdef,NULL,njets,ptmin,btag);
}									

void ATOOLS::Getter
<Analysis_Object,Argument_Matrix,Fastjet_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList    list\n"
     <<std::setw(width+7)<<" "<<"OutList   list\n"
     <<std::setw(width+7)<<" "<<"NJets     #jets (default 0 -> inclusive mode)\n"
     <<std::setw(width+7)<<" "<<"PTMin     ptmin (default 0)\n"
     <<std::setw(width+7)<<" "<<"Algorithm algorithm [kt|antikt|cambridge|siscone] (default kt)\n"
     <<std::setw(width+7)<<" "<<"Scheme    scheme [E|pt|pt2|Et|Et2|BIpt|BIpt2] (default E)\n"
     <<std::setw(width+7)<<" "<<"R         R (default 0.4)\n"
     <<std::setw(width+7)<<" "<<"p         p (default 1.0)\n"
     <<std::setw(width+7)<<" "<<"f         f (siscone only, default 0.75)\n"
     <<std::setw(width+7)<<" "<<"Strategy  strategy [N2Plain|N2Tiled|N2MinHeapTiled|NlnN|NlnNCam|Best] (default Best)\n"
     <<std::setw(width+7)<<" "<<"BTag      0|1 (default 0 -> no b-tag)\n"
     <<std::setw(width+4)<<" "<<"}";
}

#endif
