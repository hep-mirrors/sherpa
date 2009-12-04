#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "Trigger_Base.H"
#include "Primitive_Analysis.H"
#include "Kt_Algorithm.H"
#include "MyStrStream.H"
#include "Exception.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace ANALYSIS;
using namespace ATOOLS;

class Fastjet_Interface: public Trigger_Base {
private:

  fastjet::JetDefinition m_jdef;

  size_t m_njets;
  double m_ptmin;

public:

  // constructor
  Fastjet_Interface(const std::string &inlist,
		    const std::string &outlist,
		    const fastjet::JetDefinition &jdef,
		    const size_t &njets,const double &ptmin):
    Trigger_Base(inlist,outlist), m_jdef(jdef),
    m_njets(njets), m_ptmin(ptmin) {}

  // member functions
  Analysis_Object *GetCopy() const 
  {
    return new Fastjet_Interface
      (m_inlist,m_outlist,m_jdef,m_njets,m_ptmin);
  }

  void Evaluate(const ATOOLS::Particle_List &plist,
		ATOOLS::Particle_List &outlist,
		double value,int ncount)
  {
    std::vector<fastjet::PseudoJet> input(plist.size()), jets;
    for (size_t i(0);i<input.size();++i) {
      Vec4D p(plist[i]->Momentum());
      input[i]=fastjet::PseudoJet(p[1],p[2],p[3],p[0]);
    }
    fastjet::ClusterSequence cs(input,m_jdef);
    if (m_njets>0) {
      jets=cs.exclusive_jets((int)m_njets);
    }
    else {
      jets=cs.inclusive_jets(m_ptmin);
    }
    outlist.resize(jets.size());
    for (size_t i(0);i<outlist.size();++i) {
      outlist[i] = new Particle
	(1,Flavour(kf_jet),Vec4D
	 (jets[i][3],jets[i][0],jets[i][1],jets[i][2]));
    }
    std::sort(outlist.begin(),outlist.end(),Order_PT());
  }

};

DECLARE_GETTER(FastJet_Getter,"FastJets",
	       Analysis_Object,Argument_Matrix);	

Analysis_Object *FastJet_Getter::
operator()(const Argument_Matrix &parameters) const	
{
  if (parameters.size()<1) return NULL;
  fastjet::JetAlgorithm algo(fastjet::kt_algorithm);
  fastjet::RecombinationScheme recom(fastjet::E_scheme);
  fastjet::Strategy strategy(fastjet::Best);
  double R=0.4, p=1.0, ptmin=0.0;
  size_t njets=0;
  std::string inlist="FinalState", outlist="FastJets";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Algorithm") {
      if (parameters[i][1]=="kt") algo=fastjet::kt_algorithm;
      if (parameters[i][1]=="cambridge") algo=fastjet::cambridge_algorithm;
      if (parameters[i][1]=="antikt") algo=fastjet::antikt_algorithm;
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
  }
  fastjet::JetDefinition jdef(algo,R,recom,strategy);
  return new Fastjet_Interface(inlist,outlist,jdef,njets,ptmin);
}									

void FastJet_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList    list\n"
     <<std::setw(width+7)<<" "<<"OutList   list\n"
     <<std::setw(width+7)<<" "<<"NJets     #jets (0 for inclusive mode)\n"
     <<std::setw(width+7)<<" "<<"PTMin     ptmin (default is 0)\n"
     <<std::setw(width+7)<<" "<<"Algorithm algorithm\n"
     <<std::setw(width+7)<<" "<<"Scheme    scheme\n"
     <<std::setw(width+7)<<" "<<"R         R\n"
     <<std::setw(width+7)<<" "<<"p         p\n"
     <<std::setw(width+7)<<" "<<"Strategy  strategy\n"
     <<std::setw(width+7)<<"}";
}

#endif
