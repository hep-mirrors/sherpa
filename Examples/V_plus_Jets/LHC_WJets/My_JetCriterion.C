#include "PDF/Main/Jet_Criterion.H"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

namespace MYSTUFF {

  class FastJet_Jet_Criterion: public Jet_Criterion {
  private:

    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;

    double m_y;

  public:

    FastJet_Jet_Criterion(const std::string &args): p_siscplug(NULL)
    {
      std::string jtag(args);
      size_t pos(jtag.find("FASTJET["));
      if (pos==std::string::npos)
	THROW(fatal_error,"Invalid scale '"+args+"'");
      jtag=jtag.substr(pos+8);
      pos=jtag.find(']');
      if (pos==std::string::npos)
	THROW(fatal_error,"Invalid scale '"+args+"'");
      jtag=jtag.substr(0,pos);
      Data_Reader read(" ",",","#","=");
      read.AddIgnore(":");
      read.SetAddCommandLine(false);
      read.SetString(jtag);
      m_y=read.StringValue<double>("Y",100.0);
      double R(read.StringValue<double>("R",0.4));
      double f(read.StringValue<double>("f",0.75));
      std::string algo(read.StringValue<std::string>("A","antikt"));
      fastjet::JetAlgorithm ja(fastjet::kt_algorithm);
      if (algo=="cambridge") ja=fastjet::cambridge_algorithm;
      if (algo=="antikt") ja=fastjet::antikt_algorithm;
      if (algo=="siscone") p_siscplug=new fastjet::SISConePlugin(R,f);
      std::string reco(read.StringValue<std::string>("C","E"));
      fastjet::RecombinationScheme recom(fastjet::E_scheme);
      if (reco=="pt") recom=fastjet::pt_scheme;
      if (reco=="pt2") recom=fastjet::pt2_scheme;
      if (reco=="Et") recom=fastjet::Et_scheme;
      if (reco=="Et2") recom=fastjet::Et2_scheme;
      if (reco=="BIpt") recom=fastjet::BIpt_scheme;
      if (reco=="BIpt2") recom=fastjet::BIpt2_scheme;
      if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
      else p_jdef=new fastjet::JetDefinition(ja,R,recom);
    }

    ~FastJet_Jet_Criterion()
    {
      if (p_siscplug) delete p_siscplug;
      delete p_jdef;
    }

    bool Jets(Cluster_Amplitude *ampl,int mode)
    {
      int nj=ampl->NIn();
      std::vector<fastjet::PseudoJet> input,jets;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
	Vec4D p(ampl->Leg(i)->Mom());
	if (Flavour(kf_jet).Includes(ampl->Leg(i)->Flav()))
	  input.push_back(fastjet::PseudoJet(p[1],p[2],p[3],p[0])); 
	else ++nj;
      }
      double pt2(ampl->JF<Jet_Finder>()->Ycut()*sqr(rpa->gen.Ecms()));
      fastjet::ClusterSequence cs(input,*p_jdef);
      jets=fastjet::sorted_by_pt(cs.inclusive_jets());
      for (size_t i(0);i<jets.size();++i) {
	Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
	if (pj.PPerp2()>pt2 && (m_y==100 || dabs(pj.Y())<m_y)) ++nj;
      }
      return nj==ampl->Legs().size();
    }

  };// end of class FastJet_Jet_Criterion

  DECLARE_GETTER(FastJet_Jet_Criterion_Getter,"FASTJET",
		 Jet_Criterion,JetCriterion_Key);
  Jet_Criterion *FastJet_Jet_Criterion_Getter::
  operator()(const JetCriterion_Key &args) const
  { return new FastJet_Jet_Criterion(args.m_key); }
  void FastJet_Jet_Criterion_Getter::
  PrintInfo(std::ostream &str,const size_t width) const
  { str<<"The FastJet jet criterion"; }

}
