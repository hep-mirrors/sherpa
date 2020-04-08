#include "PDF/Main/Jet_Criterion.H"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PDF;
using namespace ATOOLS;

namespace PHASIC {

  class FastJet_Jet_Criterion: public Jet_Criterion {
  private:

    fjcore::JetDefinition * p_jdef;

    double m_y;

  public:

    FastJet_Jet_Criterion(const std::string &args)
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
      read.SetString(jtag);
      m_y=read.StringValue<double>("Y",100.0);
      double R(read.StringValue<double>("R",0.4));
      double f(read.StringValue<double>("f",0.75));
      std::string algo(read.StringValue<std::string>("A","antikt"));
      fjcore::JetAlgorithm ja(fjcore::kt_algorithm);
      if (algo=="cambridge") ja=fjcore::cambridge_algorithm;
      if (algo=="antikt") ja=fjcore::antikt_algorithm;
      std::string reco(read.StringValue<std::string>("C","E"));
      fjcore::RecombinationScheme recom(fjcore::E_scheme);
      if (reco=="pt") recom=fjcore::pt_scheme;
      if (reco=="pt2") recom=fjcore::pt2_scheme;
      if (reco=="Et") recom=fjcore::Et_scheme;
      if (reco=="Et2") recom=fjcore::Et2_scheme;
      if (reco=="BIpt") recom=fjcore::BIpt_scheme;
      if (reco=="BIpt2") recom=fjcore::BIpt2_scheme;
      p_jdef=new fjcore::JetDefinition(ja,R,recom);
    }

    ~FastJet_Jet_Criterion()
    {
      delete p_jdef;
    }

    double Value(Cluster_Amplitude *ampl,int mode)
    {
      int nj=ampl->NIn();
      std::vector<fjcore::PseudoJet> input,jets;
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
	Vec4D p(ampl->Leg(i)->Mom());
	if (Flavour(kf_jet).Includes(ampl->Leg(i)->Flav()))
	  input.push_back(fjcore::PseudoJet(p[1],p[2],p[3],p[0])); 
	else ++nj;
      }
      double pt2(sqr(ampl->JF<Jet_Finder>()->Qcut()));
      fjcore::ClusterSequence cs(input,*p_jdef);
      jets=fjcore::sorted_by_pt(cs.inclusive_jets());
      for (size_t i(0);i<jets.size();++i) {
	Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
	if (pj.PPerp2()>pt2 && (m_y==100 || dabs(pj.Y())<m_y)) ++nj;
      }
      return nj+mode>=ampl->Legs().size()?2.0*pt2:0.0;
    }

  };// end of class FastJet_Jet_Criterion

}

using namespace PHASIC;

DECLARE_GETTER(FastJet_Jet_Criterion,"FASTJET",
	       Jet_Criterion,JetCriterion_Key);
Jet_Criterion *ATOOLS::Getter
<Jet_Criterion,JetCriterion_Key,FastJet_Jet_Criterion>::
operator()(const JetCriterion_Key &args) const
{ return new FastJet_Jet_Criterion(args.m_key); }
void ATOOLS::Getter
<Jet_Criterion,JetCriterion_Key,FastJet_Jet_Criterion>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"The FastJet jet criterion"; }
