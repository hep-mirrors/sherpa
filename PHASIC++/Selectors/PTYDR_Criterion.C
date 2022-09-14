#include "PDF/Main/Jet_Criterion.H"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

namespace PHASIC {

  class PTYDR_Jet_Criterion: public Jet_Criterion {
  private:

    double m_y, m_r;

  public:

    PTYDR_Jet_Criterion(const std::string &args)
    {
      std::string jtag(args);
      size_t pos(jtag.find("PTYDR["));
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
      m_y=read.StringValue<double>("Y",6.0);
      m_r=read.StringValue<double>("R",0.4);
    }

    bool Jets(Cluster_Amplitude *ampl,int mode)
    {
      double pt2(ampl->JF<Jet_Finder>()->Ycut()*sqr(rpa->gen.Ecms()));
      for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
	if (Flavour(kf_jet).Includes(ampl->Leg(i)->Flav())) {
	  Vec4D pi(ampl->Leg(i)->Mom());
	  if (pi.PPerp2()<pt2 || dabs(pi.Y())>m_y) return false;
	  for (size_t j(i+1);j<ampl->Legs().size();++j)
	    if (Flavour(kf_jet).Includes(ampl->Leg(j)->Flav())) {
	      Vec4D pj(ampl->Leg(j)->Mom());
	      if (pi.DR(pj)<m_r) return false;
	    }
	}
      return true;
    }

  };// end of class PTYDR_Jet_Criterion

}

using namespace PHASIC;

DECLARE_GETTER(PTYDR_Jet_Criterion,"PTYDR",
	       Jet_Criterion,JetCriterion_Key);
Jet_Criterion *ATOOLS::Getter
<Jet_Criterion,JetCriterion_Key,PTYDR_Jet_Criterion>::
operator()(const JetCriterion_Key &args) const
{ return new PTYDR_Jet_Criterion(args.m_key); }
void ATOOLS::Getter
<Jet_Criterion,JetCriterion_Key,PTYDR_Jet_Criterion>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"The PTYDR jet criterion"; }
