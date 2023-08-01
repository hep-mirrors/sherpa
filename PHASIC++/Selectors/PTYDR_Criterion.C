#include "PDF/Main/Jet_Criterion.H"

#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"

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
      read.SetString(jtag);
      m_y=read.StringValue<double>("Y",6.0);
      m_r=read.StringValue<double>("R",0.4);
    }

    double Value(Cluster_Amplitude *ampl,int mode)
    {
      int nj=ampl->NIn();
      std::vector<Vec4D> input;
      for (size_t i(ampl->NIn()); i < ampl->Legs().size(); ++i) {
        Vec4D p(ampl->Leg(i)->Mom());
        if (Flavour(kf_jet).Includes(ampl->Leg(i)->Flav()))
          input.push_back(p);
        else
          ++nj;
      }
      double pt2(sqr(ampl->JF<Jet_Finder>()->Qcut()));
      for (size_t i(0);i<input.size();++i) {
        if (input[i].PPerp2() > pt2 && dabs(input[i].Y()) < m_y) {
          for (size_t j(i + 1); j < input.size(); ++j) {
            if (input[i].DR(input[j]) > m_r)
              nj++;
          }
        }
      }
      return nj+mode>=ampl->Legs().size()?2.0*pt2:0.0;
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
