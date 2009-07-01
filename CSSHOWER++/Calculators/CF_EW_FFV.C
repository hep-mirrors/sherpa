#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"

namespace CSSHOWER {
  
  class CF_EW_FFZ: public SF_Coupling {
  protected:

    ATOOLS::Function_Base *p_cpl;
    ATOOLS::Flavour m_cfl, m_ffl;

    double m_q[2];

  public:

    inline CF_EW_FFZ(const SF_Key &key):
      SF_Coupling(key), m_cfl(key.p_v->in[0]) 
    {
      m_ffl=m_cfl.Kfcode()==kf_Z?key.p_v->in[1]:m_cfl;
    }

    bool SetCoupling(MODEL::Model_Base *md,const double &k0sq,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int mode);
    bool AllowSpec(const ATOOLS::Flavour &fl);

  };

  class CF_EW_FFW: public SF_Coupling {
  protected:

    ATOOLS::Function_Base *p_cpl;
    ATOOLS::Flavour m_cfl, m_ffl;

    double m_q[2];

  public:

    inline CF_EW_FFW(const SF_Key &key):
      SF_Coupling(key), m_cfl(key.p_v->in[0]) 
    {
      m_ffl=m_cfl.Kfcode()==kf_Wplus?key.p_v->in[1]:m_cfl;
    }

    bool SetCoupling(MODEL::Model_Base *md,const double &k0sq,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int mode);
    bool AllowSpec(const ATOOLS::Flavour &fl);

  };

}

using namespace CSSHOWER;
using namespace ATOOLS;

bool CF_EW_FFZ::SetCoupling(MODEL::Model_Base *md,const double &k0sq,
			    const double &isfac,const double &fsfac)
{
  double stw(md->GetInteractionModel()->ScalarConstant("sin2_thetaW"));
  double af(m_ffl.IsoWeak()), vf(af-2.0*m_ffl.Charge()*stw);
  m_q[0]=0.25/(stw*(1.0-stw))*(sqr(vf)+sqr(af));
  m_q[1]=2.0/stw*sqr(af*m_ffl.Mass()/Flavour(kf_Wplus).Mass());
  m_cplfac=(m_type/10==1)?fsfac:isfac;
  p_cpl=md->GetScalarFunction("alpha_QED");
  double cqed((*p_cpl)(m_cplfac*rpa.gen.CplScale()));
  m_cplmax.push_back(cqed*m_q[0]);
  m_cplmax.push_back(cqed*m_q[1]);
  return true;
}

double CF_EW_FFZ::Coupling(const double &scale,const int mode)
{
  if (mode>1) return 0.0;
  return (*p_cpl)(m_cplfac*scale)*m_q[mode];
}

bool CF_EW_FFZ::AllowSpec(const ATOOLS::Flavour &fl) 
{
  if (m_cfl.IntCharge()==0) return fl.Charge();
  return fl.IntCharge()*m_cfl.IntCharge()<0;
}

bool CF_EW_FFW::SetCoupling(MODEL::Model_Base *md,const double &k0sq,
			    const double &isfac,const double &fsfac)
{
  double stw(md->GetInteractionModel()->ScalarConstant("sin2_thetaW"));
  Complex vij(Complex(1.0,0.0));
  Flavour f1(p_lf->FlB()), f2(p_lf->FlC());
  if (!f1.IsFermion()) f1=p_lf->FlA();
  else if (!f2.IsFermion()) f2=p_lf->FlA();
  if (f1.IsQuark()) {
    int i((int)(f1.Kfcode())), j((int)(f2.Kfcode()));
    if (f1.IsDowntype()) std::swap<int>(i,j);
    vij=md->ComplexMatrixElement("CKM",i/2-1,(j-1)/2);
  }
  double vf(sqr(std::abs(vij)));
  m_q[0]=0.5/stw*sqr(vf);
  m_q[1]=1.0/stw*sqr(vf*m_ffl.Mass()/Flavour(kf_Wplus).Mass());
  m_cplfac=(m_type/10==1)?fsfac:isfac;
  p_cpl=md->GetScalarFunction("alpha_QED");
  double cqed((*p_cpl)(m_cplfac*rpa.gen.CplScale()));
  m_cplmax.push_back(cqed*m_q[0]);
  m_cplmax.push_back(cqed*m_q[1]);
  return m_q[0]>0.0;
}

double CF_EW_FFW::Coupling(const double &scale,const int mode)
{
  if (mode>1) return 0.0;
  return (*p_cpl)(m_cplfac*scale)*m_q[mode];
}

bool CF_EW_FFW::AllowSpec(const ATOOLS::Flavour &fl) 
{
  return true;
}

DECLARE_CPL_GETTER(CF_EW_FFZ_Getter);

SF_Coupling *CF_EW_FFZ_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_EW_FFZ(args);
}

void CF_EW_FFZ_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ffz coupling";
}

DECLARE_CPL_GETTER(CF_EW_FFW_Getter);

SF_Coupling *CF_EW_FFW_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_EW_FFW(args);
}

void CF_EW_FFW_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ffw coupling";
}

DECLARE_GETTER(CF_EW_FFV_Filler,"SF_EW_FFV_Fill",
	       void,const MODEL::Model_Base *);

void *CF_EW_FFV_Filler::operator()
  (const MODEL::Model_Base *const &model) const
{
  if (!Flavour(kf_Z).IsOn()) return NULL;
  std::string ptag("{"+Flavour(kf_Z).IDName()+"}");
  std::string wptag("{"+Flavour(kf_Wplus).IDName()+"}");
  std::string wmtag("{"+Flavour(kf_Wplus).Bar().IDName()+"}");
  // ffz couplings
  if (Flavour(kf_Z).IsOn()) {
    for (int i(1);i<=16;++i) {
      if (i==7) i=11;
      Flavour f((kf_code)i);
      if (!f.IsOn()) continue;
      std::string qtag("{"+f.IDName()+"}");
      std::string qbtag ("{"+f.Bar().IDName()+"}");
      new CF_EW_FFZ_Getter(ptag+qtag+qbtag);
      new CF_EW_FFZ_Getter(qbtag+qbtag+ptag);
      new CF_EW_FFZ_Getter(qtag+qtag+ptag);
    }
  }
  if (Flavour(kf_Wplus).IsOn()) {
    // qqw couplings
    for (int i(1);i<=5;i+=2) {
      for (int j(2);j<=6;j+=2) {
	Flavour f1((kf_code)i), f2((kf_code)j);
	if (!f1.IsOn() || !f2.IsOn()) continue;
	std::string f1tag("{"+f1.IDName()+"}");
	std::string f2tag("{"+f2.IDName()+"}");
	std::string f1btag("{"+f1.Bar().IDName()+"}");
	std::string f2btag("{"+f2.Bar().IDName()+"}");
	new CF_EW_FFW_Getter(f1tag+f2tag+wmtag);
	new CF_EW_FFW_Getter(f2tag+f1tag+wptag);
	new CF_EW_FFW_Getter(f2btag+wmtag+f1btag);
	new CF_EW_FFW_Getter(f1btag+wptag+f2btag);
	new CF_EW_FFW_Getter(wptag+f1btag+f2tag);
	new CF_EW_FFW_Getter(wmtag+f2btag+f1tag);
      }
    }
    // llw couplings
    for (int i(11);i<=16;i+=2) {
      Flavour f1((kf_code)i), f2((kf_code)(i+1));
      if (!f1.IsOn() || !f2.IsOn()) continue;
      std::string f1tag("{"+f1.IDName()+"}");
      std::string f2tag("{"+f2.IDName()+"}");
      std::string f1btag("{"+f1.Bar().IDName()+"}");
      std::string f2btag("{"+f2.Bar().IDName()+"}");
      new CF_EW_FFW_Getter(f1tag+f2tag+wmtag);
      new CF_EW_FFW_Getter(f2tag+f1tag+wptag);
      new CF_EW_FFW_Getter(f2btag+wmtag+f1btag);
      new CF_EW_FFW_Getter(f1btag+wptag+f2btag);
      new CF_EW_FFW_Getter(wptag+f1btag+f2tag);
      new CF_EW_FFW_Getter(wmtag+f2btag+f1tag);
    }
  }
  return NULL;
}

void CF_EW_FFV_Filler::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ew ffv coupling filler";
}
