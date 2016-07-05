#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"

namespace CSSHOWER {
  
  class CF_EW_FFZ: public SF_Coupling {
  protected:

    ATOOLS::Function_Base *p_cpl;
    ATOOLS::Flavour m_cfl;

    double m_q[2];

  public:

    inline CF_EW_FFZ(const SF_Key &key):
      SF_Coupling(key), m_cfl(key.p_v->in[0].Bar()) 
    {
      if (key.m_type==cstp::IF || key.m_type==cstp::II)
	m_cfl=key.p_v->in[key.m_mode==0?1:2];
    }

    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int pol);
    bool AllowSpec(const ATOOLS::Flavour &fl);

  };

  class CF_EW_FFW: public SF_Coupling {
  protected:

    ATOOLS::Function_Base *p_cpl;
    ATOOLS::Flavour m_cfl;

    double m_q[2];

  public:

    inline CF_EW_FFW(const SF_Key &key):
      SF_Coupling(key), m_cfl(key.p_v->in[0].Bar()) 
    {
      if (key.m_type==cstp::IF || key.m_type==cstp::II)
	m_cfl=key.p_v->in[key.m_mode==0?1:2];
    }

    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int pol);
    bool AllowSpec(const ATOOLS::Flavour &fl);

  };

}

using namespace CSSHOWER;
using namespace ATOOLS;

bool CF_EW_FFZ::SetCoupling(MODEL::Model_Base *md,
			    const double &k0sqi,const double &k0sqf,
			    const double &isfac,const double &fsfac)
{
  double stw(std::abs(md->ComplexConstant("csin2_thetaW")));
  Flavour ffl(p_lf->FlB().IsFermion()?p_lf->FlB():p_lf->FlC());
  if (!ffl.IsFermion()) THROW(fatal_error,"Internal error");
  if (ffl.IsAnti()) ffl=ffl.Bar();
  double af(ffl.IsoWeak()), vf(af-2.0*ffl.Charge()*stw);
  m_q[0]=0.25/(stw*(1.0-stw))*(sqr(vf)+sqr(af));
  m_q[1]=2.0/stw*sqr(af*ffl.Mass()/Flavour(kf_Wplus).Mass());
  p_cpl=md->GetScalarFunction("alpha_QED");
  m_cplfac=1.0;
  double cqed((*p_cpl)(sqr(rpa->gen.Ecms())));
  m_cplmax.push_back(cqed*m_q[0]);
  m_cplmax.push_back(cqed*m_q[1]);
  return true;
}

double CF_EW_FFZ::Coupling(const double &scale,const int pol)
{
  if (pol>1) return 0.0;
  if (scale<0.0) return m_last = m_cplmax.front()*m_q[pol];
  double scl(CplFac(scale)*scale);
  return m_last = (*p_cpl)(scl)*m_q[pol];
}

bool CF_EW_FFZ::AllowSpec(const ATOOLS::Flavour &fl) 
{
  if (m_cfl.IntCharge()==0) return fl.Charge();

  switch (m_type) {
  case cstp::FF:
  case cstp::II:
    return fl.IntCharge()*m_cfl.IntCharge()<0;
  case cstp::FI:
  case cstp::IF:
    return fl.IntCharge()*m_cfl.IntCharge()>0;
  default:
    return false;
  }
}

bool CF_EW_FFW::SetCoupling(MODEL::Model_Base *md,
			    const double &k0sqi,const double &k0sqf,
			    const double &isfac,const double &fsfac)
{
  double stw(std::abs(md->ComplexConstant("csin2_thetaW")));
  Complex vij(Complex(1.0,0.0));
  Flavour f1(p_lf->FlB()), f2(p_lf->FlC());
  if (!f1.IsFermion()) f1=p_lf->FlA();
  else if (!f2.IsFermion()) f2=p_lf->FlA();
  if (f1.IsQuark()) {
    if (f1.IsDowntype()) std::swap<Flavour>(f1,f2);
    //int i((int)(f1.Kfcode())), j((int)(f2.Kfcode()));
    if (md->Name().find("SM")==std::string::npos) vij=1.0;
    // else vij=md->ComplexMatrixElement("CKM",i/2-1,(j-1)/2);
  }
  else {
    if (f1.Kfcode()%2==0) std::swap<Flavour>(f1,f2);
  }
  double vf(sqr(std::abs(vij)));
  m_q[0]=0.5/stw*vf;
  m_q[1]=1.0/stw*vf*sqr(f1.Mass()/Flavour(kf_Wplus).Mass());
  p_cpl=md->GetScalarFunction("alpha_QED");
  m_cplfac=1.0;
  double cqed((*p_cpl)(sqr(rpa->gen.Ecms())));
  m_cplmax.push_back(cqed*m_q[0]);
  m_cplmax.push_back(cqed*m_q[1]);
  return m_q[0]>0.0;
}

double CF_EW_FFW::Coupling(const double &scale,const int pol)
{
  if (pol>1) return 0.0;
  if (scale<0.0) return m_cplmax.front()*m_q[pol];
  double scl(CplFac(scale)*scale);
  return (*p_cpl)(scl)*m_q[pol];
}

bool CF_EW_FFW::AllowSpec(const ATOOLS::Flavour &fl) 
{
  return true;
}

namespace CSSHOWER {

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

}

DECLARE_GETTER(CF_EW_FFW_Getter,"SF_EW_FFV_Fill",
	       void,SFC_Filler_Key);

void *ATOOLS::Getter<void,SFC_Filler_Key,CF_EW_FFW_Getter>::
operator()(const SFC_Filler_Key &key) const
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
      key.p_gets->push_back(new CF_EW_FFZ_Getter(ptag+qtag+qbtag));
      key.p_gets->push_back(new CF_EW_FFZ_Getter(qbtag+qbtag+ptag));
      key.p_gets->push_back(new CF_EW_FFZ_Getter(qtag+qtag+ptag));
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
	key.p_gets->push_back(new CF_EW_FFW_Getter(f1tag+f2tag+wmtag));
	key.p_gets->push_back(new CF_EW_FFW_Getter(f2tag+f1tag+wptag));
	key.p_gets->push_back(new CF_EW_FFW_Getter(f2btag+wmtag+f1btag));
	key.p_gets->push_back(new CF_EW_FFW_Getter(f1btag+wptag+f2btag));
	key.p_gets->push_back(new CF_EW_FFW_Getter(wptag+f1btag+f2tag));
	key.p_gets->push_back(new CF_EW_FFW_Getter(wmtag+f2btag+f1tag));
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
      key.p_gets->push_back(new CF_EW_FFW_Getter(f1tag+f2tag+wmtag));
      key.p_gets->push_back(new CF_EW_FFW_Getter(f2tag+f1tag+wptag));
      key.p_gets->push_back(new CF_EW_FFW_Getter(f2btag+wmtag+f1btag));
      key.p_gets->push_back(new CF_EW_FFW_Getter(f1btag+wptag+f2btag));
      key.p_gets->push_back(new CF_EW_FFW_Getter(wptag+f1btag+f2tag));
      key.p_gets->push_back(new CF_EW_FFW_Getter(wmtag+f2btag+f1tag));
    }
  }
  return NULL;
}

void ATOOLS::Getter<void,SFC_Filler_Key,CF_EW_FFW_Getter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"ew ffv coupling filler";
}
