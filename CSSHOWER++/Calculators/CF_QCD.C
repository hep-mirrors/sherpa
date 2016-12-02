#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include <algorithm>
#include <cassert>

using namespace ATOOLS;

namespace CSSHOWER {
  
  const double s_Nc = 3.;
  const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  const double s_CA = s_Nc;
  const double s_TR = 1./2.;

  class CF_QCD: public SF_Coupling {
  protected:

    class QCD_Coupling_Info {

    public:

      QCD_Coupling_Info():
      p_cpl(NULL), m_rsf(1.0), p_cplmax(NULL)
      {};
      QCD_Coupling_Info(MODEL::Running_AlphaS * cpl, double rsf, std::vector<double> * cplmax=NULL):
      p_cpl(cpl), m_rsf(rsf), p_cplmax(cplmax)
      {};

      MODEL::Running_AlphaS * const Coupling() const { return p_cpl; }
      double RSF() const { return m_rsf; }
      std::vector<double> * const MaxCoupling() const { return p_cplmax; }
      void SetMaxCoupling(std::vector<double> * const cplmax) { p_cplmax = cplmax; }
      bool IsValid() const { return (p_cpl != NULL); }

    private:

      MODEL::Running_AlphaS * p_cpl;
      double m_rsf;
      std::vector<double> * p_cplmax;
    };

    /*!
     * Underlying couplings set by SetCoupling and SetAlternativeUnderlyingCoupling, respectively
     *
     * If the alternative coupling is set, it takes precedence. This can be used for reweighting purposes.
     */
    QCD_Coupling_Info m_maincplinfo, m_altcplinfo;
    const QCD_Coupling_Info & CurrentCouplingInfo() const
    { return (m_altcplinfo.IsValid() ? m_altcplinfo : m_maincplinfo); }

    //! Buffer of max alphas values to avoid re-calculations
    std::map<MODEL::Running_AlphaS *, std::vector<double> > m_altcplmax;

    double m_q, m_k0sq;
    int m_scvmode;

  public:

    CF_QCD(const SF_Key &key):
      SF_Coupling(key),
      m_maincplinfo(QCD_Coupling_Info()), m_altcplinfo(QCD_Coupling_Info()),
      m_altcplmax(),
      m_q(0.), m_k0sq(0.0), m_scvmode(0)
    {
      if (key.p_v->in[0].StrongCharge()==8 &&
	  key.p_v->in[1].StrongCharge()==8 &&
	  key.p_v->in[2].StrongCharge()==8) m_q=s_CA;
      else m_q=(key.p_v->in[0].StrongCharge()==8)?s_TR:s_CF;
      if (key.m_type==cstp::FF || key.m_type==cstp::FI) {
	if (key.p_v->in[0].StrongCharge()==8) m_q/=2.0;
      }
      else {
	if (key.m_mode==0) {
	  if (key.p_v->in[1].StrongCharge()==8) m_q/=2.0;
	}
	else {
	  if (key.p_v->in[2].StrongCharge()==8) m_q/=2.0;
	}
      }
    }

    double B0(const double &nf) const
    {
      return 11.0/6.0*s_CA-2.0/3.0*s_TR*nf;
    }

    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    double CplMax(MODEL::Running_AlphaS * as, double rsf) const;
    double Coupling(const double &scale,const int pol);
    bool AllowSpec(const ATOOLS::Flavour &fl);
    double CplFac(const double &scale) const;

    bool AllowsAlternativeCouplingUsage() const { return true; }
    void SetAlternativeUnderlyingCoupling(MODEL::Running_AlphaS *, double sf=1.0);
  };

}

using namespace CSSHOWER;
using namespace MODEL;

bool CF_QCD::SetCoupling(MODEL::Model_Base *md,
			 const double &k0sqi,const double &k0sqf,
			 const double &isfac,const double &fsfac)
{
  MODEL::Running_AlphaS * cpl = (MODEL::Running_AlphaS *)(md->GetScalarFunction("alpha_S"));

  double rsf = ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  rsf *= ToType<double>(rpa->gen.Variable("CSS_SCALE_FACTOR"));

  m_scvmode=ToType<int>(rpa->gen.Variable("CSS_SCALE_VARIATION_SCHEME"));

  // determine prefactors before calling CplMax below
  m_cplfac=((m_type/10==1)?fsfac:isfac);
  m_k0sq=(m_type/10==1)?k0sqf:k0sqi;

  m_maincplinfo = QCD_Coupling_Info(cpl, rsf);
  m_cplmax.push_back(CplMax(cpl, rsf));
  m_cplmax.push_back(0.0);
  m_maincplinfo.SetMaxCoupling(&m_cplmax);

  // invalidate alternative value and purge associated cache
  m_altcplinfo = QCD_Coupling_Info();
  m_altcplmax.clear();

  return true;
}

void CF_QCD::SetAlternativeUnderlyingCoupling(MODEL::Running_AlphaS *cpl, double sf)
{
  assert(cpl != NULL || sf == 1.0);
  if (cpl == NULL) {
    // invalidate alternative value
    m_altcplinfo = QCD_Coupling_Info();
    return;
  } else {
    m_altcplinfo = QCD_Coupling_Info(cpl, sf);
    if (m_altcplmax.find(cpl) == m_altcplmax.end()) {
      std::vector<double> altcplmax;
      altcplmax.push_back(CplMax(cpl, sf));
      altcplmax.push_back(0.0);
      m_altcplmax[cpl] = altcplmax;
    }
    m_altcplinfo.SetMaxCoupling(&(m_altcplmax[cpl]));
  }
}

double CF_QCD::CplMax(MODEL::Running_AlphaS * as, double rsf) const
{
  double minscale = Min(1.0, CplFac(m_k0sq)) * m_k0sq;
  double ct(0.);
  if (rsf > 1.) // only for f>1 cpl gets larger
    ct = -as->BoundedAlphaS(minscale) / M_PI * as->Beta0(0.) * log(rsf);
  return as->BoundedAlphaS(minscale) * (1. - ct) * m_q;
}

double CF_QCD::Coupling(const double &scale,const int pol)
{
  DEBUG_FUNC("pol="<<pol);
  if (pol!=0) return 0.0;
  QCD_Coupling_Info cplinfo = CurrentCouplingInfo();
  if (scale<0.0) return m_last = (*as)(sqr(rpa->gen.Ecms()))*m_q;
  double t(CplFac(scale)*scale), scl(CplFac(scale)*scale*cplinfo.RSF());
  double cpl=cplinfo.Coupling()->BoundedAlphaS(scl);
  msg_Debugging()<<"t="<<t<<", \\mu_R^2="<<scl<<std::endl;
  msg_Debugging()<<"as(t)="<<cplinfo.Coupling()->BoundedAlphaS(t)<<std::endl;
  if (!IsEqual(scl,t)) {
    msg_Debugging()<<"as(\\mu_R^2)="<<cpl<<std::endl;
    std::vector<double> ths(cplinfo.Coupling()->Thresholds(t,scl));
    ths.push_back((scl>t)?scl:t);
    ths.insert(ths.begin(),(scl>t)?t:scl);
    if (t<scl) std::reverse(ths.begin(),ths.end());
    msg_Debugging()<<"thresholds: "<<ths<<std::endl;
    double fac(1.),ct(0.);
    // Beta0 from One_Running_AlphaS contains extra factor 1/2
    switch (m_scvmode) {
    case 1:
      // replace as(t) -> as(t)*prod[1-as/2pi*beta(nf)*log(th[i]/th[i-1])]
      for (size_t i(1);i<ths.size();++i) {
        ct=cpl/M_PI*cplinfo.Coupling()->Beta0((ths[i]+ths[i-1])/2.0)*log(ths[i]/ths[i-1]);
        fac*=1.0-ct;
      }
      break;
    case 2:
      // replace as(t) -> as(t)*[1-sum as/2pi*beta(nf)*log(th[i]/th[i-1])]
      for (size_t i(1);i<ths.size();++i)
        ct+=cpl/M_PI*cplinfo.Coupling()->Beta0((ths[i]+ths[i-1])/2.0)*log(ths[i]/ths[i-1]);
      fac=1.-ct;
      break;
    default:
      fac=1.;
      break;
    }
    msg_Debugging()<<"ct="<<ct<<std::endl;
    if (fac<0.) {
      msg_Tracking()<<METHOD<<"(): Renormalisation term too large. Remove."
                    <<std::endl;
      fac=1.;
    }
    cpl*=fac;
    msg_Debugging()<<"as(\\mu_R^2)*(1-ct)="<<cpl<<std::endl;
  }
  cpl*=m_q;
  if (cpl>cplinfo.MaxCoupling()->front()) {
    msg_Error()<<METHOD<<"(): Value exceeds maximum at t = "
               <<sqrt(t)<<" -> \\mu_R = "<<sqrt(scl)
               <<", qmin = "<<sqrt(cplinfo.Coupling()->CutQ2())<<std::endl;
    return m_last = cplinfo.MaxCoupling()->front();
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"as weight kt = "<<sqrt(CplFac(scale))<<" * "
		 <<sqrt(scale)<<", \\alpha_s("<<sqrt(scl)<<") = "
		 <<(*cplinfo.Coupling())[scl]<<", m_q = "<<m_q<<"\n";
#endif
  return m_last = cpl;
}

bool CF_QCD::AllowSpec(const ATOOLS::Flavour &fl) 
{
  if (m_type==cstp::FF && !fl.Strong()) return true;
  if (abs(fl.StrongCharge())==3) {
    switch (m_type) {
    case cstp::FF: 
      if (abs(p_lf->FlA().StrongCharge())==3)
	return p_lf->FlA().StrongCharge()==-fl.StrongCharge();
      break;
    case cstp::FI: 
      if (abs(p_lf->FlA().StrongCharge())==3)
	return p_lf->FlA().StrongCharge()==fl.StrongCharge();
      break;
    case cstp::IF: 
      if (abs(p_lf->FlB().StrongCharge())==3)
	return p_lf->FlB().StrongCharge()==fl.StrongCharge();
      break;
    case cstp::II: 
      if (abs(p_lf->FlB().StrongCharge())==3)
	return p_lf->FlB().StrongCharge()==-fl.StrongCharge();
      break;
    case cstp::none: THROW(fatal_error,"Unknown dipole.");
    }
  }
  return fl.Strong();
}

double CF_QCD::CplFac(const double &scale) const
{
  if (m_kfmode==0) return m_cplfac;
  QCD_Coupling_Info cplinfo = CurrentCouplingInfo();
  double nf=cplinfo.Coupling()->Nf(scale);
  double kfac=exp(-(67.0-3.0*sqr(M_PI)-10.0/3.0*nf)/(33.0-2.0*nf));
  return m_cplfac*kfac;
}

namespace CSSHOWER {

DECLARE_CPL_GETTER(CF_QCD_Getter);

SF_Coupling *CF_QCD_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_QCD(args);
}

void CF_QCD_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"strong coupling";
}

}

DECLARE_GETTER(CF_QCD_Getter,"SF_QCD_Fill",
	       void,SFC_Filler_Key);

void *ATOOLS::Getter<void,SFC_Filler_Key,CF_QCD_Getter>::
operator()(const SFC_Filler_Key &key) const
{
  DEBUG_FUNC("model = "<<key.p_md->Name());
  const Vertex_Table *vtab(key.p_md->VertexTable());
  for (Vertex_Table::const_iterator
	 vlit=vtab->begin();vlit!=vtab->end();++vlit) {
    for (Vertex_List::const_iterator 
	   vit=vlit->second.begin();vit!=vlit->second.end();++vit) {
      Single_Vertex *v(*vit);
      if (v->NLegs()>3) continue;
      if (!v->PureQCD()) continue;
      msg_Debugging()<<"Add "<<v->in[0].Bar()<<" -> "<<v->in[1]<<" "<<v->in[2]<<" {\n";
      std::string atag("{"+v->in[0].Bar().IDName()+"}");
      std::string btag("{"+v->in[1].IDName()+"}");
      std::string ctag("{"+v->in[2].IDName()+"}");
      key.p_gets->push_back(new CF_QCD_Getter(atag+btag+ctag));
    }
  }
  return NULL;
}

void ATOOLS::Getter<void,SFC_Filler_Key,CF_QCD_Getter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"qcd coupling filler";
}

