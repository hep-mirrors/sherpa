#include "MCATNLO/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include <algorithm>
#include <cassert>

using namespace ATOOLS;

namespace MCATNLO {
  
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

  public:

    CF_QCD(const SF_Key &key):
      SF_Coupling(key),
      m_maincplinfo(QCD_Coupling_Info()), m_altcplinfo(QCD_Coupling_Info()),
      m_altcplmax(),
      m_q(0.0), m_k0sq(0.0)
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
    double Coupling(const double &scale,const int pol,
		    ATOOLS::Cluster_Amplitude *const sub);
    bool AllowSpec(const ATOOLS::Flavour &fl);
    double CplFac(const double &scale) const;

    bool AllowsAlternativeCouplingUsage() const { return true; }
    void SetAlternativeUnderlyingCoupling(MODEL::Running_AlphaS *, double sf=1.0);

    void ColorPoint(Parton *const p) const;
    double ColorWeight(const Color_Info &ci) const;
  };

}

using namespace MCATNLO;
using namespace MODEL;

bool CF_QCD::SetCoupling(MODEL::Model_Base *md,
			 const double &k0sqi,const double &k0sqf,
			 const double &isfac,const double &fsfac)
{
  MODEL::Running_AlphaS * cpl = (MODEL::Running_AlphaS *)md->GetScalarFunction("alpha_S");

  double rsf = ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));

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
  return as->BoundedAlphaS(minscale) * m_q;
}

double CF_QCD::Coupling(const double &scale,const int pol,
			Cluster_Amplitude *const sub)
{
  if (pol!=0) return 0.0;
  QCD_Coupling_Info cplinfo = CurrentCouplingInfo();
  double t(CplFac(scale)*scale), scl(sub?sub->MuR2():t*cplinfo.RSF());
  if (scl<cplinfo.RSF()*cplinfo.Coupling()->CutQ2()) return m_last = 0.0;
  double cpl=(*cplinfo.Coupling())(scl);
  if (sub==NULL && !IsEqual(scl,t)) {
    std::vector<double> ths(cplinfo.Coupling()->Thresholds(t,scl));
    if (scl>t) std::reverse(ths.begin(),ths.end());
    if (ths.empty() || !IsEqual(t,ths.back())) ths.push_back(t);
    if (!IsEqual(scl,ths.front())) ths.insert(ths.begin(),scl);
    for (size_t i(1);i<ths.size();++i) {
      double nf=cplinfo.Coupling()->Nf((ths[i]+ths[i-1])/2.0);
      double L=log(ths[i]/ths[i-1]), ct=cpl/(2.0*M_PI)*B0(nf)*L;
      cpl*=1.0-ct;
    }
  }
  cpl*=m_q*s_qfac;
  if (cpl>cplinfo.MaxCoupling()->front()) {
    msg_Error()<<METHOD<<"(): Value exceeds maximum at k_T = "
	       <<sqrt(scale)<<" -> q = "<<sqrt(scl)<<"."<<std::endl;
    return m_last = cplinfo.MaxCoupling()->front();
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"as weight kt = "<<(sub?1.0:sqrt(CplFac(scale)))
		 <<" * "<<(sub?sqrt(scl):sqrt(scale))<<", \\alpha_s("
		 <<sqrt(scl)<<") = "<<(*cplinfo.Coupling())[scl]
		 <<", m_q = "<<s_qfac<<" * "<<m_q<<"\n";
#endif
  return m_last = cpl;
}

bool CF_QCD::AllowSpec(const ATOOLS::Flavour &fl) 
{
  return fl.Strong();
}

double CF_QCD::CplFac(const double &scale) const
{
  if (m_kfmode==-1) return 1.0;
  if (m_kfmode==0) return m_cplfac;
  QCD_Coupling_Info cplinfo = CurrentCouplingInfo();
  double nf=cplinfo.Coupling()->Nf(scale);
  double kfac=exp(-(67.0-3.0*sqr(M_PI)-10.0/3.0*nf)/(33.0-2.0*nf));
  return m_cplfac*kfac;
}

void CF_QCD::ColorPoint(Parton *const p) const
{
  Color_Info &ci(p->Color());
  ci.m_i[0]=p->GetMEFlow(1);
  ci.m_i[1]=p->GetMEFlow(2);
  ci.m_k[0]=p->GetSpect()->GetMEFlow(1);
  ci.m_k[1]=p->GetSpect()->GetMEFlow(2);
  if (ci.m_i[0]==0 && ci.m_i[1]==0 &&
      ci.m_k[0]==0 && ci.m_k[1]==0) return;
  double disc(ran->Get());
  if (ci.m_i[0]) {
    if (ci.m_i[1]) {
      THROW(not_implemented,"1");
    }
    else {
      if (ci.m_k[0]) {
	if (ci.m_k[1]) {// q <-> g
	  THROW(not_implemented,"2");
	}
	else {// q <-> q
	  THROW(not_implemented,"3");
	}
      }
      else {// q <-> qb
	if (ci.m_i[0]==ci.m_k[1]) {
	  if (disc>5.0/14.0) {
	    ci.m_new=ci.m_i[0];
	    while (ci.m_new==ci.m_i[0])
	      ci.m_new=Min(3,(int)(3.0*ran->Get()+1.0));
	    ci.m_i[0]=ci.m_j[1]=ci.m_new;
	    ci.m_j[0]=ci.m_k[1];
	  }
	  else if (disc>1.0/14.0) {
	    ci.m_new=ci.m_j[1]=ci.m_j[0]=ci.m_i[0];
	  }
	  else {
	    ci.m_new=ci.m_i[0];
	    while (ci.m_new==ci.m_i[0])
	      ci.m_new=Min(3,(int)(3.0*ran->Get()+1.0));
	    ci.m_j[1]=ci.m_j[0]=ci.m_new;
	  }
	}
	else {
	  if (disc>4.0/7.0) {
	    ci.m_new=ci.m_j[1]=ci.m_j[0]=ci.m_k[1];
	  }
	  else if (disc>1.0/7.0) {
	    ci.m_new=ci.m_j[1]=ci.m_j[0]=ci.m_i[0];
	  }
	  else {
	    for (ci.m_new=1;ci.m_i[0]==ci.m_new ||
		   ci.m_k[1]==ci.m_new;++ci.m_new);
	    ci.m_j[1]=ci.m_j[0]=ci.m_new;
	  }
	}
	if (p_lf->FlC().Kfcode()!=kf_gluon) {
	  std::swap<int>(ci.m_i[0],ci.m_j[0]);
	  std::swap<int>(ci.m_i[1],ci.m_j[1]);
	}
      }
    }
  }
  else {
    if (ci.m_k[0]) {
      if (ci.m_k[1]) {// qb <-> g
	THROW(not_implemented,"2");
      }
      else {// qb <-> q
	if (ci.m_i[1]==ci.m_k[0]) {
	  if (disc>5.0/14.0) {
	    ci.m_new=ci.m_i[1];
	    while (ci.m_new==ci.m_i[1])
	      ci.m_new=Min(3,(int)(3.0*ran->Get()+1.0));
	    ci.m_i[1]=ci.m_j[0]=ci.m_new;
	    ci.m_j[1]=ci.m_k[0];
	  }
	  else if (disc>1.0/14.0) {
	    ci.m_new=ci.m_j[1]=ci.m_j[0]=ci.m_i[1];
	  }
	  else {
	    ci.m_new=ci.m_i[1];
	    while (ci.m_new==ci.m_i[1])
	      ci.m_new=Min(3,(int)(3.0*ran->Get()+1.0));
	    ci.m_j[1]=ci.m_j[0]=ci.m_new;
	  }
	}
	else {
	  if (disc>4.0/7.0) {
	    ci.m_new=ci.m_j[1]=ci.m_j[0]=ci.m_k[0];
	  }
	  else if (disc>1.0/7.0) {
	    ci.m_new=ci.m_j[1]=ci.m_j[0]=ci.m_i[1];
	  }
	  else {
	    for (ci.m_new=1;ci.m_i[1]==ci.m_new ||
		   ci.m_k[0]==ci.m_new;++ci.m_new);
	    ci.m_j[1]=ci.m_j[0]=ci.m_new;
	  }
	}
	if (p_lf->FlC().Kfcode()!=kf_gluon) {
	  std::swap<int>(ci.m_i[0],ci.m_j[0]);
	  std::swap<int>(ci.m_i[1],ci.m_j[1]);
	}
      }
    }
    else {// qb <-> qb
      THROW(not_implemented,"3");
    }
  }
}

double CF_QCD::ColorWeight(const Color_Info &ci) const
{
  if (ci.m_i[0]) {
    if (ci.m_i[1]) {
      if (ci.m_k[0]) {
	if (ci.m_k[1]) {// g <-> g
	  THROW(not_implemented,"1");
	}
	else {// gqb <-> q
	  if (ci.m_i[1]==ci.m_i[0]) {
	    if (ci.m_j[1]==ci.m_i[0]) {
	      if (ci.m_k[0]==ci.m_i[1]) return sqr(s_Nc-1.0)/s_Nc/s_CF;
	      return 1.0/s_CF;
	    }
	    if (ci.m_k[0]==ci.m_i[1]) return 1.0/s_CF;
	    return 1.0/s_Nc/s_CF;
	  }
	  return s_Nc/s_CF;
	}
      }
      else {// gq <-> qb
	if (ci.m_i[0]==ci.m_i[1]) {
	  if (ci.m_j[0]==ci.m_i[1]) {
	    if (ci.m_k[1]==ci.m_i[0]) return sqr(s_Nc-1.0)/s_Nc/s_CF;
	    return 1.0/s_CF;
	  }
	  if (ci.m_k[1]==ci.m_i[0]) return 1.0/s_CF;
	  return 1.0/s_Nc/s_CF;
	}
	return s_Nc/s_CF;
      }
    }
    else {
      if (ci.m_k[0]) {
	if (ci.m_k[1]) {// qg <-> g
	  THROW(not_implemented,"4");
	}
	else {// qg <-> qb
	  THROW(not_implemented,"5");
	}
      }
      else {// qg <-> q
	if (ci.m_j[0]==ci.m_j[1]) {
	  if (ci.m_i[0]==ci.m_j[1]) {
	    if (ci.m_k[1]==ci.m_j[0]) return sqr(s_Nc-1.0)/s_Nc/s_CF;
	    return 1.0/s_CF;
	  }
	  if (ci.m_k[1]==ci.m_j[0]) return 1.0/s_CF;
	  return 1.0/s_Nc/s_CF;
	}
	return s_Nc/s_CF;
      }
    }
  }
  else {
    if (ci.m_k[1]) {
      if (ci.m_k[0]) {// qbg <-> g
	THROW(not_implemented,"6");
      }
      else {// qbg <-> qb
	THROW(not_implemented,"7");
      }
    }
    else {// qbg <-> q
      if (ci.m_j[0]==ci.m_j[1]) {
	if (ci.m_i[1]==ci.m_j[0]) {
	  if (ci.m_k[0]==ci.m_j[1]) return sqr(s_Nc-1.0)/s_Nc/s_CF;
	  return 1.0/s_CF;
	}
	if (ci.m_k[0]==ci.m_j[1]) return 1.0/s_CF;
	return 1.0/s_Nc/s_CF;
      }
      return s_Nc/s_CF;
    }
  }
  return 1.0;
}

namespace MCATNLO {

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
  if (!Flavour(kf_gluon).IsOn()) return NULL;
  std::string gtag("{"+Flavour(kf_gluon).IDName()+"}");
  key.p_gets->push_back(new CF_QCD_Getter(gtag+gtag+gtag));
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)i);
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    key.p_gets->push_back(new CF_QCD_Getter(gtag+qtag+qbtag));
    key.p_gets->push_back(new CF_QCD_Getter(qbtag+qbtag+gtag));
    key.p_gets->push_back(new CF_QCD_Getter(qtag+qtag+gtag));
  }
  if (MODEL::s_model->Name().find("MSSM")==std::string::npos) return NULL;
  else {
    msg_Out()<<METHOD<<"(): MC@NLO does not shower MSSM particles yet.\n";
    return NULL;
  }
  std::string sgtag("{"+Flavour(kf_Gluino).IDName()+"}");
  key.p_gets->push_back(new CF_QCD_Getter(sgtag+sgtag+gtag));
  key.p_gets->push_back(new CF_QCD_Getter(sgtag+gtag+sgtag));
  key.p_gets->push_back(new CF_QCD_Getter(gtag+sgtag+sgtag));
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)(1000000+i));
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    key.p_gets->push_back(new CF_QCD_Getter(gtag+qtag+qbtag));
    key.p_gets->push_back(new CF_QCD_Getter(qbtag+qbtag+gtag));
    key.p_gets->push_back(new CF_QCD_Getter(qtag+qtag+gtag));
  }
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)(2000000+i));
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    key.p_gets->push_back(new CF_QCD_Getter(gtag+qtag+qbtag));
    key.p_gets->push_back(new CF_QCD_Getter(qbtag+qbtag+gtag));
    key.p_gets->push_back(new CF_QCD_Getter(qtag+qtag+gtag));
  }
  return NULL;
}

void ATOOLS::Getter<void,SFC_Filler_Key,CF_QCD_Getter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"qcd coupling filler";
}
