#include "MCATNLO/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include <algorithm>

using namespace ATOOLS;

namespace MCATNLO {
  
  const double s_Nc = 3.;
  const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  const double s_CA = s_Nc;
  const double s_TR = 1./2.;

  class CF_QCD: public SF_Coupling {
  protected:

    //! Main underlying coupling (set by SetCoupling)
    MODEL::Running_AlphaS *p_cpl;

    //! (Temporary) alternative coupling (set by SetAlternativeUnderlyingCoupling)
    MODEL::One_Running_AlphaS *p_altcpl;

    //! Alternative renormalisation scale factor (set by SetAlternativeUnderlyingCoupling)
    double m_altrsf;

    //! Buffer of max alphas values to avoid re-calculations
    std::map<MODEL::One_Running_AlphaS *, double> m_altcplmax;

    double m_q, m_rsf, m_k0sq;

    double B0(const double &nf) const
    {
      return 11.0/6.0*s_CA-2.0/3.0*s_TR*nf;
    }

  public:

    inline CF_QCD(const SF_Key &key):
      SF_Coupling(key), p_altcpl(NULL), m_altrsf(1.0), m_altcplmax(), m_k0sq(0.0)
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

    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    template<class T>
    double CplMax(T * as) const;
    double Coupling(const double &scale,const int pol,
		    ATOOLS::Cluster_Amplitude *const sub);
    bool AllowSpec(const ATOOLS::Flavour &fl);

    double CplFac(const double &scale) const;

    bool AllowsAlternativeCouplingUsage() const { return true; }
    void SetAlternativeUnderlyingCoupling(void *, double sf=1.0);

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
  p_cpl=(MODEL::Running_AlphaS*)md->GetScalarFunction("alpha_S");
  SetAlternativeUnderlyingCoupling(NULL);
  m_altcplmax.clear(); // buffered values are not valid anymore
  m_rsf=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  m_cplfac=((m_type/10==1)?fsfac:isfac);
  m_k0sq=(m_type/10==1)?k0sqf:k0sqi;
  m_cplmax.push_back(CplMax(p_cpl));
  m_cplmax.push_back(0.0);
  return true;
}

void CF_QCD::SetAlternativeUnderlyingCoupling(void *cpl, double sf)
{
  m_altrsf = sf;
  if (cpl == NULL) {
    p_altcpl = NULL;
    return;
  } else {
    p_altcpl = static_cast<MODEL::One_Running_AlphaS *>(cpl);
    if (m_altcplmax.find(p_altcpl) == m_altcplmax.end()) {
      m_altcplmax[p_altcpl] = CplMax(p_altcpl);
    }
  }
}

template<class T>
double CF_QCD::CplMax(T * as) const
{
  double minscale = Min(1.0, CplFac(m_k0sq)) * m_k0sq;
  return (*as)(Max(as->CutQ2(), minscale))*m_q;
}

double CF_QCD::Coupling(const double &scale,const int pol,
			Cluster_Amplitude *const sub)
{
  if (pol!=0) return 0.0; // we do not update m_last when polarized
  One_Running_AlphaS * const as = (p_altcpl) ? p_altcpl : p_cpl->GetAs();
  const double rsf = (p_altcpl) ? m_altrsf : m_rsf;

  double t(CplFac(scale)*scale), scl(sub?sub->MuR2():t*rsf);
  if (scl<rsf*as->CutQ2()) return m_last = 0.0;
  double cpl=(*as)(scl);
  if (sub==NULL && !IsEqual(scl,t)) {
    std::vector<double> ths(as->Thresholds(t,scl));
    if (scl>t) std::reverse(ths.begin(),ths.end());
    if (ths.empty() || !IsEqual(t,ths.back())) ths.push_back(t);
    if (!IsEqual(scl,ths.front())) ths.insert(ths.begin(),scl);
    for (size_t i(1);i<ths.size();++i) {
      double nf=as->Nf((ths[i]+ths[i-1])/2.0);
      double L=log(ths[i]/ths[i-1]), ct=cpl/(2.0*M_PI)*B0(nf)*L;
      cpl*=1.0-ct;
    }
  }
  cpl*=m_q*s_qfac;
  const double cplmax = (p_altcpl) ? s_qfac*m_altcplmax[p_altcpl] : s_qfac*m_cplmax.front();
  if (cpl>cplmax) {
    msg_Error()<<METHOD<<"(): Value exceeds maximum at k_T = "
	       <<sqrt(scale)<<" -> q = "<<sqrt(scl)<<"."<<std::endl;
    return m_last = cplmax;
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"as weight kt = "<<(sub?1.0:sqrt(CplFac(scale)))
		 <<" * "<<(sub?sqrt(scl):sqrt(scale))<<", \\alpha_s("
		 <<sqrt(scl)<<") = "<<(*as)[scl]
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
  One_Running_AlphaS * const as = (p_altcpl) ? p_altcpl : p_cpl->GetAs();
  double nf=as->Nf(scale);
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
