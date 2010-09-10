#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Exception.H"

#define CHECK__x
// #define METS__reject_unordered
// #define CHECK__stepwise

namespace PHASIC {

  struct CS_Params {
    size_t m_i, m_j, m_k;
    ATOOLS::Flavour m_fl;
    double m_kt2, m_op2, m_z, m_y;
    ATOOLS::Vec4D m_pijt, m_pkt;
    CS_Params(const size_t &i,const size_t &j,
	      const size_t &k,const ATOOLS::Flavour &fl):
      m_i(i),m_j(j), m_k(k), m_fl(fl),
      m_kt2(-1.0), m_op2(0.0), m_z(0.0), m_y(0.0) {}
    bool operator<(const CS_Params &ck) const
    { 
      if (m_i<ck.m_i) return true;
      if (m_i>ck.m_i) return false;
      if (m_j<ck.m_j) return true;
      if (m_j>ck.m_j) return false;
      if (m_k<ck.m_k) return true;
      if (m_k>ck.m_k) return false;
      return m_fl<ck.m_fl;
    }
    void SetParams(const double &kt2,const double &z,const double &y,
		   const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &pkt)
    { m_kt2=kt2, m_z=z; m_y=y; m_pijt=pijt; m_pkt=pkt; }
  };// end of struct CS_Params

  class METS_Scale_Setter: public Scale_Setter_Base {
  private:

    std::string m_muf2tag, m_mur2tag;

    ATOOLS::Algebra_Interpreter m_muf2calc, m_mur2calc;

    Tag_Setter m_muf2tagset, m_mur2tagset;

    ATOOLS::Vec4D_Vector   m_p;
    ATOOLS::Flavour_Vector m_f;

    SP(Color_Integrator) p_ci;

    size_t m_cnt, m_rej, m_mode;
    double m_lfrac, m_aqed;

    static double s_eps, s_kt2max;

    bool CheckX(const ATOOLS::Vec4D &p,const size_t &isid) const;

    double CoreScale(ATOOLS::Cluster_Amplitude *const ampl);

    bool CheckColors(const ATOOLS::Cluster_Leg *li,
		     const ATOOLS::Cluster_Leg *lj,
		     const ATOOLS::Cluster_Leg *lk,
		     const ATOOLS::Flavour &mo) const;
    ATOOLS::ColorID CombineColors(const ATOOLS::Cluster_Leg *li,
				  const ATOOLS::Cluster_Leg *lj,
				  const ATOOLS::Cluster_Leg *lk,
				  const ATOOLS::Flavour &mo) const;

    double Lam(const double &s,
	       const double &sb,const double &sc) const;

    void KT2(const ATOOLS::Cluster_Leg *li,
	     const ATOOLS::Cluster_Leg *lj,
	     const ATOOLS::Cluster_Leg *lk,CS_Params &cs) const;

    bool Combine(ATOOLS::Cluster_Amplitude &ampl,
		 int i,int j,int k,const CS_Params &cs) const;

    double SetScales(const double &scale);

    double CalculateStrict(const ATOOLS::Vec4D_Vector &momenta);

  public:

    METS_Scale_Setter(const Scale_Setter_Arguments &args,
		      const int mode=1);

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);
    double CalculateScale2(const std::vector<ATOOLS::Vec4D> &p);

    ATOOLS::Vec4D Momentum(const size_t &i) const;

    void SetScale(const std::string &mu2tag,Tag_Setter &mu2tagset,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class METS_Scale_Setter

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Loose_METS_Scale_Setter_Getter,"METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *Loose_METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args,0);
}

void Loose_METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"loose mets scale scheme\n";
}

DECLARE_GETTER(METS_Scale_Setter_Getter,"SEMI_STRICT_METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args,1);
}

void METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mets scale scheme\n";
}

DECLARE_GETTER(Strict_METS_Scale_Setter_Getter,"STRICT_METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *Strict_METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args,2);
}

void Strict_METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"strict mets scale scheme\n";
}

double METS_Scale_Setter::s_eps=1.0e-6;
double METS_Scale_Setter::s_kt2max=
       sqrt(std::numeric_limits<double>::max());

METS_Scale_Setter::METS_Scale_Setter
(const Scale_Setter_Arguments &args,const int mode):
  Scale_Setter_Base(args), m_muf2tagset(this), m_mur2tagset(this),
  m_cnt(0), m_rej(0), m_mode(mode), m_lfrac(0.0)
{
  m_p.resize(4);
  size_t pos(args.m_scale.find('{'));
  std::string mur2tag("MU_R2"), muf2tag("MU_F2");
  if (pos!=std::string::npos) {
    muf2tag=args.m_scale.substr(pos+1);
    pos=muf2tag.rfind('}');
    if (pos==std::string::npos)
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    muf2tag=muf2tag.substr(0,pos);
    pos=muf2tag.find("}{");
    if (pos==std::string::npos) {
      mur2tag=muf2tag;
    }
    else {
      mur2tag=muf2tag.substr(pos+2);
      muf2tag=muf2tag.substr(0,pos);
    }
  }
  SetScale(muf2tag,m_muf2tagset,m_muf2calc);
  SetScale(mur2tag,m_mur2tagset,m_mur2calc);
  m_scale.resize(p_proc->NOut());
  SetCouplings();
  m_f=p_proc->Flavours();
  m_aqed=(*MODEL::aqed)(sqr(Flavour(kf_Z).Mass()));
  for (size_t i(0);i<p_proc->NIn();++i) m_f[i]=m_f[i].Bar();
}

Vec4D METS_Scale_Setter::Momentum(const size_t &i) const
{
  if (i>m_p.size()) THROW(fatal_error,"Momentum index too large");
  return m_p[i];
}

double METS_Scale_Setter::CalculateStrict(const Vec4D_Vector &momenta)
{
  p_caller->Integrator()->SetMomenta(momenta);
  p_caller->Generator()->SetClusterDefinitions
    (p_caller->Shower()->GetClusterDefinitions());
  Cluster_Amplitude *ampl
    (p_caller->Generator()->ClusterConfiguration(p_caller));
  if (ampl==NULL) {
    ++m_rej;
    double frac(m_rej/(double)m_cnt);
    if (frac>1.25*m_lfrac && m_cnt>5000) {
      m_lfrac=frac;
      msg_Error()<<METHOD<<"(): No CSS history for '"
		 <<p_proc->Name()<<"' in >"
		 <<(int(m_lfrac*10000)/100.0)
		 <<"% of calls. Set \\hat{s}."<<std::endl;
    }
    return SetScales((m_p[0]+m_p[1]).Abs2());
  }
  Cluster_Amplitude *rampl(ampl);
  while (ampl->Next()) ampl=ampl->Next();
  double kt2max(ampl->KT2QCD());
  msg_Debugging()<<"Core = "<<*ampl<<"\n";
  m_p.resize(ampl->Legs().size());
  for (size_t i(0);i<m_p.size();++i)
    m_p[i]=ampl->Leg(i)->Mom();
  rampl->Delete();
  return SetScales(kt2max);
}

double METS_Scale_Setter::CalculateScale2(const Vec4D_Vector &momenta) 
{
  if (m_mode==2 || (m_mode==1 && !p_caller->LookUp())) {
    p_caller->Integrator()->SetMomenta(momenta);
    p_caller->Integrator()->SwapInOrder();
    double muf2(CalculateScale(p_caller->Integrator()->Momenta()));
    p_caller->Integrator()->RestoreInOrder();
    return muf2;
  }
  p_cpls->Calculate();
  return m_scale[stp::fac];
}

double METS_Scale_Setter::CalculateScale(const Vec4D_Vector &momenta) 
{
  ++m_cnt;
  m_p=momenta;
  p_ci=p_proc->Integrator()->ColorIntegrator();
  for (size_t i(0);i<p_proc->NIn();++i) m_p[i]=-m_p[i];
  if (m_mode==2 || (m_mode==1 && !p_caller->LookUp())) {
    m_scale2=p_proc->Integrator()->ISR()->On() && m_f[0]!=m_f[1];
    return CalculateStrict(momenta);
  }
  DEBUG_FUNC(p_proc->Name());
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  ampl->SetNIn(p_proc->NIn());
  if (p_ci==NULL) {
    for (size_t i(0);i<m_p.size();++i) ampl->CreateLeg(m_p[i],m_f[i]);
  }
  else {
    Int_Vector ci(p_ci->I()), cj(p_ci->J());
    for (size_t i(0);i<m_p.size();++i) 
      ampl->CreateLeg(m_p[i],m_f[i],ColorID(ci[i],cj[i]));
  }
  Single_Process *proc(p_proc->Get<Single_Process>());
  std::vector<std::set<CS_Params> > alltrials(ampl->Legs().size()-4);
  std::vector<double> ops(ampl->Legs().size()-3,0.0);
  double kt2core(ampl->Legs().size()>4?0.0:CoreScale(ampl));
  while (ampl->Legs().size()>4) {
    msg_Debugging()<<"Actual = "<<*ampl<<"\n";
    std::set<CS_Params> &trials(alltrials[ampl->Legs().size()-5]);
    size_t iw(0), jw(0), kw(0);
    CS_Params ckw(0,0,0,kf_none);
    for (size_t i(0);i<ampl->Legs().size();++i) {
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(Max((size_t)2,i+1));j<ampl->Legs().size();++j) {
	Cluster_Leg *lj(ampl->Leg(j));
	if (!proc->Combinable(li->Id(),lj->Id())) continue;
	const Flavour_Vector &cf(proc->CombinedFlavour(li->Id()+lj->Id()));
	for (size_t k(0);k<ampl->Legs().size();++k) {
	  Cluster_Leg *lk(ampl->Leg(k));
	  if (k!=i && k!=j) {
	    for (size_t f(0);f<cf.size();++f) {
	      CS_Params cs(li->Id(),lj->Id(),lk->Id(),cf[f]);
	      if (trials.find(cs)!=trials.end()) continue;
	      if (!CheckColors(li,lj,lk,cf[f])) {
		msg_Debugging()<<"Veto colours: "<<cf[f]<<" = "
			       <<ID(cs.m_i)<<" & "<<ID(cs.m_j)
			       <<" <-> "<<ID(cs.m_k)<<"\n";
		trials.insert(cs);
		continue;
	      }
	      KT2(li,lj,lk,cs);
	      if (cs.m_kt2==-1.0) {
		msg_Debugging()<<"Veto kinematics: "<<cf[f]<<" = "
			       <<ID(cs.m_i)<<" & "<<ID(cs.m_j)
			       <<" <-> "<<ID(cs.m_k)<<"\n";
		trials.insert(cs);
		continue;
	      }
	      cs.m_op2=1.0/cs.m_kt2;
	      if (cf[f].Strong() &&
		  li->Flav().Strong() &&
		  lj->Flav().Strong()) {
		// strong clustering, reweight with as
		cs.m_op2*=(*MODEL::as)(cs.m_kt2);
	      }
	      else {
		// crude: reweight with em coupling
		cs.m_op2*=m_aqed;
		if (cf[f].IsPhoton() ||
		    li->Flav().IsPhoton() ||
		    lj->Flav().IsPhoton()) {
		  // if photon, reweight with charge
		  if (li->Flav().IsPhoton() ||
		      lj->Flav().IsPhoton()) {
		    cs.m_op2*=dabs(cf[f].Charge());}
		  else {
		    cs.m_op2*=dabs(li->Flav().Charge());
		  }
		}
		else if (cf[f].Mass()) {
		  if (cf[f].Width()) {
		    // if resonance, reweight with breit-wigner
		    double s((li->Mom()+lj->Mom()).Abs2());
		    double m2(sqr(cf[f].Mass()));
		    cs.m_op2*=cs.m_kt2/
		      sqrt(sqr(s-m2)+m2*sqr(cf[f].Width()));
		  }
		  else {
		    // if non-resonant, reweight with massive prop
		    cs.m_op2*=cs.m_kt2/(cs.m_kt2+sqr(cf[f].Mass()));
		  }
		}
	      }
	      if (cs.m_op2>ckw.m_op2) {
		ckw=cs;
		iw=i;
		jw=j;
		kw=k;
	      }
	    }
	  }
	}
      }
    }
    trials.insert(ckw);
    if (iw==0 && jw==0 && kw==0) {
      if (ampl->Prev()==NULL) {
	++m_rej;
	double frac(m_rej/(double)m_cnt);
	if (frac>1.25*m_lfrac && m_cnt>5000) {
	  m_lfrac=frac;
	  msg_Error()<<METHOD<<"(): No CSS history for '"
		     <<p_proc->Name()<<"' in >"
		     <<(int(m_lfrac*10000)/100.0)
		     <<"% of calls. Set \\hat{s}."<<std::endl;
	}
	ampl->Delete();
	return SetScales((m_p[0]+m_p[1]).Abs2());
      }
      ampl=ampl->Prev();
      ampl->DeleteNext();
      trials.clear();
      continue;
    }
    msg_Debugging()<<"Cluster "<<ckw.m_fl<<" "
		   <<ID(ckw.m_i)<<" & "<<ID(ckw.m_j)
		   <<" <-> "<<ID(ckw.m_k)
		   <<" => "<<sqrt(ckw.m_kt2)
		   <<" ("<<sqrt(ckw.m_op2)<<") <-> "
		   <<sqrt(ops[ampl->Legs().size()-4])<<"\n";
    if (ops[ampl->Legs().size()-4]>ckw.m_kt2) {
      msg_Debugging()<<"unordered configuration\n";
#ifdef METS__reject_unordered
      continue;
#endif
    }
    ampl->SetKT2QCD(ckw.m_kt2);
    ampl=ampl->InitNext();
    ampl->CopyFrom(ampl->Prev());
    if (!Combine(*ampl,iw,jw,kw,ckw)) {
      ampl=ampl->Prev();
      ampl->DeleteNext();
      continue;
    }
    ops[ampl->Legs().size()-4]=ckw.m_kt2;
#ifdef CHECK__stepwise
    Vec4D psum;
    for (size_t i(0);i<ampl->Legs().size();++i) {
      psum+=ampl->Leg(i)->Mom();
    }
    if (!IsEqual(psum,Vec4D(),1.0e-3)) {
      msg_Error()<<METHOD<<"(): Momentum not conserved\n  in process '"
		 <<p_proc->Name()<<"'\n  \\sum p = "<<psum
		 <<" in step \n"<<*ampl->Prev()<<"\n"<<*ampl<<std::endl;
    }
#endif
    if (ampl->Legs().size()==4) {
      if (!proc->Combinable(ampl->Leg(0)->Id(),ampl->Leg(1)->Id()) &&
	  !proc->Combinable(ampl->Leg(0)->Id(),ampl->Leg(2)->Id()) &&
	  !proc->Combinable(ampl->Leg(0)->Id(),ampl->Leg(3)->Id())) {
	ampl=ampl->Prev();
	ampl->DeleteNext();
	continue;
      }
      kt2core=CoreScale(ampl);
      if (ops[ampl->Legs().size()-4]>kt2core) {
	msg_Debugging()<<"unordered configuration (core)\n";
#ifdef METS__reject_unordered
	ampl=ampl->Prev();
	ampl->DeleteNext();
	continue;
#endif
      }
    }
  }
  size_t idx(2);
  while (ampl->Prev()) {
    ampl=ampl->Prev();
    m_scale[idx++]=ampl->KT2QCD();
  }
  ampl->Delete();
  return SetScales(kt2core);
}

double METS_Scale_Setter::CoreScale(Cluster_Amplitude *const ampl)
{
  msg_Debugging()<<"Core = "<<*ampl<<"\n";
  m_p.resize(ampl->Legs().size());
  Vec4D psum;
  int res(0), qcd(0), csum[4]={0,0,0,0};
  size_t cid[4]={ampl->Leg(0)->Id(),ampl->Leg(1)->Id(),
		 ampl->Leg(2)->Id(),ampl->Leg(3)->Id()};
  ColorID c[4]={ampl->Leg(0)->Col(),ampl->Leg(1)->Col(),
		ampl->Leg(2)->Col(),ampl->Leg(3)->Col()};
  for (size_t i(0);i<m_p.size();++i) {
    Cluster_Leg *li(ampl->Leg(i));
    psum+=m_p[i]=li->Mom();
    ++csum[c[i].m_i];
    --csum[c[i].m_j];
    if (c[i].m_i>0 || c[i].m_j>0) qcd+=1<<i;
    if (ampl->Leg(i)->Flav().Strong() ||
	ampl->Leg(i)->Flav().Resummed()) res+=1<<i;
  }
  if (!IsEqual(psum,Vec4D(),1.0e-3)) {
    msg_Error()<<METHOD<<"(): Momentum not conserved.\n"
	       <<"\\sum p = "<<psum<<" in\n"<<*ampl<<std::endl;
  }
  if (csum[1]!=0 || csum[2]!=0 || csum[3]!=0) {
    msg_Error()<<METHOD<<"(): Colour not conserved. "<<*ampl<<std::endl;
    abort();
  }
  bool pureres(false);
  Single_Process *proc(p_proc->Get<Single_Process>());
  if ((res&7)==7 || (res&11)==11) {
    for (size_t j(1);j<4;++j) {
      if (proc->Combinable(cid[0],cid[j])) {
	const Flavour_Vector &cf(proc->CombinedFlavour(cid[0]+cid[j]));
	for (size_t i(0);i<cf.size();++i)
	  if (cf[i].Resummed() || cf[i].Strong()) {
	    pureres=true;
	    break;
	  }
      }
      if (pureres) break;
    }
  }
  double kt2cmin(s_kt2max);
  if (pureres) {
    kt2cmin=Max(m_p[2].MPerp2(),m_p[3].MPerp2());
  }
  else {
    // s-channel
    if (proc->Combinable(cid[0],cid[1])) {
      if (p_ci==NULL || qcd==0 ||
	  (c[0].m_i>0 && c[0].m_i==c[1].m_j) ||
	  (c[0].m_j>0 && c[0].m_j==c[1].m_i) ||
	  (c[2].m_i>0 && c[2].m_i==c[3].m_j) ||
	  (c[2].m_j>0 && c[2].m_j==c[3].m_i)) {
	kt2cmin=Min(kt2cmin,
		    Max((m_p[0]+m_p[1]).Abs2(),
			(m_p[2]+m_p[3]).Abs2()));
      }
    }
    // t-channel
    if (proc->Combinable(cid[0],cid[2])) {
      if (p_ci==NULL || qcd==0 ||
	  (c[0].m_i>0 && c[0].m_i==c[2].m_j) ||
	  (c[0].m_j>0 && c[0].m_j==c[2].m_i) ||
	  (c[1].m_i>0 && c[1].m_i==c[3].m_j) ||
	  (c[1].m_j>0 && c[1].m_j==c[3].m_i)) {
	kt2cmin=Min(kt2cmin,
		    Max(dabs((m_p[0]+m_p[2]).Abs2()),
			dabs((m_p[1]+m_p[3]).Abs2())));
      }
    }
    // u-channel
    if (proc->Combinable(cid[0],cid[3])) {
      if (p_ci==NULL || qcd==0 ||
	  (c[0].m_i>0 && c[0].m_i==c[3].m_j) ||
	  (c[0].m_j>0 && c[0].m_j==c[3].m_i) ||
	  (c[1].m_i>0 && c[1].m_i==c[2].m_j) ||
	  (c[1].m_j>0 && c[1].m_j==c[2].m_i)) {
	kt2cmin=Min(kt2cmin,
		    Max(dabs((m_p[0]+m_p[3]).Abs2()),
			dabs((m_p[1]+m_p[2]).Abs2())));
      }
    }
  }
  return kt2cmin;
}

double METS_Scale_Setter::SetScales(const double &scale)
{
  m_scale[stp::ren]=m_scale[stp::fac]=scale;
  msg_Debugging()<<"QCD scale = "<<sqrt(m_scale[stp::ren])<<"\n";
  m_scale[stp::ren]=m_mur2calc.Calculate()->Get<double>();
  m_scale[stp::fac]=m_muf2calc.Calculate()->Get<double>();
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n";
  for (size_t i(2);i<m_scale.size();++i)
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_scale[i])<<"\n";
  msg_Debugging()<<"} <- "<<p_proc->Name()<<"\n";
  p_cpls->Calculate();
  return m_scale[stp::fac];
}

void METS_Scale_Setter::SetScale
(const std::string &mu2tag,Tag_Setter &mu2tagset,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Name()<<"' {\n";
  msg_Indent();
  mu2tagset.SetCalculator(&mu2calc);
  mu2calc.SetTagReplacer(&mu2tagset);
  mu2tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}

double METS_Scale_Setter::Lam
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

void METS_Scale_Setter::KT2
(const Cluster_Leg *li,const Cluster_Leg *lj,
 const Cluster_Leg *lk,CS_Params &cs) const
{
  if ((li->Id()&3)<(lj->Id()&3)) std::swap<const Cluster_Leg*>(li,lj);
  Vec4D pi(li->Mom()), pj(lj->Mom()), pk(lk->Mom()), Q(pi+pj+pk);
  double Q2=Q.Abs2(), mij2=sqr(cs.m_fl.Mass()), mk2=sqr(lk->Flav().Mass());
  double lrat=Lam(Q2,mij2,mk2)/Lam(Q2,(pi+pj).Abs2(),mk2);
  Vec4D pkt(sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-mij2)/(2.*Q2)*Q), pijt(Q-pkt);
  if (lrat<0.0 || Sign(pkt[0])!=Sign(pk[0]) || Sign(pijt[0])!=Sign(pi[0]) ||
      IsZero(pkt[0],1.0e-6) || IsZero(pijt[0],1.0e-6)) return;
  if ((li->Id()&3)==0) {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
 	double pipj=pi*pj, pipk=pi*pk, pjpk=pj*pk, Q2=Q*Q;
  	double mi2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
 	double yijk=pipj/(pipj+pipk+pjpk), zi=pipk/(pipk+pjpk);
  	double kt2=(Q2-mi2-mj2-mk2)*yijk*zi*(1.-zi)
  	  -(1.0-zi)*(1.0-zi)*mi2-zi*zi*mj2;
 	cs.SetParams(kt2,zi,yijk,pijt,pkt);
      }
      else {
 	double pipj=pi*pj, pipa=pi*pk, pjpa=pj*pk;
  	double mi2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
 	double xija=(pipa+pjpa+pipj)/(pipa+pjpa), zi=pipa/(pipa+pjpa);
  	double kt2=-2.0*(pipa+pjpa)*(1.0-xija)*zi*(1.0-zi)
  	  -sqr(1.0-zi)*mi2-zi*zi*mj2;
 	cs.SetParams(kt2,zi,1.0-xija,pijt,pkt);
      }
    }
  }
  else {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
 	double pjpa=pj*pi, pkpa=pk*pi, pjpk=pj*pk;
  	double ma2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
 	double xjka=(pjpa+pkpa+pjpk)/(pjpa+pkpa), uj=pjpa/(pjpa+pkpa);
  	double kt2=-2.*(pjpa+pkpa)*(1.-xjka)*uj-mj2-sqr(1.0-xjka)*ma2;
 	cs.SetParams(kt2,xjka,uj,pijt,pkt);
      }
      else {
 	double papb=pi*pk, pjpa=pj*pi, pjpb=pj*pk;
  	double mj2=sqr(lj->Flav().Mass()), ma2=sqr(li->Flav().Mass());
 	double xjab=(papb+pjpa+pjpb)/papb, vj=-pjpa/papb;
  	double kt2=2.0*papb*vj*(1.0-xjab)-mj2-sqr(1.0-xjab)*ma2;
 	cs.SetParams(kt2,xjab,vj,pijt,pkt);
      }
    }
  }
}
  
bool METS_Scale_Setter::Combine(Cluster_Amplitude &ampl,int i,int j,int k,
				const CS_Params &cs) const
{
  if (i>j) std::swap<int>(i,j);
  Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j)), *lk(ampl.Leg(k));
  li->SetCol(CombineColors(li,lj,lk,cs.m_fl));
  li->SetFlav(cs.m_fl);
  bool ii(i<2 && k<2);
  if (ii) li->SetMom(li->Mom()+lj->Mom());
  else {
    li->SetMom(cs.m_pijt);
    lk->SetMom(cs.m_pkt);
  }
  if (i<2 || k<2) {
    Cluster_Leg *la(ampl.Leg(i<2?i:k)), *lb(ampl.Leg(1-(i<2?i:k)));
    ZAlign lt(-la->Mom(),-lb->Mom(),ii?li->Mom().Abs2():
	      sqr(la->Flav().Mass()),sqr(lb->Flav().Mass()));
    for (size_t m(0);m<ampl.Legs().size();++m) {
      ampl.Leg(m)->SetMom(lt.Align(ampl.Leg(m)->Mom()));
      ampl.Leg(m)->SetK(0);
    }
    if (ii) {
      li->SetMom(lt.Align(cs.m_pijt));
      lk->SetMom(lt.Align(cs.m_pkt));
    }
#ifdef CHECK__x
    if (!CheckX(la->Mom(),la->Id()&3)) return false;
    if (!CheckX(lb->Mom(),lb->Id()&3)) return false;
#endif
  }
  li->SetId(li->Id()+lj->Id());
  li->SetK(lk->Id());
  std::vector<Cluster_Leg*>::iterator lit(ampl.Legs().begin());
  for (int l(0);l<j;++l) ++lit;
  (*lit)->Delete();
  ampl.Legs().erase(lit);
  return true;
}

bool METS_Scale_Setter::CheckX
(const ATOOLS::Vec4D &p,const size_t &isid) const
{
  double x=0.0;
  if (isid==1) x=-p.PPlus()/rpa.gen.PBeam(0).PPlus();
  else if (isid==2) x=-p.PMinus()/rpa.gen.PBeam(1).PMinus();
  else THROW(fatal_error,"Invalid index");
  return x>1.0e-6 && x<1.0;
}

bool METS_Scale_Setter::CheckColors
(const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
 const ATOOLS::Cluster_Leg *lk,const ATOOLS::Flavour &mo) const
{
  if (mo.StrongCharge()==8) {
    if (!lk->Flav().Strong()) return false;
  }
  else if (mo.Strong()) {
    if (!(lk->Flav().StrongCharge()==8 ||
	  lk->Flav().StrongCharge()==-mo.StrongCharge())) return false;
  }
  else {
    if (lk->Flav().StrongCharge()==8) return false;
    if (li->Col().m_i==-1 && lj->Col().m_i==-1 &&
	lk->Col().m_i==-1) return true;
    ColorID ci(li->Col()), cj(lj->Col());
    if (ci.m_i==cj.m_j && ci.m_j==0 && cj.m_i==0) return true;
    if (ci.m_j==cj.m_i && ci.m_i==0 && cj.m_j==0) return true;
    return false;
  }
  if (li->Col().m_i==-1 && lj->Col().m_i==-1 &&
      lk->Col().m_i==-1) return true;
  ColorID ci(li->Col()), cj(lj->Col()), ck(lk->Col());
  if (ci.m_i<0 && cj.m_i<0 && ck.m_i<0) return true;
  if (li->Flav().StrongCharge()==3) {
    if (lj->Flav().StrongCharge()==-3) {
      if (lk->Flav().StrongCharge()==0) return true;
      if (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	  (ci.m_i==cj.m_j && (ck.m_i>0 || ck.m_j>0))) return true;
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (lk->Flav().StrongCharge()==0) return false;
      if (ci.m_i==cj.m_j && 
	  (cj.m_i==ck.m_j || ck.Singlet())) return true;
      if ((ci.m_i==ck.m_j || ck.Singlet()) && 
	  cj.Singlet()) return true;
    }
    else {
      if (lk->Flav().StrongCharge()==8) return false;
      return true;
    }
  }
  else if (li->Flav().StrongCharge()==-3) {
    if (lj->Flav().StrongCharge()==3) {
      if (lk->Flav().StrongCharge()==0) return true;
      if (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	  (ci.m_j==cj.m_i && (ck.m_i>0 || ck.m_j>0))) return true;
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (lk->Flav().StrongCharge()==0) return false;
      if (ci.m_j==cj.m_i && 
	  (cj.m_j==ck.m_i || ck.Singlet())) return true;
      if ((ci.m_j==ck.m_i || ck.Singlet()) && 
	  cj.Singlet()) return true;
    }
    else {
      if (lk->Flav().StrongCharge()==8) return false;
      return true;
    }
  }
  else if (li->Flav().StrongCharge()==8) {
    if (lk->Flav().StrongCharge()==0) return false;
    if (lj->Flav().StrongCharge()==8) {
      if (ci.m_i==cj.m_j && 
	  (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	   (ci.m_j==cj.m_i && lk->Flav().StrongCharge()!=8))) 
	return true;
      if (ci.m_j==cj.m_i && 
	  (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	   (ci.m_i==cj.m_j && lk->Flav().StrongCharge()!=8)))
	return true;
    }
    else if (lj->Flav().StrongCharge()==3) {
      if (ci.m_j==cj.m_i &&
	  (ci.m_i==ck.m_j || ck.Singlet())) return true;
      if ((cj.m_i==ck.m_j || ck.Singlet()) &&
	  ci.Singlet()) return true;
    }
    else if (lj->Flav().StrongCharge()==-3) {
      if (ci.m_i==cj.m_j &&
	  (ci.m_j==ck.m_i || ck.Singlet())) return true;
      if ((cj.m_j==ck.m_i || ck.Singlet()) &&
	  ci.Singlet()) return true;
    }
    else {
      return false;
    }
  }
  else {
    if (lj->Flav().StrongCharge()==8 ||
	lk->Flav().StrongCharge()==8) {
      return false;
    }
    return true;
  }
  return false;
}

ColorID METS_Scale_Setter::CombineColors
(const Cluster_Leg *li,const Cluster_Leg *lj,const Cluster_Leg *lk,
 const ATOOLS::Flavour &mo) const
{
  ColorID ci(li->Col()), cj(lj->Col()), ck(lk->Col());
  if (ci.m_i==-1 && cj.m_i==-1 && ck.m_i==-1) return ColorID();
  if (!mo.Strong()) return ColorID(0,0);
  if (li->Flav().StrongCharge()==3) {
    if (lj->Flav().StrongCharge()==-3) {
      return ColorID(ci.m_i,cj.m_j);
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (cj.Singlet()) return ColorID(ci.m_i,0);
      return ColorID(cj.m_i,0);
    }
    else {
      return ColorID(ci.m_i,0);
    }
  }
  else if (li->Flav().StrongCharge()==-3) {
    if (lj->Flav().StrongCharge()==3) {
      return ColorID(cj.m_i,ci.m_j);
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (cj.Singlet()) return ColorID(0,ci.m_j);
      return ColorID(0,cj.m_j);
    }
    else {
      return ColorID(0,ci.m_j);
    }
  }
  else if (li->Flav().StrongCharge()==8) {
    if (lj->Flav().StrongCharge()==8) {
      if (ci.m_i==cj.m_j && 
	  (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	   (ci.m_j==cj.m_i && lk->Flav().StrongCharge()!=8))) 
	return ColorID(cj.m_i,ci.m_j);
      if (ci.m_j==cj.m_i && 
	  (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	   (ci.m_i==cj.m_j && lk->Flav().StrongCharge()!=8)))
	return ColorID(ci.m_i,cj.m_j);
      THROW(fatal_error,"Invalid clustering");
    }
    else if (lj->Flav().StrongCharge()==3) {
      if (ci.Singlet()) return ColorID(cj.m_i,0);
      return ColorID(ci.m_i,0);
    }
    else if (lj->Flav().StrongCharge()==-3) {
      if (ci.Singlet()) return ColorID(0,cj.m_j);
      return ColorID(0,ci.m_j);
    }
    else {
      THROW(fatal_error,"Invalid combination");
    }
  }
  else {
    if (lj->Flav().StrongCharge()==3) {
      return ColorID(cj.m_i,0);
    }
    else if (lj->Flav().StrongCharge()==-3) {
      return ColorID(0,cj.m_j);
    }
    else {
      return ColorID(0,0);
    }
  }
  return ColorID();
}
