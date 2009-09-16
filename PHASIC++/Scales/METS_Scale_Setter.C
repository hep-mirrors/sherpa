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
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"

#define CHECK__x

namespace PHASIC {

  struct CS_Params {
    double m_kt2, m_op2, m_z, m_y;
    ATOOLS::Flavour m_fl;
    int m_mode;
    CS_Params(const double &kt2,const double &z,
	      const double &y,const ATOOLS::Flavour &fl,
	      const int mode=-1):
      m_kt2(kt2), m_op2(0.0), m_z(z), m_y(y), m_fl(fl), m_mode(mode) {}
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

    CS_Params KT2(const ATOOLS::Cluster_Leg *li,
		  const ATOOLS::Cluster_Leg *lj,
		  const ATOOLS::Cluster_Leg *lk,
		  const ATOOLS::Flavour &mo) const;

    bool Combine(ATOOLS::Cluster_Amplitude &ampl,
		 int i,int j,int k,const ATOOLS::Flavour &mo) const;

    double SetScales(const double &scale);

    double CalculateStrict(const ATOOLS::Vec4D_Vector &momenta);

  public:

    METS_Scale_Setter(Process_Base *const proc,
		      const std::string &scale,const int mode=1);

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);
    double CalculateScale2(const std::vector<ATOOLS::Vec4D> &p);

    ATOOLS::Vec4D Momentum(const size_t &i) const;

    void SetScale(const std::string &mu2tag,Tag_Setter &mu2tagset,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class METS_Scale_Setter

  struct MCKey {
    size_t m_i, m_j, m_k;
    ATOOLS::Flavour m_fl;
    MCKey(const size_t &i,const size_t &j,const size_t &k,
	  const ATOOLS::Flavour &fl):
      m_i(i),m_j(j), m_k(k), m_fl(fl) {}
    bool operator<(const MCKey &ck) const
    { 
      if (m_i<ck.m_i) return true;
      if (m_i>ck.m_i) return false;
      if (m_j<ck.m_j) return true;
      if (m_j>ck.m_j) return false;
      if (m_k<ck.m_k) return true;
      if (m_k>ck.m_k) return false;
      return m_fl<ck.m_fl;
    }
  };// end of struct MCKey

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Loose_METS_Scale_Setter_Getter,"LOOSE_METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *Loose_METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args.p_proc,args.m_scale,0);
}

void Loose_METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"loose mets scale scheme\n";
}

DECLARE_GETTER(METS_Scale_Setter_Getter,"METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args.p_proc,args.m_scale);
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
  return new METS_Scale_Setter(args.p_proc,args.m_scale,2);
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
(Process_Base *const proc,const std::string &scale,const int mode):
  Scale_Setter_Base(proc), m_muf2tagset(this), m_mur2tagset(this),
  m_cnt(0), m_rej(0), m_mode(mode), m_lfrac(0.0)
{
  m_p.resize(4);
  size_t pos(scale.find('{'));
  std::string mur2tag("MU_R2"), muf2tag("MU_F2");
  if (pos!=std::string::npos) {
    muf2tag=scale.substr(pos+1);
    pos=muf2tag.rfind('}');
    if (pos==std::string::npos)
      THROW(fatal_error,"Invalid scale '"+scale+"'");
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
  p_proc->Integrator()->SetMomenta(momenta);
  p_proc->Generator()->SetClusterDefinitions
    (p_proc->Shower()->GetClusterDefinitions());
  Cluster_Amplitude *ampl
    (p_proc->Generator()->ClusterConfiguration(p_proc));
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
  while (ampl->Next()) ampl=ampl->Next();
  double kt2max(ampl->KT2QCD());
  msg_Debugging()<<"Core = "<<*ampl<<"\n";
  m_p.resize(ampl->Legs().size());
  for (size_t i(0);i<m_p.size();++i)
    m_p[i]=ampl->Leg(i)->Mom();
  return SetScales(kt2max);
}

double METS_Scale_Setter::CalculateScale2(const Vec4D_Vector &momenta) 
{
  if (m_mode==2 || (m_mode==1 && !p_proc->LookUp())) {
    p_proc->Integrator()->SetMomenta(momenta);
    p_proc->Integrator()->SwapInOrder();
    double muf2(CalculateScale(p_proc->Integrator()->Momenta()));
    p_proc->Integrator()->RestoreInOrder();
    return muf2;
  }
  return m_scale[stp::fac];
}

double METS_Scale_Setter::CalculateScale(const Vec4D_Vector &momenta) 
{
  if (!m_kfkey.Assigned()) {
    std::string kfinfo("O(QCD)="+ToString(p_proc->OrderQCD()));
    msg_Debugging()<<"Assign '"<<p_proc->Name()
		   <<"' '"<<kfinfo<<"'\n";
    m_kfkey.Assign(p_proc->Name(),3,0,p_proc->
		   Integrator()->PSHandler()->GetInfo());
    m_kfkey.SetInfo(kfinfo);
    p_ci=p_proc->Integrator()->ColorIntegrator();
  }
  ++m_cnt;
  m_p=momenta;
  for (size_t i(0);i<p_proc->NIn();++i) m_p[i]=-m_p[i];
  if (m_mode==2 || (m_mode==1 && !p_proc->LookUp()))
    return CalculateStrict(momenta);
  if (p_proc->IsMapped() && p_proc->LookUp()) {
    m_kfkey[0]=m_scale[stp::ren]=
      p_proc->MapProc()->ScaleSetter()->Scale(stp::ren);
    m_kfkey[2]=m_kfkey[1]=m_scale[stp::fac]=
      p_proc->MapProc()->ScaleSetter()->Scale(stp::fac);
    return m_scale[stp::fac];
  }
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
  ampl->SetX1(p_proc->Integrator()->ISR()->X1());
  ampl->SetX2(p_proc->Integrator()->ISR()->X2());
  Single_Process *proc(p_proc->Get<Single_Process>());
  std::set<MCKey> trials;
  std::vector<double> ops(ampl->Legs().size()-3,0.0);
  while (ampl->Legs().size()>4) {
    double op2w(0.0), kt2w(0.0);
    size_t iw(0), jw(0), kw(0);
    MCKey ckw(0,0,0,kf_none);
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
	      MCKey ck(li->Id(),lj->Id(),lk->Id(),cf[f]);
	      if (trials.find(ck)!=trials.end()) {
		msg_Debugging()<<"Vetoed "<<cf[f]<<" = "
			       <<ID(ck.m_i)<<" & "<<ID(ck.m_j)
			       <<" <-> "<<ID(ck.m_k)<<"\n";
		continue;
	      }
	      if (!CheckColors(li,lj,lk,cf[f])) {
		trials.insert(ck);
		continue;
	      }
	      CS_Params cs(KT2(li,lj,lk,cf[f]));
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
		    cs.m_op2*=dabs(s)/
		      sqrt(sqr(s-m2)+m2*sqr(cf[f].Width()));
		  }
		  else {
		    // if non-resonant, reweight with massive prop
		    cs.m_op2*=cs.m_kt2/(cs.m_kt2+sqr(cf[f].Mass()));
		  }
		}
	      }
	      if (cs.m_op2>op2w) {
		op2w=cs.m_op2;
		kt2w=cs.m_kt2;
		ckw=ck;
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
      continue;
    }
    msg_Debugging()<<"Actual = "<<*ampl<<"\n";
    msg_Debugging()<<"Cluster "<<ckw.m_fl<<" "
		   <<ID(ckw.m_i)<<" & "<<ID(ckw.m_j)
		   <<" <-> "<<ID(ckw.m_k)
		   <<" => "<<sqrt(kt2w)
		   <<" ("<<sqrt(op2w)<<") <-> "
		   <<sqrt(ops[ampl->Legs().size()-4])<<"\n";
    if (ops[ampl->Legs().size()-4]>kt2w) {
      msg_Debugging()<<"unordered configuration\n";
#ifdef METS__reject_unordered
      continue;
#endif
    }
    ampl->SetKT2QCD(kt2w);
    ampl=ampl->InitNext();
    ampl->CopyFrom(ampl->Prev());
    if (!Combine(*ampl,iw,jw,kw,ckw.m_fl)) {
      ampl=ampl->Prev();
      ampl->DeleteNext();
      continue;
    }
    ops[ampl->Legs().size()-4]=kt2w;
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
      size_t qcd(0), qc(0);
      for (size_t i(0);i<4;++i) {
	int cc(ampl->Leg(i)->Flav().StrongCharge());
	qc+=abs(cc)==3;
	qcd+=cc;
      }
      if (qcd%8!=0 || (qc==0 && qcd<32)) {
	ampl=ampl->Prev();
	ampl->DeleteNext();
	continue;
      }
    }
  }
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
  while (ampl->Prev()) ampl=ampl->Prev();
  ampl->Delete();
  return SetScales(kt2cmin);
}

double METS_Scale_Setter::SetScales(const double &scale)
{
  m_scale[stp::ren]=m_scale[stp::fac]=scale;
  msg_Debugging()<<"QCD scale = "<<sqrt(m_scale[stp::ren])<<"\n";
  m_scale[stp::ren]=m_mur2calc.Calculate()->Get<double>();
  m_scale[stp::fac]=m_muf2calc.Calculate()->Get<double>();
  msg_Debugging()<<"Set \\mu_r = "
		 <<sqrt(m_scale[stp::ren])<<", \\mu_f = "
		 <<sqrt(m_scale[stp::fac])<<"\n";
  m_kfkey[0]=m_scale[stp::ren];
  m_kfkey[2]=m_kfkey[1]=m_scale[stp::fac];
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
  mu2calc.AddTag("MU_F2","1.0");
  mu2calc.AddTag("MU_R2","1.0");
  mu2calc.AddTag("H_T2","1.0");
  mu2calc.AddTag("Q2_CUT","1.0");
  mu2calc.AddTag("Q2_MIN","1.0");
  Process_Integrator *ib(p_proc->Integrator());
  for (size_t i=0;i<ib->NIn()+ib->NOut();++i) 
    mu2calc.AddTag("p["+ToString(i)+"]",ToString(ib->Momenta()[i]));
  mu2calc.Interprete(mu2tag);
  msg_Debugging()<<"}\n";
}

double METS_Scale_Setter::Lam
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

CS_Params METS_Scale_Setter::KT2
(const Cluster_Leg *li,const Cluster_Leg *lj,
 const Cluster_Leg *lk,const Flavour &mo) const
{
  if ((li->Id()&3)<(lj->Id()&3)) std::swap<const Cluster_Leg*>(li,lj);
  if ((li->Id()&3)==0) {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	Vec4D pk(lk->Mom()), Q(li->Mom()+lj->Mom()+pk);
	double pipj=li->Mom()*lj->Mom(), pipk=li->Mom()*lk->Mom();
	double pjpk=lj->Mom()*lk->Mom(), Q2=Q*Q, yijk=pipj/(pipj+pipk+pjpk);
	double mi2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
	double mij2=sqr(mo.Mass()), mk2=sqr(lk->Flav().Mass());
	double zi=pipk/(pipk+pjpk);
	double kt2=(Q2-mi2-mj2-mk2)*yijk*zi*(1.-zi)
	  -(1.0-zi)*(1.0-zi)*mi2-zi*zi*mj2;
	double lrat=Lam(Q2,mij2,mk2)/Lam(Q2,(li->Mom()+lj->Mom()).Abs2(),mk2);
	Vec4D pkt(sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-mij2)/(2.0*Q2)*Q);
	if (lrat<0.0 || IsZero(lrat,1.0e-6) ||
	    pkt[0]<0.0 || IsZero(pkt[0],1.0e-6))
	  kt2=s_kt2max;
	return CS_Params(kt2,zi,yijk,mo,0);
      }
      else {
	double pipj=li->Mom()*lj->Mom(), pipa=li->Mom()*lk->Mom();
	double pjpa=lj->Mom()*lk->Mom(), xija=(pipa+pjpa+pipj)/(pipa+pjpa);
	double mi2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
	double ma2=sqr(lk->Flav().Mass()), zi=pipa/(pipa+pjpa);
	double arg=sqr(pipa+pjpa)-ma2*(li->Mom()+lj->Mom()).Abs2();
	double kt2=-2.0*(pipa+pjpa)*(1.0-xija)*zi*(1.0-zi)
	  -sqr(1.0-zi)*mi2-zi*zi*mj2;
	if (kt2<0.0 || IsZero(kt2,1.0e-6) ||
	    xija<0.0 || IsZero(xija,1.0e-6) || arg<0.0)
	  kt2=s_kt2max;
	return CS_Params(kt2,zi,1.0-xija,mo,2);
      }
    }
  }
  else {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	double pjpa=lj->Mom()*li->Mom(), pkpa=lk->Mom()*li->Mom();
	double pjpk=lj->Mom()*lk->Mom(), xjka=(pjpa+pkpa+pjpk)/(pjpa+pkpa);
	double ma2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
	double uj=pjpa/(pjpa+pkpa);
	double arg=sqr(pjpa+pkpa)-ma2*(lj->Mom()+lk->Mom()).Abs2();
	double kt2=-2.*(pjpa+pkpa)*(1.-xjka)*uj-mj2-sqr(1.0-xjka)*ma2;
	if (kt2<0.0 || IsZero(kt2,1.0e-6) ||
	    xjka<0.0 || IsZero(xjka,1.0e-6) || arg<0.0)
	  kt2=s_kt2max;
	return CS_Params(kt2,xjka,uj,mo,1);
      }
      else {
	double papb=li->Mom()*lk->Mom(), pjpa=lj->Mom()*li->Mom();
	double pjpb=lj->Mom()*lk->Mom(), xjab=(papb+pjpa+pjpb)/papb;
	double mj2=sqr(lj->Flav().Mass()), ma2=sqr(li->Flav().Mass());
	double maj2=sqr(mo.Mass()), mb2=sqr(lk->Flav().Mass());
	double Q2=(li->Mom()+lj->Mom()+lk->Mom()).Abs2(), vj=-pjpa/papb;
	double kt2=2.0*papb*vj*(1.0-xjab)-mj2-sqr(1.0-xjab)*ma2;
	double ttau=Q2-maj2-mb2, tau=Q2-ma2-mj2-mb2;
	if (xjab<0.0 || IsZero(xjab,1.0e-6) ||
	    ttau<0.0 || IsZero(ttau,1.0e-6) ||
	    tau<0.0 || IsZero(tau,1.0e-6))
	  kt2=s_kt2max;
	return CS_Params(kt2,xjab,vj,mo,3);
      }
    }
  }
  THROW(fatal_error,"Unknown CS dipole configuration");  
}

bool METS_Scale_Setter::Combine
  (Cluster_Amplitude &ampl,int i,int j,int k,const Flavour &mo) const
{
  if (i>j) std::swap<int>(i,j);
  Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j)), *lk(ampl.Leg(k));
  if (i>1 && j>1 && k>1) {
    Vec4D pi(li->Mom()), pj(lj->Mom()), pk(lk->Mom()), Q(pi+pj+pk);
    double Q2=Q*Q, mij2=sqr(mo.Mass()), mk2=sqr(lk->Flav().Mass());
    double lrat=Lam(Q2,mij2,mk2)/Lam(Q2,(pi+pj)*(pi+pj),mk2);
    Vec4D pkt(sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-mij2)/(2.*Q2)*Q);
    Vec4D pijt(Q-pkt); 
    if (pkt[0]<0.0 || pijt[0]<0.0) return false;
    li->SetMom(pijt);
    lk->SetMom(pkt);
  }
  if (i>1 && j>1 && k<2) {
    Vec4D pi(li->Mom()), pj(lj->Mom()), pa(lk->Mom()), Q(pa+pi+pj);
    double mij2=sqr(mo.Mass()), ma2=sqr(lk->Flav().Mass());
    double mb2=sqr(ampl.Leg(1-k)->Flav().Mass()), Q2=Q.Abs2();
    double lrat=Lam(Q2,mij2,ma2)/Lam(Q2,(pi+pj).Abs2(),ma2);
    if (lrat<0.0) return false;
    Vec4D pat(sqrt(lrat)*(pa-(Q*pa/Q2)*Q)+(Q2+ma2-mij2)/(2.*Q2)*Q);
    Vec4D pijt(Q-pat), pb(ampl.Leg(1-k)->Mom()), path(pat);
    if (pijt[0]<0.0) return false;
    double patpb=pat*pb, sb=Sign(pb[3]), ea=0.0, s=(pat+pb).Abs2();
    if (IsZero(mb2)) ea=0.5*(patpb+ma2*sqr(pb[3])/patpb)/pb[0];
    else ea=(pb[0]*patpb+dabs(pb[3])*sqrt(patpb*patpb-ma2*mb2))/mb2;
    Vec4D pan(ea,0.0,0.0,-sb*sqrt(ea*ea-ma2)), pam(ea,0.0,0.0,-pan[3]);
    if (dabs((pam+pb).Abs2()-s)<dabs((pan+pb).Abs2()-s)) pan=pam;
    if (ea>0.0 || IsZero(ea,1.0e-6) ||
	patpb*patpb<ma2*mb2) return false;
#ifdef CHECK__x
    if (!CheckX(pan,lk->Id()&3)) return false;
#endif
    Poincare cmso(-pat-pb), cmsn(-pan-pb);
    cmso.Boost(pat);
    Poincare zrot(pat,-sb*Vec4D::ZVEC);
    for (size_t m(0);m<ampl.Legs().size();++m) {
      if (m==(size_t)j) continue;
      {
	Vec4D cm(ampl.Leg(m)->Mom());
	if (m==(size_t)i) cm=pijt;
	else if (m==(size_t)k) cm=path;
	cmso.Boost(cm);
	zrot.Rotate(cm);
	cmsn.BoostBack(cm);
	ampl.Leg(m)->SetMom(cm);
      }
    }
    double x(1.0+(pi*pj)/((pi+pj)*pa));
    if (k==0) ampl.SetX1(ampl.X1()*x);
    else ampl.SetX2(ampl.X2()*x);
  }
  if (i<2 && j>1 && k>1) {
    Vec4D pa(li->Mom()), pj(lj->Mom()), pk(lk->Mom()), Q(pa+pj+pk);
    double pjpa=pj*pa, pkpa=pk*pa, pjpk=pj*pk;
    double xjka=(pjpa+pkpa+pjpk)/(pjpa+pkpa);
    double mj2=sqr(lj->Flav().Mass()), mk2=sqr(lk->Flav().Mass());
    double ma2=sqr(li->Flav().Mass()), maj2=sqr(mo.Mass());
    double mb2=sqr(ampl.Leg(1-i)->Flav().Mass());
    double Q2=Q.Abs2(), ttau=Q2-maj2-mk2, tau=Q2-ma2-mj2-mk2;
    double sjk=-((1.0-xjka)*(Q2-ma2)-(mj2+mk2))/xjka;
    if (ttau*ttau<4.*maj2*mk2 || ttau>0.0 ||
	tau*tau<4.*ma2*sjk*sqr(xjka) || tau>0.0) return false;
    double xijka=xjka*(ttau-sqrt(ttau*ttau-4.*maj2*mk2))/
      (tau-sqrt(tau*tau-4.*ma2*sjk*sqr(xjka)));
    double pjkpa=pjpa+pkpa, gam=-pjkpa+sqrt(pjkpa*pjkpa-ma2*sjk);
    if (IsZero(xijka,1.0e-6) || IsZero(gam,1.0e-6)) return false;
    double bet=1.0-ma2*sjk/(gam*gam), gamt=gam*xijka;
    Vec4D l((-pa-ma2/gam*(pj+pk))/bet), n(((pj+pk)+sjk/gam*pa)/bet);
    l*=(1.0-sjk/gam)/(1.0-mk2/gamt);
    n*=(1.0-ma2/gam)/(1.0-maj2/gamt);
    Vec4D pat(-l-maj2/gamt*n), pjkt(n+mk2/gamt*l), pb(ampl.Leg(1-i)->Mom());
    if (pat[3]*pb[3]>0.0 || pjkt[0]<0.0) return false;
    double patpb=pat*pb, sb=Sign(pb[3]), ea=0.0, s=(pat+pb).Abs2();
    if (IsZero(mb2)) ea=0.5*(patpb+maj2*sqr(pb[3])/patpb)/pb[0];
    else ea=(pb[0]*patpb+dabs(pb[3])*sqrt(patpb*patpb-maj2*mb2))/mb2;
    Vec4D pan(ea,0.0,0.0,-sb*sqrt(ea*ea-maj2)), pam(ea,0.0,0.0,-pan[3]);
    if (dabs((pam+pb).Abs2()-s)<dabs((pan+pb).Abs2()-s)) pan=pam;
    if (ea>0.0 || IsZero(ea,1.0e-6) ||
	patpb*patpb<maj2*mb2) return false;
#ifdef CHECK__x
    if (!CheckX(pan,li->Id()&3)) return false;
#endif
    Vec4D path(pat);
    Poincare cmso(-pat-pb), cmsn(-pan-pb);
    cmso.Boost(pat);
    Poincare zrot(pat,-sb*Vec4D::ZVEC);
    for (size_t m(0);m<ampl.Legs().size();++m) {
      if (m==(size_t)j) continue;
      {
	Vec4D cm(ampl.Leg(m)->Mom());
	if (m==(size_t)i) cm=path;
	else if (m==(size_t)k) cm=pjkt;
	cmso.Boost(cm);
	zrot.Rotate(cm);
	cmsn.BoostBack(cm);
	ampl.Leg(m)->SetMom(cm);
      }
    }
    double x(1.0+(pj*pk)/((pj+pk)*pa));
    if (i==0) ampl.SetX1(ampl.X1()*x);
    else ampl.SetX2(ampl.X2()*x);
  }
  if (i<2 && j>1 && k<2) {
    Vec4D pa(li->Mom()), pj(lj->Mom()), pb(lk->Mom());
    double papb=pa*pb, pjpa=pj*pa, pjpb=pj*pb;
    double mj2=sqr(lj->Flav().Mass()), ma2=sqr(li->Flav().Mass());
    double mb2=sqr(ampl.Leg(1-i)->Flav().Mass());
    double maj2=sqr(mo.Mass()), Q2=(pa+pj+pb).Abs2();
    double xjab=(papb+pjpa+pjpb)/papb;
    double ttau=Q2-maj2-mb2, tau=Q2-ma2-mj2-mb2;
    if (ttau*ttau<4.0*maj2*mb2 ||
	tau*tau<4.0*ma2*mb2*xjab*xjab) return false;
    double xijab=xjab*(ttau+sqrt(ttau*ttau-4.0*maj2*mb2))
      /(tau+sqrt(tau*tau-4.0*ma2*mb2*xjab*xjab));
    double gam=papb+sqrt(papb*papb-ma2*mb2);
    if (IsZero(xijab,1.0e-6) || IsZero(gam,1.0e-6)) return false;
    Vec4D pajt=xijab
      *(1.0-maj2*mb2/sqr(gam*xijab))/(1.0-ma2*mb2/sqr(gam))
      *(pa-ma2/gam*pb)+maj2/(xijab*gam)*pb;
    if (pajt[3]*pb[3]>0.0) return false;
#ifdef CHECK__x
    if (!CheckX(pajt,li->Id()&3)) return false;
#endif
    Vec4D K(-pa-pb-pj), Kt(-pajt-pb), KpKt(K+Kt);
    for (size_t m(0);m<ampl.Legs().size();++m) {
      if (m==(size_t)j) continue;
      if (m==(size_t)i) ampl.Leg(m)->SetMom(pajt);
      else if (m==(size_t)k) ampl.Leg(m)->SetMom(pb);
      else {
	Vec4D km = ampl.Leg(m)->Mom();
	km=km-2.*km*KpKt/(KpKt*KpKt)*KpKt+2.*km*K/(K*K)*Kt;
	ampl.Leg(m)->SetMom(km);
      }
    }
    double x(1.0+(pj*(pa+pb))/(pa*pb));
    if (i==0) ampl.SetX1(ampl.X1()*x);
    else ampl.SetX2(ampl.X2()*x);
  }
  li->SetCol(CombineColors(li,lj,lk,mo));
  li->SetId(li->Id()+lj->Id());
  li->SetFlav(mo);
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
  if (li->Col().m_i==-1 && lj->Col().m_i==-1 &&
      lk->Col().m_i==-1) return true;
  if (!mo.Strong()) {
    if (lk->Flav().StrongCharge()==8) return false;
    ColorID ci(li->Col()), cj(lj->Col());
    if (ci.m_i==cj.m_j && ci.m_j==0 && cj.m_i==0) return true;
    if (ci.m_j==cj.m_i && ci.m_i==0 && cj.m_j==0) return true;
    return false;
  }
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
