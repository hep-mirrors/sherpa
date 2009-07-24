#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"

namespace PHASIC {

  class QCD_Scale_Setter: public Scale_Setter_Base {
  private:

    std::string m_muf2tag, m_mur2tag;

    ATOOLS::Algebra_Interpreter m_muf2calc, m_mur2calc;

    Tag_Setter m_muf2tagset, m_mur2tagset;

    ATOOLS::Vec4D_Vector   m_p;
    ATOOLS::Flavour_Vector m_f;

    SP(Color_Integrator) p_ci;

  public:

    QCD_Scale_Setter(Process_Base *const proc,
		     const std::string &scale);

    double CalculateScale(const std::vector<ATOOLS::Vec4D> &p);

    ATOOLS::Vec4D Momentum(const size_t &i) const;

    void SetScale(const std::string &mu2tag,Tag_Setter &mu2tagset,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class Scale_Setter_Base

  struct CS_Params {
    double m_kt2, m_z, m_y;
    int m_mode;
    CS_Params(const double &kt2,const double &z,
		  const double &y,const int mode=-1):
      m_kt2(kt2), m_z(z), m_y(y), m_mode(mode) {}
  };// end of struct CS_Params

  class QCD_Setter_CS_CD {
  private:

    Single_Process *p_proc;

    bool CheckColors(const ATOOLS::Cluster_Leg *li,
		     const ATOOLS::Cluster_Leg *lj,
		     const ATOOLS::Cluster_Leg *lk);
    ATOOLS::ColorID CombineColors(const ATOOLS::Cluster_Leg *li,
				  const ATOOLS::Cluster_Leg *lj,
				  const ATOOLS::Cluster_Leg *lk);

  public:

    inline QCD_Setter_CS_CD(Single_Process *const proc):
      p_proc(proc) {}

    CS_Params KT2(const ATOOLS::Cluster_Amplitude &ampl,
		  int i,int j,int k);

    void Combine(ATOOLS::Cluster_Amplitude &ampl,
		 int i,int j,int k);
    
  };// end of class QCD_Setter_CS_CD

  struct CKey {
    size_t m_i, m_j, m_k;
    CKey(const size_t &i,const size_t &j,const size_t &k):
      m_i(i),m_j(j), m_k(k) {}
    bool operator<(const CKey &ck) const
    { 
      if (m_i<ck.m_i) return true;
      if (m_i>ck.m_i) return false;
      if (m_j<ck.m_j) return true;
      if (m_j>ck.m_j) return false;
      return m_k<ck.m_k;
    }
  };// end of struct CKey

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(QCD_Scale_Setter_Getter,"QCD",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *QCD_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new QCD_Scale_Setter(args.p_proc,args.m_scale);
}

void QCD_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"qcd scale scheme\n";
}

QCD_Scale_Setter::QCD_Scale_Setter
(Process_Base *const proc,const std::string &scale): 
  Scale_Setter_Base(proc), m_muf2tagset(this), m_mur2tagset(this)
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
  for (size_t i(0);i<p_proc->NIn();++i) m_f[i]=m_f[i].Bar();
}

Vec4D QCD_Scale_Setter::Momentum(const size_t &i) const
{
  if (i>m_p.size()) THROW(fatal_error,"Momentum index too large");
  return m_p[i];
}

double QCD_Scale_Setter::CalculateScale(const std::vector<ATOOLS::Vec4D> &momenta) 
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
  if (p_proc->IsMapped() && p_proc->LookUp()) {
    m_kfkey[0]=m_scale[stp::ren]=
      p_proc->MapProc()->ScaleSetter()->Scale(stp::ren);
    m_kfkey[2]=m_kfkey[1]=m_scale[stp::fac]=
      p_proc->MapProc()->ScaleSetter()->Scale(stp::fac);
    return m_scale[stp::fac];
  }
  m_p=momenta;
  for (size_t i(0);i<p_proc->NIn();++i) m_p[i]=-m_p[i];
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
  double kt2max(0.0);
  QCD_Setter_CS_CD cd(p_proc->Get<Single_Process>());
  std::set<CKey> trials;
  while (ampl->Legs().size()>4) {
    double kt2w(std::numeric_limits<double>::max());
    size_t iw(0), jw(0), kw(0);
    CKey ckw(0,0,0);
    for (size_t i(0);i<ampl->Legs().size();++i)
      for (size_t j(Max((size_t)2,i+1));j<ampl->Legs().size();++j)
	for (size_t k(0);k<ampl->Legs().size();++k) 
	  if (k!=i && k!=j) {
	    CKey ck(ampl->Leg(i)->Id(),ampl->Leg(j)->Id(),
		    ampl->Leg(k)->Id());
	    if (trials.find(ck)!=trials.end()) {
	      msg_Debugging()<<"Vetoed "<<ID(ck.m_i)
			     <<" & "<<ID(ck.m_j)
			     <<" <-> "<<ID(ck.m_k)<<"\n";
	      continue;
	    }
	    CS_Params cs(cd.KT2(*ampl,i,j,k));
	    if (cs.m_kt2<kt2w) {
	      kt2w=cs.m_kt2;
	      ckw=ck;
	      iw=i;
	      jw=j;
	      kw=k;
	    }
	  }
    trials.insert(ckw);
    if (iw==0 && jw==0 && kw==0) {
      if (ampl->Prev()==NULL) {
	msg_Error()<<*ampl<<std::endl;
	THROW(fatal_error,"No valid clustering");
      }
      ampl=ampl->Prev();
      ampl->DeleteNext();
      continue;
    }
    kt2max=Max(kt2max,kt2w);
    msg_Debugging()<<"Actual = "<<*ampl<<"\n";
    msg_Debugging()<<"Cluster "<<ID(ampl->Leg(iw)->Id())
		   <<" & "<<ID(ampl->Leg(jw)->Id())
		   <<" <-> "<<ID(ampl->Leg(kw)->Id())
		   <<" => "<<sqrt(kt2w)<<"\n";
    ampl=ampl->InitNext();
    ampl->CopyFrom(ampl->Prev());
    cd.Combine(*ampl,iw,jw,kw);
  }
  msg_Debugging()<<"Core = "<<*ampl<<"\n";
  m_p.resize(ampl->Legs().size());
  Vec4D psum;
  int csum[4]={0,0,0,0};
  size_t qcd(0);
  ColorID c[4]={ampl->Leg(0)->Col(),ampl->Leg(1)->Col(),
		ampl->Leg(2)->Col(),ampl->Leg(3)->Col()};
  for (size_t i(0);i<m_p.size();++i) {
    Cluster_Leg *li(ampl->Leg(i));
    psum+=m_p[i]=li->Mom();
    ++csum[c[i].m_i];
    --csum[c[i].m_j];
    if (c[i].m_i>0 || c[i].m_j>0) qcd+=1<<i;
  }
  /*
  if (!IsEqual(psum,Vec4D(),1.0e-6))
    msg_Error()<<METHOD<<"(): Momentum not conserved. "<<*ampl<<std::endl;
  if (csum[1]!=0 || csum[2]!=0 || csum[3]!=0)
    msg_Error()<<METHOD<<"(): Colour not conserved. "<<*ampl<<std::endl;
  */
  double kt2cmin(std::numeric_limits<double>::max());
  if (qcd!=15) {
    if (p_ci==NULL) {
      bool s[4]={ampl->Leg(0)->Flav().Strong(),
		 ampl->Leg(1)->Flav().Strong(),
		 ampl->Leg(2)->Flav().Strong(),
		 ampl->Leg(3)->Flav().Strong()};
      if ((s[0] && s[1]) || (s[2] && s[3])) {
	kt2cmin=Min(kt2cmin,(m_p[0]+m_p[1]).Abs2());
      }
      if ((s[0] && s[2]) || (s[1] && s[3])) {
	kt2cmin=Min(kt2cmin,dabs((m_p[0]+m_p[2]).Abs2()));
      }
      if ((s[0] && s[3]) || (s[1] && s[2])) {
	kt2cmin=Min(kt2cmin,dabs((m_p[0]+m_p[3]).Abs2()));
      }
    }
    else {
      if ((c[0].m_i>0 && c[0].m_i==c[1].m_j) ||
	  (c[0].m_j>0 && c[0].m_j==c[1].m_i) ||
	  (c[2].m_i>0 && c[2].m_i==c[3].m_j) ||
	  (c[2].m_j>0 && c[2].m_j==c[3].m_i)) {
	kt2cmin=Min(kt2cmin,(m_p[0]+m_p[1]).Abs2());
      }
      if ((c[0].m_i>0 && c[0].m_i==c[2].m_j) ||
	  (c[0].m_j>0 && c[0].m_j==c[2].m_i) ||
	  (c[1].m_i>0 && c[1].m_i==c[3].m_j) ||
	  (c[1].m_j>0 && c[1].m_j==c[3].m_i)) {
	kt2cmin=Min(kt2cmin,dabs((m_p[0]+m_p[2]).Abs2()));
      }
      if ((c[0].m_i>0 && c[0].m_i==c[3].m_j) ||
	  (c[0].m_j>0 && c[0].m_j==c[3].m_i) ||
	  (c[1].m_i>0 && c[1].m_i==c[2].m_j) ||
	  (c[1].m_j>0 && c[1].m_j==c[2].m_i)) {
	kt2cmin=Min(kt2cmin,dabs((m_p[0]+m_p[3]).Abs2()));
      }
    }
  }
  if (kt2cmin==std::numeric_limits<double>::max()) {
    if (ampl->Leg(2)->Flav().IsMassive()) {
      if (ampl->Leg(3)->Flav().IsMassive()) {
	kt2cmin=sqrt(m_p[2].MPerp2()*m_p[3].MPerp2());
      }
      else {
	kt2cmin=m_p[2].MPerp2();
      }
    }
    else {
      if (ampl->Leg(3)->Flav().IsMassive()) {
	kt2cmin=m_p[3].MPerp2();
      }
      else {
	kt2cmin=m_p[3].PPerp2();
      }
    }
  }
  while (ampl->Prev()) ampl=ampl->Prev();
  ampl->Delete();
  m_scale[stp::ren]=m_scale[stp::fac]=Max(kt2max,kt2cmin);
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

void QCD_Scale_Setter::SetScale
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

CS_Params QCD_Setter_CS_CD::KT2
(const ATOOLS::Cluster_Amplitude &ampl,int i,int j,int k)
{
  static const CS_Params nd(std::numeric_limits<double>::max(),0.0,0.0,-1);
  const Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j)), *lk(ampl.Leg(k));
  if (!p_proc->Combinable(li->Id(),lj->Id())) return nd;
  if (!CheckColors(li,lj,lk)) return nd;
  if ((li->Id()&3)<(lj->Id()&3)) std::swap<const Cluster_Leg*>(li,lj);
  if ((li->Id()&3)==0) {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	Vec4D  Q(li->Mom()+lj->Mom()+lk->Mom());
	double pipj=li->Mom()*lj->Mom(), pipk=li->Mom()*lk->Mom();
	double pjpk=lj->Mom()*lk->Mom(), Q2=Q*Q, yijk=pipj/(pipj+pipk+pjpk);
	double zi=pipk/(pipk+pjpk), kt2=Q2*yijk*zi*(1.-zi);
	return CS_Params(kt2,zi,yijk,0);
      }
      else {
	double pipj=li->Mom()*lj->Mom(), pipa=li->Mom()*lk->Mom();
	double pjpa=lj->Mom()*lk->Mom(), xija=(pipa+pjpa+pipj)/(pipa+pjpa);
	double zi=pipa/(pipa+pjpa), kt2=-2.*(pipa+pjpa)*(1.-xija)*zi*(1.-zi);
	return CS_Params(kt2,zi,1.0-xija,2);
      }
    }
  }
  else {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	double pjpa=lj->Mom()*li->Mom(), pkpa=lk->Mom()*li->Mom();
	double pjpk=lj->Mom()*lk->Mom(), xjka=(pjpa+pkpa+pjpk)/(pjpa+pkpa);
	double uj=pjpa/(pjpa+pkpa), kt2=-2.*(pjpa+pkpa)*(1.-xjka)*uj; 
	return CS_Params(kt2,xjka,uj,1);
      }
      else {
	double papb=li->Mom()*lk->Mom(), pjpa=lj->Mom()*li->Mom();
	double pjpb=lj->Mom()*lk->Mom(), xjab=(papb+pjpa+pjpb)/papb;
	double vj=-pjpa/papb, kt2=2.*papb*vj*(1.-xjab);
	return CS_Params(kt2,xjab,vj,3);
      }
    }
  }
  THROW(fatal_error,"Unknown CS dipole configuration");  
}

void QCD_Setter_CS_CD::Combine
(Cluster_Amplitude &ampl,int i,int j,int k)
{
  if (i>j) std::swap<int>(i,j);
  Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j)), *lk(ampl.Leg(k));
  if (i>1 && j>1 && k>1) {
    Vec4D pi(li->Mom()), pj(lj->Mom()), pk(lk->Mom()), Q(pi+pj+pk);
    double Q2=Q*Q;
    Vec4D pkt(Q2/(Q2-(pi+pj).Abs2())*(pk-(Q*pk/Q2)*Q)+0.5*Q), pijt(Q-pkt);
    li->SetMom(pijt);
    lk->SetMom(pkt);
  }
  if (i>1 && j>1 && k<2) {
    Vec4D pi(li->Mom()), pj(lj->Mom()), pa(lk->Mom());
    double pipa=pi*pa, pjpa=pj*pa, pipj=pi*pj;
    double xija=(pipa+pjpa+pipj)/(pipa+pjpa);
    Vec4D pijt(pi+pj+(1.-xija)*pa), pat(xija*pa);
    li->SetMom(pijt);
    lk->SetMom(pat);
  }
  if (i<2 && j>1 && k>1) {
    Vec4D pa(li->Mom()), pj(lj->Mom()), pk(lk->Mom());
    double pjpa=pj*pa, pkpa=pk*pa, pjpk=pj*pk;
    double xjka=(pjpa+pkpa+pjpk)/(pjpa+pkpa);
    Vec4D pajt(xjka*pa), pkt(pj+pk+(1.-xjka)*pa);
    li->SetMom(pajt);
    lk->SetMom(pkt);
  }
  if (i<2 && j>1 && k<2) {
    Vec4D pa(li->Mom()), pj(lj->Mom()), pb(lk->Mom());
    double papb=pa*pb, pjpa=pj*pa, pjpb=pj*pb;
    double xjab=(papb+pjpa+pjpb)/papb;
    Vec4D pajt(xjab*pa), K(-pa-pb-pj), Kt(-pajt-pb), KpKt(K+Kt);
    for (size_t m(0);m<ampl.Legs().size();++m) {
      if (m==(size_t)j) continue;
      Vec4D km = ampl.Leg(m)->Mom();
      km=km-2.*km*KpKt/(KpKt*KpKt)*KpKt+2.*km*K/(K*K)*Kt;
      ampl.Leg(m)->SetMom(km);
    }
    li->SetMom(pajt);
    lk->SetMom(pb);
  }
  li->SetCol(CombineColors(li,lj,lk));
  li->SetId(li->Id()+lj->Id());
  const Flavour_Vector &cf(p_proc->CombinedFlavour(li->Id()));
  li->SetFlav(cf.front());
  for (size_t i(0);i<cf.size();++i)
    if (cf[i].Strong()) {
      li->SetFlav(cf[i]);
      break;
    }
  std::vector<Cluster_Leg*>::iterator lit(ampl.Legs().begin());
  for (int l(0);l<j;++l) ++lit;
  (*lit)->Delete();
  ampl.Legs().erase(lit);
}

bool QCD_Setter_CS_CD::CheckColors
(const Cluster_Leg *li,const Cluster_Leg *lj,const Cluster_Leg *lk)
{
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

ColorID QCD_Setter_CS_CD::CombineColors
(const Cluster_Leg *li,const Cluster_Leg *lj,const Cluster_Leg *lk)
{
  ColorID ci(li->Col()), cj(lj->Col()), ck(lk->Col());
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
