#include "Higgs_Virtual.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"

#include "Wrappers.H"
#include "Ahiggs.h"
#include "A_higgs.h"
#include "A_spin2.h"
#include "Acont.h"
#include "A_cont.h"

using namespace HIGGS;
using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;

MODEL::Model_Base *HIGGS::Higgs_Virtual::s_model=NULL;

Higgs_Virtual::Higgs_Virtual(const Process_Info &pi,
			     const Flavour_Vector &flavs,
			     int mode,int io,int spin,
			     double kg, double kq):
  Virtual_ME2_Base(pi,flavs), m_int(mode), m_io(io),
  m_spin(spin), m_kg(kg), m_kq(kq)
{
  m_mh=Flavour(kf_h0).Mass();
  m_gh=Flavour(kf_h0).Width();
  m_b=std::vector<int>(4,1);
  m_b[1]=m_b[0]=-1;
  p_bs = new Basic_Sfuncs
    (m_flavs.size(),m_flavs.size(),
     (Flavour*)&m_flavs.front(),&m_b.front());  
  p_bs->Initialize();
  m_proc=1;
  if (m_flavs[0].IsQuark() &&
      m_flavs[1]==m_flavs[0].Bar()) {
    if (m_flavs[0].IsAnti()) m_proc=5;
    else m_proc=4;
  }
}

Higgs_Virtual::~Higgs_Virtual()
{
  delete p_bs;
}

void Higgs_Virtual::Calc(const Vec4D_Vector &p)
{
  DEBUG_FUNC(this<<", m_proc = "<<m_proc);
  mu_sq=m_mur2;
  double muR=sqrt(m_mur2);
  alpha0=s_model->ScalarConstant("alpha_QED(0)");
  msg_Debugging()<<"\\mu_R = "<<muR<<"\n";
  for (size_t i(0);i<p.size();++i)
    msg_Debugging()<<"p["<<i<<"]="<<p[i]<<"\n";
  double m_kgr=kgr(muR,m_kg,m_kq,m_mh);
  double m_kqr=kqr(muR,m_kg,m_kq,m_mh);
  s_bs=p_bs;
  p_bs->Setk0(11);
  p_bs->CalcEtaMu((Vec4D*)&p.front());
  double s=(p[2]+p[3]).Abs2(), rts=sqrt(s);
  Complex lo=0.0, nlo=0.0;
  Complex los=0.0, nlos=0.0;
  Complex lob=0.0, nlob=0.0;
  Complex fslo=A_prod_1l(rts,muR)/s*A_dec_1l(rts,muR)/s/
    ((s-m_mh*m_mh)+I*m_mh*m_gh);
  if (m_spin!=0)
    fslo=A_prod_1l(rts,muR).real()/s*A_dec_1l(rts,muR).real()/s/
      ((s-m_mh*m_mh)+I*m_mh*m_gh);
  Complex fsnlo=(A_prod_2l(rts,muR)/s*A_dec_1l(rts,muR)/s+
		 A_prod_1l(rts,muR)/s*A_dec_2l(rts,muR)/s)/
    ((s-m_mh*m_mh)+I*m_mh*m_gh)/(alpha_s(muR)/(2.0*M_PI));
  Complex fblo=-4.0*sumQsq*alpha0*alpha_s(muR);
  double qc=m_flavs[m_flavs[0].IsQuark()?0:1].Charge();
  Complex ft=-2.0*(4.0*M_PI*alpha0)*sqr(qc);
  for (int i(1);i>=-1;i-=2) {
    for (int j(1);j>=-1;j-=2) {
      for (int k(1);k>=-1;k-=2) {
	for (int l(1);l>=-1;l-=2) {
	  Complex clos(0.0), cnlos(0.0);
	  Complex clob(0.0), cnlob(0.0);
	  Complex clo(0.0), cnlo(0.0);
	  if (m_proc==1) {
	    if (m_int&1) {
	      if (m_spin!=0) {
		Complex met=fslo*ggXgamgam(i,j,k,l,m_kgr);
		Complex me1l=fslo*ggXgamgam1l(i,j,k,l,muR,m_kgr,m_kqr);
		clos+=met; clo+=met;
		cnlos+=me1l; cnlo+=me1l;
	      }
	      else {
		Complex me1l=fslo*ggH(i,j)*Hgamgam(k,l);
		Complex me2l=fsnlo*ggH(i,j)*Hgamgam(k,l);
		clos+=me1l; clo+=me1l;
		cnlos+=me2l; cnlo+=me2l;
	      }
	    }
	    if (m_int&2) {
	      Complex me1l=gggamgam(i,j,k,l);
	      Complex me2l=fblo*gggamgam2l(i,j,k,l)*
		(me1l==0.0?1.0:me1l/gggamgam1l(i,j,k,l));
	      me1l*=fblo;
	      clob+=me1l; clo+=me1l;
	      cnlob+=me2l; cnlo+=me2l;
	    }
	  }
	  if (m_proc==4) {
	    if ((m_int&1) && i!=j) {
	      if (m_spin!=0) {
		Complex met=fslo*qqbXgamgam(i,k,l,m_kqr);
		Complex me1l=fslo*qqbXgamgam1l(i,k,l,muR,m_kgr,m_kqr);
		clos+=met; clo+=met;
		cnlos+=me1l; cnlo+=me1l;
	      }
	    }
	    if ((m_int&4) && i!=j) {
	      Complex met=ft*qqbgamgam_tree(i,k,l);
	      Complex me1l=ft*qqbyy1l(i,k,l).f/C_F;
	      clob+=met; clo+=met;
	      cnlob+=me1l; cnlo+=me1l;
	    }
	  }
	  if (m_proc==5) {
	    if ((m_int&1) && j!=i) {
	      if (m_spin!=0) {
		Complex met=fslo*qbqXgamgam(j,k,l,m_kqr);
		Complex me1l=fslo*qbqXgamgam1l(j,k,l,muR,m_kgr,m_kqr);
		clos+=met; clo+=met;
		cnlos+=me1l; cnlo+=me1l;
	      }
	    }
	    if ((m_int&4) && j!=i) {
	      Complex met=ft*qbqgamgam_tree(j,k,l);
	      Complex me1l=ft*qbqyy1l(j,k,l).f/C_F;
	      clob+=met; clo+=met;
	      cnlob+=me1l; cnlo+=me1l;
	    }
	  }
	  if (m_io==1) {
	    los+=clos*std::conj(clos);
	    nlos+=2.0*clos*std::conj(cnlos);
	    lob+=clob*std::conj(clob);
	    nlob+=2.0*clob*std::conj(cnlob);
	  }
	  if (m_io==2) {
	    lob+=clob*std::conj(clob);
	    nlob+=2.0*clob*std::conj(cnlob);
	  }
	  lo+=clo*std::conj(clo);
	  nlo+=2.0*clo*std::conj(cnlo);
	}
      }
    }
  }
  if (m_io==1) {
    lo-=los+lob;
    nlo-=nlos+nlob;
  }
  if (m_io==2) {
    lo-=lob;
    nlo-=nlob;
  }
  if (nlo.real()==0.0 && lo.real()==0.0) {
    m_res.IR2()=m_res.IR()=0.0;
    m_res.Finite()=m_born=0.0;
    return;
  }
  if (m_proc==4 || m_proc==5) {
    double lmur=log(m_mur2/s);
    m_res.IR2()=-2.0*4.0/3.0;
    m_res.IR()=-(3.0+2.0*lmur)*4.0/3.0;
    // extra single power of lmur term here compared to the gluon case
    m_res.Finite()=(nlo/lo).real()*4.0/3.0+(sqr(M_PI)-lmur*lmur-3.*lmur)*4.0/3.0;
    m_born=lo.real()/24.0;
    return;
  }
  double b0=(11.0*3.0-2.0*(Flavour(kf_quark).Size()/2))/6.0;
  double lmur=log(m_mur2/s);
  m_res.IR2()=-2.0*3.0;
  m_res.IR()=-2.0*(b0+3.0*lmur);
  // no extra single power of lmur term here compared to the quark case
  m_res.Finite()=(nlo/lo).real()+3.0*sqr(M_PI)-3.0*lmur*lmur;
  m_born=lo.real()/64.0;
}

double Higgs_Virtual::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

DECLARE_VIRTUALME2_GETTER(Higgs_Virtual,"Higgs_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,Higgs_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Higgs") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype==nlo_type::loop) {
    Data_Reader read(" ",";","#","=");
    int io=read.GetValue<int>("HIGGS_INTERFERENCE_ONLY",0);
    int mode=read.GetValue<int>("HIGGS_INTERFERENCE_MODE",7);
    int spin=read.GetValue<int>("HIGGS_INTERFERENCE_SPIN",0);
    double kg=read.GetValue<double>("HIGGS_INTERFERENCE_KAPPAG",1.0);
    double kq=read.GetValue<double>("HIGGS_INTERFERENCE_KAPPAQ",1.0);
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()==4) {
      if (((fl[0].IsGluon() && fl[1].IsGluon()) ||
	   (fl[0].IsQuark() && fl[1]==fl[0].Bar())) &&
	  fl[2].IsPhoton() && fl[3].IsPhoton()) {
	msg_Info()<<"!";
	return new Higgs_Virtual(pi,fl,mode,io,spin,kg,kq);
      }
    }
  }
  return NULL;
}
