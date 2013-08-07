#include "Higgs_Virtual.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"

#include "Wrappers.H"
#include "Ahiggs.h"
#include "Acont.h"

Complex ggH(int h1, int h2);
Complex Hgamgam(int h3, int h4);
Complex gggamgam(int h1, int h2, int h3, int h4);

using namespace HIGGS;
using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;

MODEL::Model_Base *HIGGS::Higgs_Virtual::s_model=NULL;

Higgs_Virtual::Higgs_Virtual(const Process_Info &pi,
			     const Flavour_Vector &flavs,
			     int mode,int io):
  Virtual_ME2_Base(pi,flavs), m_int(mode), m_io(io)
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

double Higgs_Virtual::htheta(const double &x)
{
  return x>0.0?1.0:0.0;
}

Complex Higgs_Virtual::lnrat(const double &x,const double &y)
{
  return Complex(log(dabs(x/y)),-M_PI*(htheta(-x)-htheta(-y)));
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
  s_bs=p_bs;
  p_bs->Setk0(11);
  p_bs->CalcEtaMu((Vec4D*)&p.front());
  double s=(p[2]+p[3]).Abs2(), rts=sqrt(s);
  Complex lo=0.0, nlo=0.0;
  Complex los=0.0, nlos=0.0;
  Complex lob=0.0, nlob=0.0;
  if (m_proc>1) {
    double deltar=1.0, cf=4.0/3.0;
    double s12=(p[2]+p[3]).Abs2();
    double s13=(p[2]-p[1]).Abs2();
    double s23=(p[3]-p[1]).Abs2();
    Complex l12=lnrat(-s12,mu_sq);
    Complex l13=lnrat(-s13,mu_sq);
    Complex l23=lnrat(-s23,mu_sq);
    lo=s13/s23+s23/s13;
    nlo=(-7.0-l12*l12)*lo
      +l13*(s23-2.0*s12)/s13+l23*(s13-2.0*s12)/s23
      +(s12*s12+s23*s23)/s13/s23*(sqr(l12-l23)+sqr(M_PI))
      +(s12*s12+s13*s13)/s13/s23*(sqr(l12-l13)+sqr(M_PI))
      -4.0*l12;
    double lmur=log(m_mur2/s);
    m_res.IR2()=-2.0*4.0/3.0;
    m_res.IR()=-(3.0+2.0*lmur)*4.0/3.0;
    m_res.Finite()=(nlo/lo).real()*4.0/3.0;
    double ft=sqr(-(4.0*M_PI*alpha0)*sqr(m_flavs[0].Charge()))/3.0;
    m_born=ft*lo.real();
    return;
  }
  Complex fslo=A_prod_1l(rts,muR)/s*A_dec_1l(rts,muR)/s/
    ((s-m_mh*m_mh)+I*m_mh*m_gh);
  Complex fsnlo=(A_prod_2l(rts,muR)/s*A_dec_1l(rts,muR)/s+
		 A_prod_1l(rts,muR)/s*A_dec_2l(rts,muR)/s)/
    ((s-m_mh*m_mh)+I*m_mh*m_gh)/(alpha_s(muR)/(2.0*M_PI));
  Complex fblo=-4.0*sumQsq*alpha0*alpha_s(muR);
  Complex bloset[5];
  // precalculated massive version for continuum bkg at LO
  // for (int i(0);i<5;++i) bloset[i]=-A_cont_1l(i+1,1.+2.*sij(1,4)/sij(1,2),sij(1,2),muR);
  for (int i(1);i>=-1;i-=2) {
    for (int j(1);j>=-1;j-=2) {
      for (int k(1);k>=-1;k-=2) {
	for (int l(1);l>=-1;l-=2) {
	  Complex clos(0.0), cnlos(0.0);
	  Complex clob(0.0), cnlob(0.0);
	  Complex clo(0.0), cnlo(0.0);
	  int hel;
          if (i==j && j==k && k==l) hel=1; //pppp
          else if(-i==-j && -j==k && k==l) hel=3; // mmpp
          else if(-i==j && j==-k && -k==l) hel=4; // mpmp
          else if(-i==j && j==k && k==-l) hel=5; // mppm
          else hel=2; //mppp, pmpp, ppmp, pppm
	  if (m_int&1) {
	    if (i==j && k==l) {
	      clos+=fslo*s*s;
	      clo+=fslo*s*s;
	    }
	    if (i==j && k==l) {
	      cnlos+=fsnlo*s*s;
	      cnlo+=fsnlo*s*s;
	    }
	  }
	  if (m_int&2) {
	    clob+=fblo*gggamgam1l(i,j,k,l);
	    cnlob+=fblo*gggamgam2l(i,j,k,l);
	    clo+=fblo*gggamgam1l(i,j,k,l);
	    cnlo+=fblo*gggamgam2l(i,j,k,l);
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
  double b0=(11.0*3.0-2.0*(Flavour(kf_quark).Size()/2))/6.0;
  double lmur=log(m_mur2/s);
  m_res.IR2()=-2.0*3.0;
  m_res.IR()=-2.0*(b0+3.0*lmur);
  if (m_io==1) {
    lo-=los+lob;
    nlo-=nlos+nlob;
  }
  if (m_io==2) {
    lo-=lob;
    nlo-=nlob;
  }
  m_res.Finite()=(nlo/lo).real()+3.0*sqr(M_PI)+2.0*b0*lmur;
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
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()==4) {
      if (((fl[0].IsGluon() && fl[1].IsGluon()) ||
	   (fl[0].IsQuark() && fl[1]==fl[0].Bar())) &&
	  fl[2].IsPhoton() && fl[3].IsPhoton()) {
	msg_Info()<<"!";
	return new Higgs_Virtual(pi,fl,mode,io);
      }
    }
  }
  return NULL;
}
