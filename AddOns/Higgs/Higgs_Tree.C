#include "Higgs_Tree.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include "Wrappers.H"
#include "dilog.h"
#include "Ahiggs.C"
#include "A_higgs.C"
#include "A_spin2.C"
#include "Acont.C"
#include "A_cont.C"

using namespace HIGGS;
using namespace AMEGIC;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

MODEL::Model_Base *HIGGS::Higgs_Tree::s_model=NULL;

Higgs_Tree::Higgs_Tree(const Process_Info &pi,
		       const Flavour_Vector &flavs,
		       int mode,int io,int spin):
  Tree_ME2_Base(pi,flavs), m_int(mode), m_io(io), m_spin(spin)
{
  m_mh=Flavour(kf_h0).Mass();
  m_gh=Flavour(kf_h0).Width();
  m_oqcd=m_oew=2;
  m_n=m_flavs.size();
  m_b=std::vector<int>(m_n,1);
  m_b[1]=m_b[0]=-1;
  p_bs = new Basic_Sfuncs
    (m_n,m_n,(Flavour*)&m_flavs.front(),&m_b.front());  
  p_bs->Initialize();
  m_ress.resize(1<<m_n,0.0);
  m_resb.resize(1<<m_n,0.0);
  m_resa.resize(1<<m_n,0.0);
  m_rest.resize(1<<m_n,0.0);
  m_namps=3;
  m_mode=0;
  if (m_n==4) {
    if (m_flavs[0].IsGluon() && m_flavs[1].IsGluon()) m_mode=1;
    if (m_flavs[0].IsQuark() && m_flavs[1]==m_flavs[0].Bar()) {
      if (m_flavs[0].IsAnti()) m_mode=5;
      else m_mode=4;
    }
  }
  else {
    if (m_flavs[0].IsGluon() && m_flavs[1].IsGluon() &&
	m_flavs[4].IsGluon()) m_mode=1;
    if (m_flavs[0].IsGluon() && m_flavs[1].IsQuark() &&
	m_flavs[4]==m_flavs[1]) m_mode=2;
    if (m_flavs[0].IsQuark() && m_flavs[1].IsGluon() &&
	m_flavs[4]==m_flavs[0]) m_mode=3;
    if (m_flavs[0].IsQuark() && m_flavs[4].IsGluon() && 
	m_flavs[1]==m_flavs[0].Bar() ) {
      if (m_flavs[0].IsAnti()) m_mode=5;
      else m_mode=4;
    }
    if (m_mode!=1 && (m_int&4)) m_namps=4;
  }
}

Higgs_Tree::~Higgs_Tree()
{
  delete p_bs;
}

double Higgs_Tree::Calc(const Vec4D_Vector &p)
{
  DEBUG_FUNC(this<<", m_mode = "<<m_mode);
  double muR=p_aqcd->Scale();
  if (muR>0.0) muR=sqrt(muR);
  else muR=sqrt(rpa->gen.CplScale());
  mu_sq=sqr(muR);
  alpha0=s_model->ScalarConstant("alpha_QED(0)");
  msg_Debugging()<<"\\mu_R = "<<muR<<" -> alpha = "
		 <<alpha0<<", alpha_s = "<<alpha_s(muR)<<"\n";
  for (size_t i(0);i<p.size();++i)
    msg_Debugging()<<"p["<<i<<"]="<<p[i]<<"\n";
  s_bs=p_bs;
  p_bs->Setk0(11);
  p_bs->CalcEtaMu((Vec4D*)&p.front());
  double s=(p[2]+p[3]).Abs2(), rts=sqrt(s);
  Complex fs=A_prod_1l(rts,muR)/s*A_dec_1l(rts,muR)/s/
    ((s-m_mh*m_mh)+I*m_mh*m_gh);
  Complex fb=-4.0*sumQsq*alpha0*alpha_s(muR);
  double qc=m_flavs[m_flavs[0].IsQuark()?0:1].Charge();
  Complex ft=-2.0*(4.0*M_PI*alpha0)*sqr(qc);
  if (m_n==5) {
    double cpl=2.0*sqrt(2.0*PI*alpha_s(muR));
    fs*=cpl;
    fb*=cpl;
    ft*=cpl;
    double res=0.0;
    size_t n(0);
    for (int i(1);i>=-1;i-=2) {
      for (int j(1);j>=-1;j-=2) {
	for (int k(1);k>=-1;k-=2) {
	  for (int l(1);l>=-1;l-=2) {
	    for (int m(1);m>=-1;m-=2) {
	      m_ress[n]=m_resb[n]=m_resa[n]=m_rest[n]=0.0;
	      if (m_mode==1) {
		if (m_int&1) {
		  if (m_spin==0) m_ress[n]+=fs*sqrt(3.0)*ggHg(i,j,m)*Hgamgam(k,l);
		  else m_ress[n]+=fs*sqrt(3.0)*ggXgamgamg(i,j,k,l,m);
		}
		// p_lab[0]=-Vec4D(-3066.22278407256,0,0,-3066.22278407256);
		// p_lab[1]=-Vec4D(-3392.650805216435,0,0,3392.650805216435);
		// p_lab[2]=Vec4D(2796.639070684271,-39.75655359831894,-42.22478302139858,-2796.037656367053);
		// p_lab[3]=Vec4D(597.1602026406933,39.25491506394217,23.44625670287839,-595.4071147608329);
		// p_lab[4]=Vec4D(3065.074315963385,0.5016385343768671,18.77852631852376,3065.016749984655);
		// original term is:43.42257556794053
		// 1st dipole / original is:1.000008474288235
		// 2nd dipole / original is:8.487253850919057e-06
		if (m_int&2) m_resb[n]+=fb*sqrt(3.0)*gggamgamg(i,j,k,l,m);
		// p_lab[0]=-Vec4D(-32.23418530748044,0.,0.,-32.23418530748044);
		// p_lab[1]=-Vec4D(-1660.706232403784,0.,0.,1660.706232403784);
		// p_lab[2]=Vec4D(712.2231111636451,62.25933681401244,-0.236953896794764,-709.4966377002445);
		// p_lab[3]=Vec4D(950.8479315163148,-61.77272159808692,-0.06968860974897866,-948.8392460679088);
		// p_lab[4]=Vec4D(29.86937503128269,-0.4866152159255259,0.3066425065437467,29.86383667187162);
		// original term is:0.1068900378010579
		// 1st dipole / original is:0.9992348714535241
		// 2nd dipole / original is:1.796378435464722e-06
		// alpha_s=0.1146279368135488, 1/alpha=137.036
	      }
	      if (m_mode==2) {
		if ((m_int&1) && j!=m) {
		  if (m_spin==0) m_ress[n]+=fs*gqHq(i,j)*Hgamgam(k,l)/sqrt(2.0);
		  else m_ress[n]+=fs*gqXgamgamq(i,j,k,l)/sqrt(2.0);
		}
		if ((m_int&2) && j!=m) m_resb[n]+=fb*gqgamgamq(i,j,k,l)/sqrt(2.0);
		if ((m_int&4) && j!=m) m_rest[n]+=ft*gqgamgamq_tree(i,j,k,l)/sqrt(2.0);
	      }
	      if (m_mode==3) {
		if ((m_int&1) && i!=m) {
		  if (m_spin==0) m_ress[n]+=fs*qgHq(i,j)*Hgamgam(k,l)/sqrt(2.0);
		  else m_ress[n]+=fs*qgXgamgamq(i,j,k,l)/sqrt(2.0);
		}
		// p_lab[0]=-Vec4D(-32.11250970836028,0,0,-32.11250970836028)
		// p_lab[1]=-Vec4D(-1659.561441954768,0,0,1659.561441954768)
		// p_lab[2]=Vec4D(711.7248741577324,62.05695068106898,-0.236335893196161,-709.0142280041727)
		// p_lab[3]=Vec4D(950.1876162610141,-61.57134343344105,-0.0696714381307442,-948.1906300472918)
		// p_lab[4]=Vec4D(29.76146124438999,-0.4856072476279195,0.3060073313269036,29.75592580504852)
		// original term is:7.42145225688603
		// 1st dipole / original is:1.000001662922538
		if ((m_int&2) && i!=m) m_resb[n]+=fb*qggamgamq(i,j,k,l)/sqrt(2.0);
		// p_lab[0]=-Vec4D(-501.2244608524589,0.,0.,-501.2244608524589);
		// p_lab[1]=-Vec4D(-15.8309789562357,0.,0.,15.8309789562357);
		// p_lab[2]=Vec4D(91.72724351622065,-7.686408643444808,58.82278327525236,69.96203608388414 );
		// p_lab[3]=Vec4D(184.6125678122652,7.626009750358896,-59.14347007428552,174.7160385227921 );
		// p_lab[4]=Vec4D(240.715628480209,0.06039889308591277,0.3206867990331586,240.7154072895473 );
		// original term is:0.001221596879623276
		// 1st dipole / original is:0.9974613722643826
		// alpha_s=0.1146279368135488, 1/alpha=137.036
		if ((m_int&4) && i!=m) m_rest[n]+=ft*qggamgamq_tree(i,j,k,l)/sqrt(2.0);
	      }
	      if (m_mode==4) {
		if ((m_int&1) && i!=j) {
		  if (m_spin==0) m_ress[n]+=fs*qqbHg(i,m)*Hgamgam(k,l)/sqrt(2.0);
		  else m_ress[n]+=fs*qqbXgamgamg(i,k,l,m)/sqrt(2.0);
		}
		if ((m_int&2) && i!=j) m_resb[n]+=fb*qqbgamgamg(i,k,l,m)/sqrt(2.0);
		if ((m_int&4) && i!=j) m_rest[n]+=ft*qqbgamgamg_tree(i,k,l,m)/sqrt(2.0);
	      }
	      if (m_mode==5) {
		if ((m_int&1) && j!=i) {
		  if (m_spin==0) m_ress[n]+=fs*qbqHg(j,m)*Hgamgam(k,l)/sqrt(2.0);
		  else m_ress[n]+=fs*qbqXgamgamg(j,k,l,m)/sqrt(2.0);
		}
		if ((m_int&2) && j!=i) m_resb[n]+=fb*qbqgamgamg(j,k,l,m)/sqrt(2.0);
		if ((m_int&4) && j!=i) m_rest[n]+=ft*qbqgamgamg_tree(j,k,l,m)/sqrt(2.0);
	      }
	      res+=2.0*(m_ress[n]*std::conj(m_rest[n])).real();
	      m_resa[n]=m_ress[n]+m_resb[n];
	      if (m_io==0) m_ress[n]=m_resb[n]=0.0;
	      if (m_io==2) m_ress[n]=m_rest[n]=0.0;
	      if (m_io==1) m_rest[n]=0.0;
	      res+=(m_resa[n]*std::conj(m_resa[n])).real();
	      res-=(m_ress[n]*std::conj(m_ress[n])).real();
	      res-=(m_resb[n]*std::conj(m_resb[n])).real();
	      res+=(m_rest[n]*std::conj(m_rest[n])).real();
	      ++n;
	    }
	  }
	}
      }
    }
    return res*8.0;
  }
  double res(0.0); 
  size_t n(0);
  for (int i(1);i>=-1;i-=2) {
    for (int j(1);j>=-1;j-=2) {
      for (int k(1);k>=-1;k-=2) {
	for (int l(1);l>=-1;l-=2) {
	  m_ress[n]=m_resb[n]=m_resa[n]=m_rest[n]=0.0;
	  if (m_mode==1) {
	    if (m_int&1) {
	      if (m_spin==0) m_ress[n]+=fs*ggH(i,j)*Hgamgam(k,l);
	      else m_ress[n]+=fs*ggXgamgam(i,j,k,l);
	    }
	    if (m_int&2) m_resb[n]+=fb*gggamgam(i,j,k,l);
	  }
	  if (m_mode==4) {
	    if (m_int&1) {
	      if (m_spin!=0 && i!=j) m_ress[n]+=fs*qqbXgamgam(i,k,l);
	    }
	    if ((m_int&4) && i!=j) m_resb[n]+=ft*qqbgamgam_tree(i,k,l);
	  }
	  if (m_mode==5) {
	    if (m_int&1) {
	      if (m_spin!=0 && i!=j) m_ress[n]+=fs*qbqXgamgam(j,k,l);
	    }
	    if ((m_int&4) && j!=i) m_resb[n]+=ft*qbqgamgam_tree(j,k,l);
	  }
	  m_resa[n]=m_ress[n]+m_resb[n];
	  if (m_io==0) m_ress[n]=m_resb[n]=0.0;
	  if (m_io==2) m_ress[n]=0.0;
	  res+=(m_resa[n]*std::conj(m_resa[n])).real();
	  res-=(m_ress[n]*std::conj(m_ress[n])).real();
	  res-=(m_resb[n]*std::conj(m_resb[n])).real();
	  ++n;
	}
      }
    }
  }
  if (m_mode==4 || m_mode==5) return res*3.0;
  return res*8.0;
}

std::vector<Complex> Higgs_Tree::GetAmplitudes(const size_t &id)
{
  if (id==3) return m_rest;
  if (id==2) return m_resb;
  if (id==1) return m_ress;
  return m_resa;
}

Complex Higgs_Tree::GetPhase(const size_t &id)
{
  if (id==3) return Complex(1.0,0.0);
  if (id==2) return Complex(-1.0,0.0);
  if (id==1) return Complex(-1.0,0.0);
  return Complex(1.0,0.0);
}

void Higgs_Tree::FillCombinations
(std::set<std::pair<size_t,size_t> > &combs,
 std::map<size_t,ATOOLS::Flavour_Vector> &fls)
{
  DEBUG_FUNC("n="<<m_n<<",mode="<<m_mode<<",int="<<m_int);
  if (m_n==4) {
    if (m_mode>=4) {
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<2));
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<3));
      fls[(1<<1)|(1<<3)].push_back(Flavour(m_flavs[1].Bar()));
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<3));
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<2));
      fls[(1<<1)|(1<<2)].push_back(Flavour(m_flavs[1].Bar()));
    }
    else {
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<1));
      combs.insert(std::pair<size_t,size_t>(1<<2,1<<3));
      fls[(1<<2)|(1<<3)].push_back(Flavour(kf_h0));
    }
  }
  else {
    if (m_mode==1) {
      combs.insert(std::pair<size_t,size_t>(1<<2,1<<3));
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<1));
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<4));
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<4));
      combs.insert(std::pair<size_t,size_t>((1<<2)|(1<<3),1<<0));
      combs.insert(std::pair<size_t,size_t>((1<<2)|(1<<3),1<<1));
      combs.insert(std::pair<size_t,size_t>((1<<2)|(1<<3),1<<4));
      fls[(1<<2)|(1<<3)].push_back(Flavour(kf_h0));
      fls[(1<<0)|(1<<1)].push_back(Flavour(kf_gluon));
      fls[(1<<0)|(1<<4)].push_back(Flavour(kf_gluon));
      fls[(1<<1)|(1<<4)].push_back(Flavour(kf_gluon));
    }
    if (m_mode==2) {
      combs.insert(std::pair<size_t,size_t>(1<<2,1<<3));
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<4));
      combs.insert(std::pair<size_t,size_t>((1<<2)|(1<<3),1<<0));
      fls[(1<<2)|(1<<3)].push_back(Flavour(kf_h0));
      fls[(1<<1)|(1<<4)].push_back(Flavour(kf_gluon));
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<2));
      fls[(1<<1)|(1<<2)].push_back(m_flavs[1].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<3));
      fls[(1<<1)|(1<<3)].push_back(m_flavs[1].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<0));
      fls[(1<<1)|(1<<0)].push_back(m_flavs[1].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<4,1<<2));
      fls[(1<<4)|(1<<2)].push_back(m_flavs[4]);
      combs.insert(std::pair<size_t,size_t>(1<<4,1<<3));
      fls[(1<<4)|(1<<3)].push_back(m_flavs[4]);
      combs.insert(std::pair<size_t,size_t>(1<<4,1<<0));
      fls[(1<<4)|(1<<0)].push_back(m_flavs[4]);
      combs.insert(std::pair<size_t,size_t>((1<<1)|(1<<2),1<<3));
      combs.insert(std::pair<size_t,size_t>((1<<1)|(1<<2),1<<0));
      combs.insert(std::pair<size_t,size_t>((1<<1)|(1<<3),1<<2));
      combs.insert(std::pair<size_t,size_t>((1<<1)|(1<<3),1<<0));
      combs.insert(std::pair<size_t,size_t>((1<<1)|(1<<0),1<<2));
      combs.insert(std::pair<size_t,size_t>((1<<1)|(1<<0),1<<3));
    }
    if (m_mode==3) {
      combs.insert(std::pair<size_t,size_t>(1<<2,1<<3));
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<4));
      combs.insert(std::pair<size_t,size_t>((1<<2)|(1<<3),1<<1));
      fls[(1<<2)|(1<<3)].push_back(Flavour(kf_h0));
      fls[(1<<0)|(1<<4)].push_back(Flavour(kf_gluon));
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<2));
      fls[(1<<0)|(1<<2)].push_back(m_flavs[0].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<3));
      fls[(1<<0)|(1<<3)].push_back(m_flavs[0].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<1));
      fls[(1<<0)|(1<<1)].push_back(m_flavs[0].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<4,1<<2));
      fls[(1<<4)|(1<<2)].push_back(m_flavs[4]);
      combs.insert(std::pair<size_t,size_t>(1<<4,1<<3));
      fls[(1<<4)|(1<<3)].push_back(m_flavs[4]);
      combs.insert(std::pair<size_t,size_t>(1<<4,1<<1));
      fls[(1<<4)|(1<<1)].push_back(m_flavs[4]);
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<2),1<<3));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<2),1<<1));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<3),1<<2));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<3),1<<1));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<0),1<<2));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<0),1<<3));
    }
    if (m_mode>=4) {
      combs.insert(std::pair<size_t,size_t>(1<<2,1<<3));
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<1));
      combs.insert(std::pair<size_t,size_t>((1<<2)|(1<<3),1<<4));
      fls[(1<<2)|(1<<3)].push_back(Flavour(kf_h0));
      fls[(1<<0)|(1<<1)].push_back(Flavour(kf_gluon));
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<2));
      fls[(1<<0)|(1<<2)].push_back(m_flavs[0].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<3));
      fls[(1<<0)|(1<<3)].push_back(m_flavs[0].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<0,1<<4));
      fls[(1<<0)|(1<<4)].push_back(m_flavs[0].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<2));
      fls[(1<<1)|(1<<2)].push_back(m_flavs[1].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<3));
      fls[(1<<1)|(1<<3)].push_back(m_flavs[1].Bar());
      combs.insert(std::pair<size_t,size_t>(1<<1,1<<4));
      fls[(1<<1)|(1<<4)].push_back(m_flavs[1].Bar());
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<2),1<<3));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<2),1<<4));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<3),1<<2));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<3),1<<4));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<4),1<<2));
      combs.insert(std::pair<size_t,size_t>((1<<0)|(1<<4),1<<3));
    }
    if (m_flavs[0].IsQuark() && m_flavs[4].IsGluon() && 
	m_flavs[1]==m_flavs[0].Bar() ) {
      if (m_flavs[0].IsAnti()) m_mode=5;
      else m_mode=4;
    }
  }
  std::set<std::pair<size_t,size_t> > tcombs(combs);
  for (std::set<std::pair<size_t,size_t> >::const_iterator
	 cit=tcombs.begin();cit!=tcombs.end();++cit) {
    size_t ida(cit->first), idb(cit->second);
    size_t idc((1<<m_n)-1-ida-idb);
    msg_Debugging()<<"comb "<<ID(ida)
		   <<" "<<ID(idb)<<" "<<ID(idc)<<"\n";
    combs.insert(std::pair<size_t,size_t>(idb,ida));
    combs.insert(std::pair<size_t,size_t>(idb,idc));
    combs.insert(std::pair<size_t,size_t>(idc,idb));
    combs.insert(std::pair<size_t,size_t>(idc,ida));
    combs.insert(std::pair<size_t,size_t>(ida,idc));
  }
  std::map<size_t,ATOOLS::Flavour_Vector> tfls(fls);
  for (std::map<size_t,ATOOLS::Flavour_Vector>::const_iterator
	 cit=tfls.begin();cit!=tfls.end();++cit) {
    msg_Debugging()<<"flav "<<ID(cit->first)
		   <<" -> "<<cit->second<<"\n";
    Flavour_Vector fl(cit->second);
    for (size_t i(0);i<fl.size();++i) fl[i]=fl[i].Bar();
    fls[(1<<m_n)-1-cit->first]=fl;
  }
}

int Higgs_Tree::OrderQCD(const int &id)
{
  return m_n==4?2:3;
}

int Higgs_Tree::OrderEW(const int &id)
{
  return 2;
}

double Higgs_Tree::TR() const
{
  return 1.0;
}

Complex Higgs_Tree::GetHelicityPhase(const Vec4D &pijt,const Vec4D &eps1)
{
  Vec4D k(Vec4D(1.0,0.0,0.0,0.0)+p_bs->Getk1());
  Vec4D p(pijt), eps(eps1+0.5/(eps1*k)*k);
  Complex eip=-std::conj(p_bs->CalcS(k,eps))*p_bs->CalcS(eps,p);
  eip/=-std::conj(p_bs->CalcS(k,p));
  return eip;
}

std::vector<Tree_ME2_Base::Map_Info> Higgs_Tree::GetFlavourHelicityMap()
{
  if (m_hmap.size()) return m_hmap;
  std::vector<int> perm(4);
  for (int i(0);i<4;++i) perm[i]=i;
  size_t n(0);
  m_hmap.resize(1<<4);
  for (int i(1);i>=-1;i-=2) {
    for (int j(1);j>=-1;j-=2) {
      for (int k(1);k>=-1;k-=2) {
	for (int l(1);l>=-1;l-=2) {
	  std::vector<int> hels;
	  hels.push_back(i);
	  hels.push_back(j);
	  hels.push_back(k);
	  hels.push_back(l);
	  m_hmap[n].m_perm=perm;
	  m_hmap[n].m_hels=hels;
	  m_hmap[n].m_id=n++;
	}
      }
    }
  }
  return m_hmap;
}

DECLARE_TREEME2_GETTER(Higgs_Tree,"Higgs_Tree")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,Higgs_Tree>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Higgs") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype==nlo_type::lo ||
      pi.m_fi.m_nloqcdtype==nlo_type::born ||
      pi.m_fi.m_nloqcdtype==nlo_type::real) {
    Data_Reader read(" ",";","#","=");
    int io=read.GetValue<int>("HIGGS_INTERFERENCE_ONLY",0);
    int mode=read.GetValue<int>("HIGGS_INTERFERENCE_MODE",7);
    int spin=read.GetValue<int>("HIGGS_INTERFERENCE_SPIN",0);
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()==4) {
      if (fl[2].IsPhoton() && fl[3].IsPhoton()) {
	if ((fl[0].IsGluon() && fl[1].IsGluon()) ||
	    ((mode&4) && fl[0].IsQuark() && fl[1]==fl[0].Bar())) {
	  msg_Info()<<"!";
	  return new Higgs_Tree(pi,fl,mode,io,spin);
	}
      }
    }
    if (fl.size()==5 &&
	fl[2].IsPhoton() && fl[3].IsPhoton()) {
      if ((fl[0].IsGluon() && fl[1].IsGluon() && fl[4].IsGluon()) ||
	  (fl[0].IsGluon() && fl[1].IsQuark() && fl[4]==fl[1]) ||
	  (fl[0].IsQuark() && fl[1].IsGluon() && fl[4]==fl[0]) ||
	  (fl[0].IsQuark() && fl[1]==fl[0].Bar() && fl[4].IsGluon())) {
	msg_Info()<<"!";
	return new Higgs_Tree(pi,fl,mode,io,spin);
      }
    }
  }
  return NULL;
}
