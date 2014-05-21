#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;
using namespace AMEGIC;
using namespace MODEL;
using namespace std;

#define SQRT_05 0.70710678118654757

DipoleSplitting_Base::DipoleSplitting_Base() 
{
  m_type  = dpt::none;
  m_spfdef = 0.0;
  p_cpl = NULL;
  m_av=sqrt(-1.0);
  m_a=-1.0;
  
  Flavour hfl(kf_quark);
  m_nf = hfl.Size()/2;

  CSC.Init();
  m_g1 = 1.5*CSC.CF;
  m_g2 = 11./6.*CSC.CA-2./3.*CSC.TR*m_nf;

  m_pfactors.clear();m_pfactors.push_back(1.);
  m_i=-1;
  m_j=-1;
  m_k=-1;
  m_tij=-1;
  m_tk=-1;
  m_m=0;
  m_alpha=1.;
  m_kt2max=std::numeric_limits<double>::max();
  m_amin = max(ATOOLS::Accu(),1.e-8);
  m_k0sqf=ToType<double>(rpa->gen.Variable("CSS_FS_PT2MIN"));
  m_k0sqi=ToType<double>(rpa->gen.Variable("CSS_IS_PT2MIN"));
  m_es=ToType<int>(rpa->gen.Variable("CSS_EVOLUTION_SCHEME"));
  double helpd;
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.SetInputPath(rpa->GetPath());
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));

  if (reader.ReadFromFile(helpd,"DIPOLE_AMIN")) {
    m_amin = helpd;
    msg_Tracking()<<"Set dipole cut alphamin="<<m_amin<<"."<<std::endl;
  }

  m_kappa=2./3.;
  if (reader.ReadFromFile(helpd,"DIPOLE_KAPPA")) {
    m_kappa = helpd;
    msg_Tracking()<<"Set massive dipole kappa="<<m_kappa<<"."<<std::endl;
  }
}

void DipoleSplitting_Base::SetCoupling(const MODEL::Coupling_Map *cpls)
{
  if (cpls->find("Alpha_QCD")!=cpls->end()) p_cpl=cpls->find("Alpha_QCD")->second;
  else THROW(fatal_error,"Coupling not found");
  msg_Tracking()<<"DipoleSplitting_Base:: alpha = "<<*p_cpl<<endl;
  m_spfdef = -8.*M_PI*p_cpl->Default();
}

void DipoleSplitting_Base::SetMomenta(const Vec4D* mom)
{
  m_mom.clear();
  for(int i=0;i<=m_m;i++) m_mom.push_back(mom[i]);
}

void DipoleSplitting_Base::CalcVectors(ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double B)
{
  m_dpollist.clear();
  m_pfactors.clear();
  if (1) {
    Vec3D pv(p2);
    Vec3D ptp=Vec3D(p1)-(p1[0]/p2[0])*pv;
    Vec3D ptt=cross(ptp,pv);

    m_dpollist.push_back(Vec4D(0.,ptt/ptt.Abs())); m_pfactors.push_back(1.);

    Vec4D vh(0.,ptp/ptp.Abs());
    m_dpollist.push_back(vh); m_pfactors.push_back((B-1.)/B);

  }
  else {
  Vec4D vh1,vh2,p=p1;//+p2;
  
  double ps=sqrt(sqr(p[1])+sqr(p[2])+sqr(p[3]));
  double pt=sqrt(sqr(p[1])+sqr(p[2]));
  if(!ATOOLS::IsZero(pt)){
    vh1 = Vec4D(0.,p[1]*p[3]/ps/pt,p[2]*p[3]/ps/pt,-pt/ps);
    vh2 = Vec4D(0.,-p[2]/pt,p[1]/pt,0.);
    if(p[1]+p[2]<0)vh2=-1.*vh2;
  }
  else {
     vh1 = SQRT_05*Vec4D(0.,1.,-1.,0.);
     vh2 = SQRT_05*Vec4D(0.,1.,1.,0.);
  }
  m_dpollist.push_back(vh1); m_pfactors.push_back(1.);
  m_dpollist.push_back(vh2); m_pfactors.push_back(1.);

   int sgn=1;
   double ph=p.Abs2();//2.*(p1*p2);
   if(!ATOOLS::IsZero(ps)){
     vh1=1./sqrt(dabs(ph))*Vec4D(ps,p[0]*p[1]/ps,p[0]*p[2]/ps,p[0]*p[3]/ps);
     if (ph<0.) sgn=-1;
   }
   else vh1=Vec4D(0.,0.,0.,1.);
   m_dpollist.push_back(vh1); m_pfactors.push_back(sgn);
   
   ph = (1.-B)/B/ph;

   if (ph>=0.) sgn = 1;
   else sgn = -1;
   vh2 = sqrt(dabs(ph))*p;
   m_dpollist.push_back(vh2); 
   m_pfactors.push_back(sgn);
  }
}

double DipoleSplitting_Base::GetR(const ATOOLS::Vec4D* mom,const ATOOLS::Vec4D* LOmom)
{
  double ptijk=2.*(m_ptij*m_pj)*(m_ptk*m_pj)/(m_ptij*m_ptk);
  double spt=0.;
  for(int a=2;a<m_m;a++) {
    for(int b=a+1;b<m_m;b++) {
      for(int c=2;c<m_m;c++)
	if(c!=a&&c!=b) spt+=0.5*sqr((LOmom[a]+LOmom[b])*LOmom[c])/
			 ((LOmom[a]*LOmom[b])*(LOmom[b]*LOmom[c])*(LOmom[c]*LOmom[a]));
    }
  }
  return 1./(1+ptijk*spt);
}



double DipoleSplitting_Base::Vijf(int type,int mode)
{
  double loga = log(m_alpha);
  switch (type) {
  case 1:
  case 2:
    return 5.-0.5*sqr(M_PI)-sqr(loga)+1.5*(m_alpha-1-loga)-(mode?0.5:0.);
  case 3:
    return -16./9.-2./3.*(m_alpha-1-loga);
  case 4:
    return 50./9.-0.5*sqr(M_PI)-sqr(loga)+11./6.*(m_alpha-1-loga)-(mode?1./6.:0.);
  }
  return 0.;
}

double DipoleSplitting_Base::Vif(int type,int mode)
{
  switch (type) {
  case 1:
  case 2:
    return Vijf(1,mode);
  case 3:
  case 4:
    return Vijf(4,mode)+m_nf*CSC.TR/CSC.CA*Vijf(3,mode);
  }
  return 0.;
}

double DipoleSplitting_Base::Vie1(int type)
{
  switch (type) {
  case 1:
  case 2:
    return 1.5;
  case 3:
  case 4:
    return m_g2/CSC.CA;
  }
  return 0.;
}

double DipoleSplitting_Base::Vie2(int type)
{
  return 1.;
}

bool DipoleSplitting_Base::Reject(const double &alpha)
{
  if (IsBad(m_av))
    msg_Error()<<METHOD<<"(): Average is "<<m_av<<" in "
	       <<Demangle(typeid(*this).name())
	       <<"[type="<<m_ft<<"]"<<std::endl;
  if (m_mcmode==1) {
    int da(m_av>0.0 && (m_kt2<m_kt2max || IsEqual(m_kt2,m_kt2max,1.0e-6))), ds(alpha<=m_alpha);
    m_mcsign=ds-da;
    return m_mcsign==0;
  }
  if (m_mcmode==2) {
    m_mcsign=m_av>0.0 && (m_kt2<m_kt2max || IsEqual(m_kt2,m_kt2max,1.0e-6));
    return m_mcsign==0;
  }
  return alpha>m_alpha || m_kt2>m_kt2max;
}

double DipoleSplitting_Base::GetF()
{
  if (Reject(m_a)) return 0.;
  else return GetValue();
}

double DipoleSplitting_Base::GetValue()
{
  THROW(fatal_error, "Virtual function not reimplemented.");
  return 0.0;
}
