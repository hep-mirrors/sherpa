#include "Hadronisation_Parameters.H"
#include "Flavour.H"
#include "MathTools.H"
#include "Data_Read.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

bool Momenta_Stretcher::MassThem(const int n,Vec4D * momenta,const double * masses)
{
  if (n==2) {
    Vec4D cms         = momenta[0]+momenta[1];
    Poincare boost(cms);
    for (int i=0;i<2;i++) boost.Boost(momenta[i]);
    double energy     = momenta[0][0]+momenta[1][0];
    if (masses[0]+masses[1]<energy) {
      double m12      = sqr(masses[0]);
      double m22      = sqr(masses[1]);
      double energy0  = (sqr(energy)+m12-m22)/(2.*energy);
      double energy1  = (sqr(energy)-m12+m22)/(2.*energy);
      Vec3D direction = Vec3D(momenta[0])/(Vec3D(momenta[0]).Abs());
      Vec3D p0        = direction*sqrt(sqr(energy0)-m12);
      Vec3D p1        = (-1.)*p0;
      momenta[0]      = Vec4D(energy0,p0);
      momenta[1]      = Vec4D(energy1,p1);
      for (int i=0;i<2;i++) boost.BoostBack(momenta[i]);
      return true; 
    }
    else {
      for (int i=0;i<2;i++) boost.BoostBack(momenta[i]);
      return false; 
    }
  }
  else {
    double xmt         = 0.;
    double * oldens2   = new double[n];
    double * ens       = new double[n];
    Vec4D cms          = Vec4D(0.,0.,0.,0.);
    for (short int k=0;k<n;k++) {
      xmt       += masses[k];
      cms       += momenta[k];
      oldens2[k] = sqr(momenta[k][0]);
    }
    if (cms[0]>xmt) {
      double ET  = sqrt(cms.Abs2()); 
      double x   = sqrt(1.-sqr(xmt/ET));
      double acc = ET*1.e-14;
      
      double f0,g0,x2;
      for (int i=0;i<10;i++) {
	f0 = -ET;g0 = 0.;x2 = x*x;
	for (short int k=0;k<n;k++) {
	  ens[k] = sqrt(sqr(masses[k])+x2*oldens2[k]);
	  f0    += ens[k];
	  g0    += oldens2[k]/ens[k];
	}
	if (dabs(f0)<acc) break; 
	x -= f0/(x*g0);  
      }
      for (short int k=0;k<n;k++) {
	momenta[k] = Vec4D(ens[k],x*Vec3D(momenta[k]));
      }
      delete [] oldens2;
      delete [] ens;
      return true;
    }
    delete [] oldens2;
    delete [] ens;
    msg.Error()<<"ERROR in Momenta_Stretcher::StretchThem: "<<endl
	       <<"   Not enough energy ("<<cms<<") for the "<<n<<" masses ("<<xmt<<"); return false"<<endl
	       <<"   Masses :";
    for (int i=0;i<n-1;i++) msg.Error()<<masses[i]<<", ";msg.Error()<<masses[n-1]<<"."<<endl;
    return false;
    abort();
  }
  return false;
}

bool Momenta_Stretcher::ZeroThem(const int n,Vec4D * momenta)
{
  if (n==2) {
    double energy   = momenta[0][0]+momenta[1][0];
    Vec3D direction = Vec3D(momenta[0])/(Vec3D(momenta[0]).Abs());
    momenta[0]      = energy/2.*Vec4D(1.,direction);
    momenta[1]      = energy/2.*Vec4D(1.,-1.*direction);
    return true; 
  }
  else {
    double xmt         = 0.;
    double * oldps2    = new double[n];
    double * ens       = new double[n];
    Vec4D cms          = Vec4D(0.,0.,0.,0.);
    for (short int i=0;i<n;i++) {
      xmt      += sqrt(Max(0.,momenta[i].Abs2()));
      oldps2[i] = sqr(Vec3D(momenta[i]).Abs());
      cms       += momenta[i];
    }
    double ET  = sqrt(cms.Abs2()); 
    double x   = 1./sqrt(1.-sqr(xmt/ET));
    double acc = ET*1.e-14;
    xmt        = 0.;

    double f0,g0,x2;    
    short int iter = 0; 
    for (int i=0;i<10;i++) {
      f0 = -ET;g0 = 0.;x2 = x*x;
      for (short int i=0;i<n;i++) {
	ens[i] = sqrt(x2*oldps2[i]);
	f0    += ens[i];
	g0    += oldps2[i]/ens[i];
      }
      if (dabs(f0)<acc) break; 
      x -= f0/(x*g0);  
    }
    for (short int k=0;k<n;k++) {
      momenta[k] = Vec4D(ens[k],x*Vec3D(momenta[k]));
    }
    delete [] oldps2;
    delete [] ens;
    return true;
  }
  return false;
}


Hadronisation_Parameters AHADIC::hadpars;


Hadronisation_Parameters::Hadronisation_Parameters() :
  p_constituents(NULL),p_multiplets(NULL),p_transitions1(NULL)
{ }

Hadronisation_Parameters::~Hadronisation_Parameters() {
  if (p_constituents!=NULL) { delete p_constituents; p_constituents=NULL; }
  if (p_multiplets!=NULL)   { delete p_multiplets;   p_multiplets=NULL;   }
  if (p_transitions1!=NULL) { delete p_transitions1; p_transitions1=NULL; }
  if (p_transitions2!=NULL) { delete p_transitions2; p_transitions2=NULL; }
}

void Hadronisation_Parameters::Init(string dir,string file)
{
  msg.Tracking()<<"In Hadronisation_Parameters::Init("<<dir<<file<<")"<<endl;
  ReadParameters(dir,file);
  p_constituents = new Constituents(false);
  //if (msg.LevelIsTracking()) 
  p_constituents->PrintConstituents();

  p_multiplets   = new All_Hadron_Multiplets();
  if (msg.LevelIsTracking()) p_multiplets->PrintWaveFunctions(); 

  p_transitions1 = new All_Single_Transitions(p_multiplets);
  if (msg.LevelIsTracking()) p_transitions1->PrintSingleTransitions(); 

  p_transitions2 = new All_Double_Transitions(p_multiplets);
  if (msg.LevelIsTracking()) p_transitions2->PrintDoubleTransitions(); 
}
  
void Hadronisation_Parameters::ReadParameters(string dir,string file)
{
  Data_Read dataread(dir+file);
  m_parametermap[string("Strange_supression")] =
    dataread.GetValue<double>("STRANGE_SUPRESSION",0.2);      
  m_parametermap[string("Baryon_supression")]  = 
    dataread.GetValue<double>("BARYON_SUPRESSION",0.2);
  m_parametermap[string("Q_breakup")]          = 
    dataread.GetValue<double>("Q_BREAKUP",1.);
  m_parametermap[string("Offset")] =
    dataread.GetValue<double>("TransitionOffset",0.75);      
  m_parametermap[string("AngularSmearing")]    = 
    dataread.GetValue<double>("AngularSmearing",Get("Q_breakup"));
  m_parametermap[string("Max_Prod_Mass")]      = 
    dataread.GetValue<double>("Max_Prod_Mass",Get("Q_breakup"));
  m_parametermap[string("P_qs_by_P_qq")]       = 
    dataread.GetValue<double>("P_{QS}/P_{QQ}",Get("Strange_supression"));
  m_parametermap[string("P_ss_by_P_qq")]       = 
    dataread.GetValue<double>("P_{SS}/P_{QQ}",sqr(Get("Strange_supression")));    
  m_parametermap[string("P_di_1_by_P_di_0")]   = 
    dataread.GetValue<double>("P_{QQ_1}/P_{QQ_0}",1.);
  m_parametermap[string("FourQ")]          = 
    dataread.GetValue<double>("Four_Q_Cluster_Treatment",1.);
  m_parametermap[string("Mass_glue")]          = 
    dataread.GetValue<double>("M_GLUE",0.75);
  m_parametermap[string("Mass_down")]          = 
    dataread.GetValue<double>("M_DOWN",0.32);
  m_parametermap[string("Mass_up")]            = 
    dataread.GetValue<double>("M_UP",0.32);
  m_parametermap[string("Mass_strange")]       = 
    dataread.GetValue<double>("M_STRANGE",0.45);
  m_parametermap[string("Mass_charm")]         = 
    dataread.GetValue<double>("M_CHARM",1.8);
  m_parametermap[string("Mass_bottom")]        = 
    dataread.GetValue<double>("M_BOTTOM",5.2);
  m_parametermap[string("Mass_dd1")]           = 
    dataread.GetValue<double>("M_DD_1",2.*Get("Mass_down"));
  m_parametermap[string("Mass_ud0")]           = 
    dataread.GetValue<double>("M_UD_0",Get("Mass_up")+Get("Mass_down"));
  m_parametermap[string("Mass_ud1")]           = 
    dataread.GetValue<double>("M_UD_1",Get("Mass_up")+Get("Mass_down"));
  m_parametermap[string("Mass_uu1")]           = 
    dataread.GetValue<double>("M_UU_1",2.*Get("Mass_up"));
  m_parametermap[string("Mass_sd0")]           = 
    dataread.GetValue<double>("M_SD_0",Get("Mass_strange")+Get("Mass_down"));
  m_parametermap[string("Mass_sd1")]           = 
    dataread.GetValue<double>("M_SD_1",Get("Mass_strange")+Get("Mass_down"));
  m_parametermap[string("Mass_su0")]           = 
    dataread.GetValue<double>("M_SU_0",Get("Mass_strange")+Get("Mass_up"));
  m_parametermap[string("Mass_su1")]           = 
    dataread.GetValue<double>("M_SU_1",Get("Mass_strange")+Get("Mass_up"));
  m_parametermap[string("Mass_ss1")]           = 
    dataread.GetValue<double>("M_SS_1",2.*Get("Mass_strange"));
  m_parametermap[string("Mass_cd0")]           = 
    dataread.GetValue<double>("M_CD_0",Get("Mass_charm")+Get("Mass_down"));
  m_parametermap[string("Mass_cd1")]           = 
    dataread.GetValue<double>("M_CD_1",Get("Mass_charm")+Get("Mass_down"));
  m_parametermap[string("Mass_cu0")]           = 
    dataread.GetValue<double>("M_CU_0",Get("Mass_charm")+Get("Mass_up"));
  m_parametermap[string("Mass_cu1")]           = 
    dataread.GetValue<double>("M_CU_1",Get("Mass_charm")+Get("Mass_up"));
  m_parametermap[string("Mass_cs0")]           = 
    dataread.GetValue<double>("M_CS_0",Get("Mass_charm")+Get("Mass_strange"));
  m_parametermap[string("Mass_cs1")]           = 
    dataread.GetValue<double>("M_CS_1",Get("Mass_charm")+Get("Mass_strange"));
  m_parametermap[string("Mass_cc1")]           = 
    dataread.GetValue<double>("M_CC_1",2.*Get("Mass_charm"));
  m_parametermap[string("Mass_bd0")]           = 
    dataread.GetValue<double>("M_BD_0",Get("Mass_bottom")+Get("Mass_down"));
  m_parametermap[string("Mass_bd1")]           = 
    dataread.GetValue<double>("M_BD_1",Get("Mass_bottom")+Get("Mass_down"));
  m_parametermap[string("Mass_bu0")]           = 
    dataread.GetValue<double>("M_BU_0",Get("Mass_bottom")+Get("Mass_up"));
  m_parametermap[string("Mass_bu1")]           = 
    dataread.GetValue<double>("M_BU_1",Get("Mass_bottom")+Get("Mass_up"));
  m_parametermap[string("Mass_bs0")]           = 
    dataread.GetValue<double>("M_BS_0",Get("Mass_bottom")+Get("Mass_strange"));
  m_parametermap[string("Mass_bs1")]           = 
    dataread.GetValue<double>("M_BS_1",Get("Mass_bottom")+Get("Mass_strange"));
  m_parametermap[string("Mass_bc0")]           = 
    dataread.GetValue<double>("M_BC_0",Get("Mass_bottom")+Get("Mass_charm"));
  m_parametermap[string("Mass_bc1")]           = 
    dataread.GetValue<double>("M_BC_1",Get("Mass_bottom")+Get("Mass_charm"));
  m_parametermap[string("Mass_bb1")]           = 
    dataread.GetValue<double>("M_BB_1",2.*Get("Mass_bottom"));
  m_parametermap[string("Mixing_Angle_0-")]    = 
    dataread.GetValue<double>("Mixing_0-",-0.3491);
  m_parametermap[string("Mixing_Angle_1+")]    = 
    dataread.GetValue<double>("Mixing_1+",0.6155);
}

double Hadronisation_Parameters::Get(string keyword) 
{
  m_piter = m_parametermap.find(keyword);
  if (m_piter!=m_parametermap.end()) return m_piter->second;
  msg.Error()<<"Error in Hadronisation_Parameters::Get("<<keyword<<") in "<<m_parametermap.size()<<endl
	     <<"   Keyword not found. Return 0 and hope for the best."<<endl;
  return 0.;
}

bool Hadronisation_Parameters::AdjustMomenta(const int n,ATOOLS::Vec4D * moms,const double * masses)
{
  if (n==1) return false;
  if (n!=2) {
    bool  prepare=false,boost=false,success=true;
    Poincare rest;
    Vec4D cms = Vec4D(0.,0.,0.,0.);
    for (int i=0;i<n;i++) {
      cms += moms[i];
      if (dabs(moms[i].Abs2())>1.e-6) prepare = true;
    } 
    if (Vec3D(cms).Abs()>1.e-3) { 
      boost = true;
      rest  = Poincare(cms);
      for (int i=0;i<n;i++) rest.Boost(moms[i]);
    }
    if (prepare) success = success && p_stretcher->ZeroThem(n,moms);
    success = success && p_stretcher->MassThem(n,moms,masses);
    if (boost) {
      for (int i=0;i<n;i++) rest.BoostBack(moms[i]);
    }
    return success;
  } 
  else return p_stretcher->MassThem(n,moms,masses);
}
