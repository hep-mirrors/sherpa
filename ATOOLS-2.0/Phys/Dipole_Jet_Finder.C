#include "Dipole_Jet_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "Blob.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ATOOLS;
using namespace std;



Dipole_Jet_Finder::Dipole_Jet_Finder(const int n,Flavour * fl,
				     const double ycut,
				     const dipjet_type::code type,
				     const dipjet_mode::code mode,
				     const bool colsorted) :
  m_stot(sqr(rpa.gen.Ecms())), m_kt2cut(m_ycut*m_stot), 
  m_diptype(type), m_mode(mode),
  m_coloursorted(colsorted)
{
  m_ycut = ycut;
  m_name = std::string("Dipole_Jetfinder");
  m_fl   = fl;
  m_n    = n;
  if (m_mode==dipjet_mode::dec) { m_nin = 1; m_nout = m_n-1; }
                           else { m_nin = 2; m_nout = m_n-2; }
  
  m_moms.resize(n);
  m_flavs.resize(n);
  for (int i=0;i<n;i++) m_flavs.push_back(fl[i]);

  m_smin    = m_ycut * m_stot;
  m_smax    = m_stot;
  m_sel_log = new Selector_Log(m_name);
}

Dipole_Jet_Finder::Dipole_Jet_Finder(const double ycut,
				     const dipjet_type::code type,
				     const dipjet_mode::code mode,
				     const bool colsorted) :
  m_stot(sqr(rpa.gen.Ecms())), m_kt2cut(m_ycut*m_stot), 
  m_diptype(type), m_mode(mode),
  m_coloursorted(colsorted)
{
  m_ycut=ycut;
}

Dipole_Jet_Finder::~Dipole_Jet_Finder() {}

bool Dipole_Jet_Finder::Trigger(const Vec4D * p) 
{
  PROFILE_HERE;
  // create copy
  for (int i(0);i<m_nin+m_nout;++i) {
    msg.Out()<<METHOD<<" "<<i<<":"<<p[i]<<"/"<<m_moms.size()<<endl;
    m_moms[i]=p[i];
  }
  Init(&m_moms.front());

  FF_Winner();
  IF_Winner();
  FI_Winner();
  II_Winner();

  m_kt2min = m_kt2minFF;
  m_rec    = m_recFF;
  if (m_kt2minIF<m_kt2min) {
    m_kt2min = m_kt2minIF;
    m_rec    = m_recIF;
  }
  if (m_kt2minFI<m_kt2min) {
    m_kt2min = m_kt2minFI;
    m_rec    = m_recFI;
  }
  if (m_kt2minII<m_kt2min) {
    m_kt2min = m_kt2minII;
    m_rec    = m_recII;
  }


  bool trigger(m_kt2min>m_kt2cut);
  m_actual_y = m_kt2min/m_stot;

  msg.Out()<<METHOD<<":"<<m_kt2min<<"("<<m_kt2minFF<<") --> "
  	   <<m_actual_y<<" --> "<<trigger<<endl
	   <<"-------------------------------------------------"<<endl;
  return (1-m_sel_log->Hit(1-trigger));
}

void Dipole_Jet_Finder::Init(const Vec4D * p) {
  // No boost needed : Everything in Lorentz invariant fashion
  msg.Out()<<METHOD<<" : "<<p[0]<<"/"<<p[1]<<endl;
  if (m_mode==dipjet_mode::dec) 
    m_kt2minII = m_kt2minIF = m_kt2minFI = m_kt2minFF = m_shat = (p[0]).Abs2();
  else
    m_kt2minII = m_kt2minIF = m_kt2minFI = m_kt2minFF = m_shat = (p[0]+p[1]).Abs2();
}

//-----------------------------------------------------------------------------
// Building the cuts
//-----------------------------------------------------------------------------

void Dipole_Jet_Finder::BuildCuts(Cut_Data * cuts)
{
  int hadron;
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = Max(cuts->energymin[i],m_fl[i].SelMass());
    if (m_fl[i].Strong()) {                
      /* 
	 minimal energies : 
	 either   :  E^2 > kt^2 > y_cut s      
	             (hadron-hadron collisions)
	 or       :  4 min{E_i^2,E_j^2} > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	             ycut s' > ycut s_min   
	             (lepton-lepton collisions)
      */
      switch (m_mode) {
      case (dipjet_mode::hh): 
	cuts->energymin[i] = Max(sqrt(m_kt2cut),cuts->energymin[i]);
	cuts->cosmax[0][i] = cuts->cosmax[1][i] = cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
	  Min(cuts->cosmax[0][i],sqrt(1.-4.*m_ycut));
	cuts->etmin[i] = Max(sqrt(m_kt2cut),cuts->etmin[i]);
	break;
      case (dipjet_mode::DIS): 
	cuts->energymin[i] = Max(sqrt(m_kt2cut),cuts->energymin[i]);
	hadron=m_fl[0].Strong()?0:1;
	cuts->cosmax[hadron][i] = cuts->cosmax[i][hadron] = 
	  Min(cuts->cosmax[hadron][i],sqrt(1.-4.*m_ycut));
	cuts->cosmin[hadron][i] = cuts->cosmin[i][hadron] = 
	  Max(cuts->cosmin[hadron][i],-sqrt(1.-4.*m_ycut));
	cuts->etmin[i] = Max(sqrt(m_kt2cut),cuts->etmin[i]);
	break;
      case (dipjet_mode::ee): 
      default:
	cuts->energymin[i] = Max(sqrt(m_kt2cut),cuts->energymin[i]);
	msg.Out()<<METHOD<<" : E_min("<<i<<") = "<<cuts->energymin[i]<<endl;
	break;
      }      
      for (int j=i+1; j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong()) {
	  cuts->scut[j][i] = cuts->scut[i][j] = Max(cuts->scut[i][j],m_kt2cut);	
	  msg.Out()<<METHOD<<" : s_min("<<i<<","<<j<<") = "<<cuts->scut[i][j]<<endl;
	}
      }
    }
  }
}

void Dipole_Jet_Finder::UpdateCuts(double,double,Cut_Data *)
{
}



//-----------------------------------------------------------------------------
// Getting the smallest kt and its dipole from given momenta
//-----------------------------------------------------------------------------

void Dipole_Jet_Finder::II_Winner()
{
  m_kt2minII = m_shat;
  m_recII.i=m_recII.j=m_recII.k=-1;
  if (m_mode==dipjet_mode::hh) {
    double pt2;
    for (int i=m_nin;i<m_nin+m_nout;i++) {
      if (!m_flavs[i].Strong()) continue;
      pt2 = KT_min_II(1,2,i);
      if (pt2<m_kt2minII) m_kt2minII=pt2;
    }
  }
}

void Dipole_Jet_Finder::IF_Winner()
{
  m_kt2minIF = m_shat;
  m_recIF.i=m_recIF.j=m_recIF.k=-1;
  if (m_mode!=dipjet_mode::ee &&
      !(m_mode==dipjet_mode::dec && !m_flavs[0].Strong())) {
    double pt2;
    int perm, start = (m_mode==dipjet_mode::dec)?1:2;
    for (int i=0;i<m_nin;i++) {
      if (!m_flavs[i].Strong()) continue;
      for (int j=start;j<m_nin+m_nout-1;j++) {
	if (!m_flavs[j].Strong()) continue;
	for (int k=j+1;k<m_nin+m_nout;k++) {
	  if (!m_flavs[k].Strong()) continue;
	  pt2 = KT_min_IF(i,j,k,perm);
	  if (pt2<m_kt2minIF) {
	    m_kt2minIF=pt2;
	    if (perm==1) { m_recIF.i=i;m_recIF.j=j;m_recIF.k=k; }
	    if (perm==2) { m_recIF.i=i;m_recIF.j=k;m_recIF.k=j; }
	  }
	}
      }
    }
  }
}

void Dipole_Jet_Finder::FI_Winner() {}

void Dipole_Jet_Finder::FF_Winner()
{
  m_kt2minFF = m_shat;
  m_recFF.i=m_recFF.j=m_recFF.k=-1;
  double pt2;
  int perm, start = (m_mode==dipjet_mode::dec)?1:2;
  for (int i=start;i<m_nin+m_nout-2;i++) {
    //if (!m_flavs[i].Strong()) continue;
    for (int j=i+1;j<m_nin+m_nout-1;j++) {
      //if (!m_flavs[j].Strong()) continue;
      for (int k=j+1;k<m_nin+m_nout;k++) {
	//if (!m_flavs[k].Strong()) continue;
	pt2 = KT_min_FF(i,j,k,perm);
	if (pt2<m_kt2minFF) {
	  m_kt2minFF=pt2;
	  if (perm==1) { m_recIF.i=i;m_recIF.j=j;m_recIF.k=k; }
	  if (perm==2) { m_recIF.i=j;m_recIF.j=k;m_recIF.k=i; }
	  if (perm==3) { m_recIF.i=k;m_recIF.j=i;m_recIF.k=j; }
	}
      }
    }
  }
}

double Dipole_Jet_Finder::KT_min_1I()
{
  double min1I(m_shat); 
  if (m_mode==dipjet_mode::DIS || m_mode==dipjet_mode::hh) {
    double pt2;
    for (int i=2;i<m_nin+m_nout;i++) {
      pt2 = m_moms[i].PPerp2();
      if (pt2<min1I) min1I = pt2;
    }
  }
  return min1I;
}

double Dipole_Jet_Finder::KT_min_II(const int i,const int j,const int k)
{
  return -1.0;
}

double Dipole_Jet_Finder::KT_min_IF(const int i,const int j,const int k,
				    int & perm)
{
  return -1.0;
}

double Dipole_Jet_Finder::KT_min_FI(const int i,const int j,const int k,
				    int & perm)
{
  return -1.0;
}

double Dipole_Jet_Finder::KT_min_FF(const int i,const int j,const int k,
				    int &perm)
{
  double sprime((m_moms[i]+m_moms[j]+m_moms[k]).Abs2()), pt2min(sprime);
  perm = -1;

  double 
    sij((m_moms[i]+m_moms[j]).Abs2()),
    sjk((m_moms[j]+m_moms[k]).Abs2()),
    ski((m_moms[k]+m_moms[i]).Abs2());

  // ij emitter
  double pt2(sjk*ski/sprime);
  if (pt2<pt2min) {
    pt2min = pt2;
    perm   = 1;
  }
  // jk emitter
  pt2 = ski*sij/sprime;
  if (pt2<pt2min) {
    pt2min = pt2;
    perm   = 2;
  }
  // ki emitter
  pt2 = sij*sjk/sprime;
  if (pt2<pt2min) {
    pt2min = pt2;
    perm   = 3;
  }
  msg.Out()<<METHOD<<" ("
	   <<" 0"<<m_moms[0]<<endl
	   <<"                               1"<<m_moms[1]<<endl
	   <<"                               "<<i<<m_moms[i]<<endl
	   <<"                               "<<j<<m_moms[j]<<endl
	   <<"                               "<<k<<m_moms[k]<<") : "
	   <<sij<<" / "<<sjk<<" / "<<ski<<endl;
  return pt2min;
}

//-----------------------------------------------------------------------------
// Analyse process name and prepare for cuts
//-----------------------------------------------------------------------------




