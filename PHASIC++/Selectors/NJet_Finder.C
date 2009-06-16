#ifndef PHASIC_Selectors_NJet_Finder_h
#define PHASIC_Selectors_NJet_Finder_h

#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Poincare.H"

namespace PHASIC {
  class NJet_Finder : public Selector_Base {
    double m_ycut,m_delta_r,m_r2min;
    int    m_type;
    bool   m_pt_def;
    double m_ene, m_s, m_sprime; //, m_smin, m_smax;
    int    m_nstrong;

    double ** p_ktij;
    int    *  p_imap;
    double *  p_kis;
    std::vector<double> m_jetpt;
    std::vector<double> m_kts;


    double DEta12(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double CosDPhi12(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double DPhi12(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double R2(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double Y12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &) const;
    double DCos12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &) const;

    void   AddToKtlist(double );
    void   AddToJetlist(double );

    double Ktmin(ATOOLS::Vec4D * p, int n);
    void   Ymin(ATOOLS::Vec4D * p, int n);
  public:
    NJet_Finder(int nin, int nout,ATOOLS::Flavour * fl, double ycut, double dr, int nn, int type);

    ~NJet_Finder();


    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    bool   Trigger(const ATOOLS::Vec4D_Vector &);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,const ATOOLS::Flavour_Vector &,int);

    void   BuildCuts(Cut_Data *);
    void   UpdateCuts(double,double,Cut_Data *);
    void   PrintPT();
  };
}

#endif

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"


using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

NJet_Finder::NJet_Finder(int nin, int nout,Flavour* fl, double ycut, double dr, int nn,int type) : 
  Selector_Base("NJetfinder"), m_ycut(ycut), m_delta_r(dr), m_type(type)
{
  m_fl         = fl;
  m_ene        = rpa.gen.Ecms()/2.;
  m_sprime     = m_s = sqr(2.*m_ene); 
  m_smin       = m_ycut * m_s;
  m_smax       = m_s;

  m_r2min      = sqr(m_delta_r);

  m_nin        = nin;
  m_nout       = nout;
  m_n          = nn;


  p_kis  = new double[m_nout];
  p_imap = new int[m_nout];
  p_ktij = new double*[m_nout];
  for (int i=0;i<m_nout;++i) p_ktij[i] = new double[m_nout];
  for (int i=0;i<m_nout;++i) p_imap[i] = i;
  
  m_nstrong = 0;
  for (int i=m_nin;i<m_nout+m_nin;i++) {
    if (fl[i].Strong()) m_nstrong++;
  }

  m_sel_log    = new Selector_Log(m_name);
}


NJet_Finder::~NJet_Finder() {
      for (int i=0;i<m_nout;++i) delete [] p_ktij[i];
      delete [] p_ktij;
      delete [] p_imap;
      delete [] p_kis;
}


void NJet_Finder::AddToKtlist(double kt2) {
  m_kts.push_back(kt2);
}

void NJet_Finder::AddToJetlist(double kt2) {
  m_jetpt.push_back(kt2);
}

void NJet_Finder::PrintPT() {
for(size_t i=0;i<m_kts.size();i++) std::cout<<"kt "<<i<<": "<<m_kts[i]<<std::endl;
}

bool NJet_Finder::NoJetTrigger(const Vec4D_Vector &p)
{
  double s=(p[0]+p[1]).Abs2();
  return (s>m_smin*4.);
}

bool NJet_Finder::Trigger(const Vec4D_Vector &p)
{
  if (m_n<1) return true;

  // create copy
  m_jetpt.clear();
  m_kts.clear();
  int n=0;
  Vec4D * moms = new Vec4D[m_nout];
  for (int i=m_nin;i<m_nout+m_nin;i++) {
//      if (m_fl[i].Strong()&&!(m_fl[i].IsMassive())) {
       if (m_fl[i].Strong()) {
      moms[n]=p[i];
      n++;
    }
  }
  // cluster
  for (int i=0;i<n;++i) p_imap[i] = i;
  if (m_type==1) Ymin(moms,n);
  else Ktmin(moms,n);

  delete [] moms;
  bool trigger(true);
  if (n<m_n) return 0;
  if (m_type>1) {
    if (m_jetpt.size()<m_n) trigger=false;
    else if (m_jetpt[m_jetpt.size()-m_n]/m_s<m_ycut) trigger=false;
  } 
  else {
    if (m_kts[n-m_n]/m_s<m_ycut) trigger=false;
    if (m_kts[0]<=0.) trigger=false;
  }

  return (1-m_sel_log->Hit(1-trigger));
}

bool NJet_Finder::JetTrigger(const Vec4D_Vector &p,const Flavour_Vector &fl,int ns)
{
  if (m_n<1) return true;

  // create copy
  m_jetpt.clear();
  m_kts.clear();
  int n=0;
  Vec4D * moms = new Vec4D[ns];
  for (int i=m_nin;i<ns+m_nin;i++) {
//      if (m_fl[i].Strong()&&!(m_fl[i].IsMassive())) {
    if (fl[i].Strong()) {
      moms[n]=p[i];
      n++;
    }
  }

  // cluster
  for (int i=0;i<n;++i) p_imap[i] = i;

  if (m_type==1) Ymin(moms,n);
  else Ktmin(moms,n);

  delete [] moms;

  bool trigger(true);
  if (n<m_n) return 0;
  if (1&&m_type>1) {
    if (m_jetpt.size()<m_n) trigger=false;
    else if (m_jetpt[m_jetpt.size()-m_n]/m_s<m_ycut) trigger=false;
  } 
  else {
    if (m_kts[n-m_n]/m_s<m_ycut) trigger=false;
    if (m_kts[0]<=0.) trigger=false;
  }
  
  return (1-m_sel_log->Hit(1-trigger));
}

double NJet_Finder::Ktmin(Vec4D * p, int n)
{
  if (n==0) return 0.;
  if (n==1) {
    double dmin=(m_type>1)?p[0].PPerp2():sqr(p[0][0]);
    AddToJetlist(dmin);
    AddToKtlist(dmin);
    return dmin;
  }

  //cal first matrix
  int ii=0, jj=0;
  double dmin=(m_type>1)?p[0].PPerp2():sqr(p[0][0]);
  
  {
    
    for (int i=0;i<n;++i) {
      double di = 0.;
      if (m_type>1) di = p_ktij[i][i] = p[i].PPerp2();
      else di = p_ktij[i][i] = sqr(p[i][0]);
      if (di<dmin) { dmin=di; ii=i; jj=i;}
      for (int j=0;j<i;++j) {
	double dj  = p_ktij[j][j]; 
	double dij = p_ktij[i][j] = Min(di,dj)*R2(p[i],p[j]) /m_r2min;
	if (dij<dmin) {dmin=dij; ii=i; jj=j;}
      }
    }
  }
  
  // recalc matrix
  while (n>0) {
    if (ii!=jj) {
      // combine precluster
      p[p_imap[jj]]+=p[p_imap[ii]];
      AddToKtlist(dmin);
    }
    else {
      // add to jet list
      AddToJetlist(dmin);
      AddToKtlist(dmin);
    }

    --n;

    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];

    if (n==1) {
      int jjx=p_imap[jj];
      p_ktij[jjx][jjx] = (m_type>1)?p[jjx].PPerp2():sqr(p[jjx][0]);
      break;
    }
    // update map (remove precluster)
    {

    
    // update matrix (only what is necessary)
    int jjx=p_imap[jj];
    p_ktij[jjx][jjx] = (m_type>1)?p[jjx].PPerp2():sqr(p[jjx][0]);
    for (int j=0;j<jj;++j)   p_ktij[jjx][p_imap[j]] = 
			       Min(p_ktij[jjx][jjx],p_ktij[p_imap[j]][p_imap[j]])
			       *R2(p[jjx],p[p_imap[j]])/m_r2min;
    for (int i=jj+1;i<n;++i) p_ktij[p_imap[i]][jjx] = 
			       Min(p_ktij[jjx][jjx],p_ktij[p_imap[i]][p_imap[i]])
			       *R2(p[p_imap[i]],p[jjx])/m_r2min;
    }
    // redetermine rmin and dmin
    ii=jj=0;
    {

    dmin=p_ktij[p_imap[0]][p_imap[0]];
    for (int i=0;i<n;++i) {
      int ix=p_imap[i];
      double di = p_ktij[ix][ix];
      if (di<dmin) { dmin=di; ii=jj=i;}
      for (int j=0;j<i;++j) {
	int jx=p_imap[j];
	double dij = p_ktij[ix][jx];
	if (dij<dmin) { dmin=dij; ii=i; jj=j;}
      }
    }
    }
  }

  // add remaining preclusters to jetlist
  for (int i=0;i<n;++i) {
    AddToJetlist(p_ktij[p_imap[i]][p_imap[i]]);
    AddToKtlist(p_ktij[p_imap[i]][p_imap[i]]);
  }
  return dmin;
}

void NJet_Finder::Ymin(Vec4D * p, int n)
{
  if (n==0) return;
  if (n==1) {
    AddToKtlist(sqr(p[0][0]));
    return;
  }

  //cal first matrix
  int ii=0, jj=0;
  double ymin=1.;
  {
    
    for (int i=1;i<n;++i) {
      for (int j=0;j<i;++j) {
	double y = p_ktij[i][j] =Y12(p[i],p[j]);
	if (y<ymin) { ymin=y; ii=i; jj=j;}
      }
    }
  }

  // recalc matrix
  while (n>1) {

    // combine precluster
    p[p_imap[jj]]+=p[p_imap[ii]];
    AddToKtlist(ymin*m_sprime);
    --n;
    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];

    if (n==1) break;
    // update map (remove precluster)
    {

      // update matrix (only what is necessary)
      int jjx=p_imap[jj];
      for (int j=0;j<jj;++j)   p_ktij[jjx][p_imap[j]] = Y12(p[jjx],p[p_imap[j]]);
      for (int i=jj+1;i<n;++i) p_ktij[p_imap[i]][jjx] = Y12(p[p_imap[i]],p[jjx]);
    }
    // redetermine rmin and dmin
    ii=jj=0;
    ymin=1.;
    {

      for (int i=1;i<n;++i) {
	int ix=p_imap[i];
	for (int j=0;j<i;++j) {
	  int jx=p_imap[j];
	  double y = p_ktij[ix][jx];
	  if (y<ymin) { ymin=y; ii=i; jj=j;}
	}
      }
    }
  }
  AddToKtlist(ymin*m_sprime);
}


void NJet_Finder::BuildCuts(Cut_Data * cuts) 
{
  return;
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = m_fl[i].SelMass();
    if (m_fl[i].Strong()) {                
      /* 
	 minimal energies : 
	 either   :  E^2 > kt^2 > y_cut s      
	             (hadron-hadron collisions)
	 or       :  4 min{E_i^2,E_j^2} > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	             ycut s' > ycut s_min   
	             (lepton-lepton collisions)
      */
      cuts->energymin[i] = Max(sqrt(1. * m_ycut * m_s),cuts->energymin[i]);
      cuts->cosmax[0][i] = cuts->cosmax[1][i] = cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
	Min(cuts->cosmax[0][i],sqrt(1.-4.*m_ycut));
      
      for (int j=i+1; j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong()) {                
	  /* 
	     minimal scut :
	     either   :  s_ij = 2 E_i E_j (1-cos(ij)) > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	     m_ycut s' > ycut s_min   
	     (lepton-lepton collisions)
	     or       :  similarly .... have to think ...
               	         (hadron-hadron collisions)
 
	  */
	  cuts->scut[j][i] = cuts->scut[i][j] 
	    = Max(cuts->scut[i][j],sqr(m_delta_r)*m_ycut*m_s);
	}
      }
    }
  }
}

void   NJet_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{
  return;
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = Max(sqrt(m_ycut * sprime/4.),cuts->energymin[i]);
    for (int j=i+1; j<m_nin+m_nout; ++j) {
      if (m_fl[j].Strong()) {                
	/* 
	   minimal scut :
	   either   :  s_ij = 2 E_i E_j (1-cos(ij)) > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	               ycut s' > ycut s_min   
	               (lepton-lepton collisions)
	   or       :  similarly .... have to think ...
	               (hadron-hadron collisions)
	   
	*/
	cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],m_ycut*sprime);
      }
    }
  }
}


/*----------------------------------------------------------------------------------

  Utilities

  ----------------------------------------------------------------------------------*/
double NJet_Finder::DEta12(Vec4D & p1,Vec4D & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  double c1=p1[3]/Vec3D(p1).Abs();
  double c2=p2[3]/Vec3D(p2).Abs();
  return  0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
}

double NJet_Finder::CosDPhi12(Vec4D & p1,Vec4D & p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return (p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2);
}

double NJet_Finder::DPhi12(Vec4D & p1,Vec4D & p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
}


double NJet_Finder::R2(Vec4D &p1, Vec4D &p2)
{
  return sqr(p1.Y()-p2.Y()) + sqr(DPhi12(p1,p2));
//   return 2.*(cosh(DEta12(p1,p2)) - CosDPhi12(p1,p2));
  return (sqr(DEta12(p1,p2)) + sqr(DPhi12(p1,p2)));
}

double NJet_Finder::Y12(const Vec4D & p1, const Vec4D & p2) const
{
  return 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2))/m_sprime;
}

double NJet_Finder::DCos12(const Vec4D & p1,const Vec4D & p2) const
{
  double s  = p1[1]*p2[1]+p1[2]*p2[2]+p1[3]*p2[3];
  double b1 = p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3];
  double b2 = p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3];
  return s/sqrt(b1*b2);
  //  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

namespace PHASIC{

DECLARE_ND_GETTER(NJet_Finder_Getter,"NJetFinder",Selector_Base,Selector_Key,true);

Selector_Base *NJet_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<2) THROW(critical_error,"Invalid syntax");
  int type(0);
  if (key.p_proc->NIn()==2) {
    int instrong(0);
    for (size_t j=0;j<key.p_proc->NIn();j++) { 
      if (key.p_proc->Process()->
	  Flavours()[j].Strong()) instrong++; 
    }
    if (instrong==0) type = 1;
    if (instrong==1) type = 2;
    if (instrong==2) type = 4;
  }
  NJet_Finder *jf(new NJet_Finder(key.p_proc->NIn(),key.p_proc->NOut(),
				  (Flavour*)&key.p_proc->Process()->Flavours().front(),
				  ToType<double>(key.p_read->Interpreter()->Interprete(key[0][0])),
				  ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1])),
				  ToType<int>(key[0][2]),type));
  jf->SetProcess(key.p_proc);
  return jf;
}

void NJet_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"n jet finder"; 
}

}
