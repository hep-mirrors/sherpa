#include "Cone_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace ATOOLS;


void Cone_Finder::Init(const Vec4D * p)
{
  //nothing has to be done for hadronic collisions
  //think about lepton-lepton collisions
}

double * Cone_Finder::ActualValue() { 
  return m_value; 
}

Cone_Finder::Cone_Finder(int _n,Flavour * _fl,double _rcone) : 
  m_rcone(_rcone)
{
  m_name    = std::string("Conefinder");
 
  m_n       = _n;
  m_nin     = 2; 
  m_nout    = m_n-2; 

  m_smin = 0.;
  m_fl   = _fl;
  
  m_value   = new double[1];
  
  m_sel_log = new Selector_Log(m_name);
}

double Cone_Finder::Rmin(Vec4D * p)
{
  double r2min = 100000.;
  double r2jk; 
    
  for (int j=m_nin;j<m_n;j++) {
    for (int k=j+1;k<m_n;k++) {
      r2jk = sqr(DEta12(p[j],p[k])) + sqr(DPhi12(p[j],p[k]));
      if (r2jk<r2min && 
	  m_fl[j].Mass()<3. && m_fl[k].Mass()<3. &&
	  !(m_fl[j].IsLepton() && m_fl[j].IntCharge()==0) &&
	  !(m_fl[k].IsLepton() && m_fl[k].IntCharge()==0))   {
	r2min = r2jk;
      }
    }
  }    
  return sqrt(r2min);
} 

bool Cone_Finder::Trigger(const Vec4D * p)
{
  // create copy
  Vec4D * moms = new Vec4D[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) moms[i]=p[i];

  Init(moms);

  bool   trigger = 1;
  double rmin   = Rmin(moms); 
  
  if (rmin<m_rcone) trigger = 0;

  delete [] moms;
    
  m_value[0] = rmin;
  
  return (1-m_sel_log->Hit(1-trigger));
}

void Cone_Finder::BuildCuts(Cut_Data * cuts) 
{
}

void   Cone_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
}

double Cone_Finder::DEta12(const Vec4D & p1,const Vec4D & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  double c1=p1[3]/Vec3D(p1).Abs();
  double c2=p2[3]/Vec3D(p2).Abs();
  return  0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
}

double Cone_Finder::DPhi12(const Vec4D & p1,const Vec4D & p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
}

