#include "Cone_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;


void Cone_Finder::Init(const AMATOOLS::Vec4D * p)
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
  AORGTOOLS::msg.Debugging()<<"Initialize the Cone_Finder : "<<std::endl
			    <<"rcone = "<<m_rcone<<std::endl;

  m_name    = std::string("Conefinder");
 
  m_n       = _n;
  m_nin     = 2; 
  m_nout    = m_n-2; 

  m_fl   = _fl;
  
  m_value   = new double[1];
  
  m_sel_log = new Selector_Log(m_name);
}

double Cone_Finder::Rmin(AMATOOLS::Vec4D * p)
{
  double r2min = 100000.;
  double r2jk, deta, dphi;
    
  for (int j=m_nin;j<m_n;j++) {
    for (int k=j+1;k<m_n;k++) {
      deta = sqr(DEta12(p[j],p[k]));
      dphi = sqr(DPhi12(p[j],p[k]));
      r2jk = sqr(DEta12(p[j],p[k])) + sqr(DPhi12(p[j],p[k]));
      //msg.Out()<<"deta "<<deta<<" dphi "<<dphi<<endl;
      //msg.Out()<<"p1 "<<p[j]<<" p[2] "<<p[k]<<" -> r2jk = "<<r2jk<<endl;
      //if (r2jk<r2min) {
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

bool Cone_Finder::Trigger(const AMATOOLS::Vec4D * p)
{
  // create copy
  Vec4D * moms = new Vec4D[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) moms[i]=p[i];

  Init(moms);

  bool   trigger = 1;
  double rmin   = Rmin(moms); 
  
  //msg.Out()<<" Trigger (Rcone vs. r_min) "<<Rcone<<" vs. "<<r_min<<endl;
  
  if (rmin<m_rcone) trigger = 0;
  
  delete [] moms;
    
  m_value[0] = rmin;
  
  return (1-m_sel_log->Hit(1-trigger));
}

void Cone_Finder::BuildCuts(Cut_Data * cuts) 
{
  /*
  msg.Debugging()<<"In Cone_Finder::BuildCuts"<<std::endl;
  // Loop over final state particles.
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = Max(cuts->energymin[i],fl[i].Mass());
  }
  */
}

void   Cone_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  /*
  
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = Max(cuts->energymin[i],fl[i].Mass());
    
  }
  */
}

double Cone_Finder::DEta12(AMATOOLS::Vec4D & p1,AMATOOLS::Vec4D & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  // => eta1 - eta2 
  //  = -log(tan1/2)+log(tan2/2) = log(tan2/tan1) = log(pt2/pt1 * pl1/pl2)
  return log(sqrt( (sqr(p2[1])+sqr(p2[2]))/(sqr(p1[1])+sqr(p1[2])) ) * dabs(p1[3]/p2[3]));
}

double Cone_Finder::DPhi12(AMATOOLS::Vec4D & p1,AMATOOLS::Vec4D & p2)
{
  // phi1-phi2 = acos(cos(phi1) cos(phi2) + sin(phi1) sin(phi2))
  //                = (p1_x p2_x + p1_y p2_y)/(p1_T * p2_T)
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  //msg.Out()<<"cos( phi1 - phi2 )="<<(p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2);
  return acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
}


