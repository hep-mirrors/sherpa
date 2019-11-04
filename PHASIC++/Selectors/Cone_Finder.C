#ifndef PHASIC_Selectors_Cone_Finder_h
#define PHASIC_Selectors_Cone_Finder_h

#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Poincare.H"

namespace PHASIC {
  class Cone_Finder : public Selector_Base {
    double m_rcone;
    double m_value;

    double Rmin(ATOOLS::Vec4D *);

    double DEta12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &);
    double DPhi12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &);

  public:
    Cone_Finder(Process_Base *const proc, double rcone);

    bool   Trigger(ATOOLS::Selector_List &);
    void   BuildCuts(Cut_Data *);
  };
}

#endif

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace PHASIC;
using namespace ATOOLS;

Cone_Finder::Cone_Finder(Process_Base *const proc, double rcone) :
  Selector_Base("Conefinder",proc),
  m_rcone(rcone),
  m_value(0.0)
{
  m_smin       = 0.0;
}

double Cone_Finder::Rmin(Vec4D * p)
{
  double r2min = 100000.;
  double r2jk; 
    
  for (int j=m_nin;j<m_n;j++) {
    for (int k=j+1;k<m_n;k++) {
      r2jk = sqr(DEta12(p[j],p[k])) + sqr(DPhi12(p[j],p[k]));
      if (r2jk<r2min && 
	  p_fl[j].Mass()<3. && p_fl[k].Mass()<3. &&
	  !(p_fl[j].IsLepton() && p_fl[j].IntCharge()==0) &&
	  !(p_fl[k].IsLepton() && p_fl[k].IntCharge()==0))   {
	if (r2jk<sqr(m_rcone)) return sqrt(r2jk);
	r2min = r2jk;
      }
    }
  }   
  return sqrt(r2min);
} 

bool Cone_Finder::Trigger(Selector_List &sl)
{
  Vec4D * moms = new Vec4D[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) moms[i]=sl[i].Momentum();

  bool   trigger = 1;
  double rmin   = Rmin(moms); 

  if (rmin<m_rcone) trigger = 0;

  delete [] moms;
    
  m_value = rmin;
  
  return (1-m_sel_log->Hit(1-trigger));
}

void Cone_Finder::BuildCuts(Cut_Data * cuts) 
{
  double rp2=1.0-cos(m_rcone);
  for (int i=m_nin;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {     
      if (p_fl[i].Mass()<3. && p_fl[j].Mass()<3. &&
	  !(p_fl[i].IsLepton() && p_fl[i].IntCharge()==0) &&
	  !(p_fl[j].IsLepton() && p_fl[j].IntCharge()==0))   {
	double mp2=Max(sqr(cuts->etmin[i])-sqr(p_fl[i].Mass()),
		       (sqr(cuts->energymin[i])-sqr(p_fl[i].Mass()))*(1-sqr(cuts->cosmax[0][i])))*
	           Max(sqr(cuts->etmin[j])-sqr(p_fl[j].Mass()),
		       (sqr(cuts->energymin[j])-sqr(p_fl[j].Mass()))*(1-sqr(cuts->cosmax[0][j])));
	cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],2.0*sqrt(mp2)*rp2);
      }
    }
  }
}


/*----------------------------------------------------------------------------

  Utilities

  ----------------------------------------------------------------------------*/
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


DECLARE_GETTER(Cone_Finder,"ConeFinder",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Cone_Finder>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings["ConeFinder"];
  auto R = s["R"].SetDefault(0.4).Get<double>();
  Cone_Finder *jf(new Cone_Finder(key.p_proc, R));
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Cone_Finder>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"ConeFinder:\n"
     <<width<<"  R: cone size";
}
