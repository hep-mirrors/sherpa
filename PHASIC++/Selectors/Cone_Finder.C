#ifndef PHASIC_Selectors_Cone_Finder_H
#define PHASIC_Selectors_Cone_Finder_H

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
    Cone_Finder(int,ATOOLS::Flavour*,double);
    
    void   Init(const ATOOLS::Vec4D *);
    bool   Trigger(const ATOOLS::Vec4D_Vector &);
    void   BuildCuts(Cut_Data *);
    void   UpdateCuts(double,double,Cut_Data *);
    int    IsConditional() { return 1; }
    double ActualValue() const { return m_value; }
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

void Cone_Finder::Init(const Vec4D * p)
{
  //nothing has to be done for hadronic collisions
  //think about lepton-lepton collisions
}

Cone_Finder::Cone_Finder(int _n,Flavour * _fl,double _rcone) : 
  Selector_Base("Condefinder"), m_rcone(_rcone), m_value(0.)
{
  m_name = std::string("Conefinder");
 
  m_n    = _n;
  m_nin  = 2; 
  m_nout = m_n-2; 
  m_smin = 0.;
  m_fl   = _fl;
  
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
	if (r2jk<sqr(m_rcone)) return sqrt(r2jk);
	r2min = r2jk;
      }
    }
  }   
  return sqrt(r2min);
} 

bool Cone_Finder::Trigger(const Vec4D_Vector &p)
{
  // create copy
  Vec4D * moms = new Vec4D[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) moms[i]=p[i];

  Init(moms);

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
      if (m_fl[i].Mass()<3. && m_fl[j].Mass()<3. &&
	  !(m_fl[i].IsLepton() && m_fl[i].IntCharge()==0) &&
	  !(m_fl[j].IsLepton() && m_fl[j].IntCharge()==0))   {
	double mp2=Max(sqr(cuts->etmin[i])-sqr(m_fl[i].Mass()),
		       (sqr(cuts->energymin[i])-sqr(m_fl[i].Mass()))*(1-sqr(cuts->cosmax[0][i])))*
	           Max(sqr(cuts->etmin[j])-sqr(m_fl[j].Mass()),
		       (sqr(cuts->energymin[j])-sqr(m_fl[j].Mass()))*(1-sqr(cuts->cosmax[0][j])));
	cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],2.0*sqrt(mp2)*rp2);
      }
    }
  }
}

void   Cone_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) {}

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

namespace PHASIC{

DECLARE_ND_GETTER(Cone_Finder_Getter,"ConeFinder",Selector_Base,Selector_Key,true);

Selector_Base *Cone_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<1) THROW(critical_error,"Invalid syntax");
  Cone_Finder *jf(new Cone_Finder(key.p_proc->NIn()+key.p_proc->NOut(),
				  (Flavour*)&key.p_proc->Process()->
				  Flavours().front(),ToType<double>
				  (key.p_read->Interpreter()->Interprete(key[0][0]))));
  jf->SetProcess(key.p_proc);
  return jf;
}

void Cone_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"cone jet finder"; 
}

}

