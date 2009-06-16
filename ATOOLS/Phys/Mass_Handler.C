#include "ATOOLS/Phys/Mass_Handler.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;
using namespace std;

Mass_Handler::Mass_Handler(Flavour flav)
{
  double peak = flav.HadMass();
  double width = flav.Width();
  switch(flav.Kfcode()) {
  default: p_mass = new Relativistic_Breit_Wigner(peak, width);
  }
}

Mass_Handler::~Mass_Handler()
{
  if(p_mass) { delete p_mass; p_mass=NULL; }
}

double Mass_Handler::GetMass(double min, double max)
{
  if(min<0.0 || max<0.0 || min>max) {
    msg_Error()<<METHOD<<" range not valid: min="<<min<<" max="<<max<<endl;
    return 0.0;
  }
  double m = p_mass->GetMass(min, max);
  if (!(m>0) && !(m<0) && m!=0) {
  msg_Error()<<METHOD<<"produced a nan:"<<endl
      <<"  peak="<<p_mass->Peak()<<" width="<<p_mass->Width()
      <<" min="<<min<<" max="<<max<<endl;
  }
  return m;
}

// ========================================================

Breit_Wigner::Breit_Wigner(double peak,double width)
{
  m_peak=peak; m_width=width;
}

double Breit_Wigner::GetMass(double min, double max)
{
  msg_Error()<<METHOD<<" not yet implemented."<<endl;
  abort();
}

// ========================================================

Relativistic_Breit_Wigner::Relativistic_Breit_Wigner(double peak,double width)
{
  m_peak=peak; m_width=width;
}

double Relativistic_Breit_Wigner::GetMass(double min, double max)
{
  if( m_peak<1.e-6 || m_width/m_peak < 1.e-8) return m_peak;
  double random = ran.Get();
  double peak2 = m_peak*m_peak;
  double mw    = m_peak*m_width;
  double s;
  if (min==0.0 && max==MYINF) s = peak2+mw*tan(M_PI*(random-0.5));
  else {
    double smin = sqr(min); double smax = sqr(max);
    double ymax=atan((smin-peak2)/mw);
    double ymin=atan((smax-peak2)/mw);
    s = peak2+mw*tan(ymin + random*(ymax-ymin));
  }
  return sqrt(s);
}
