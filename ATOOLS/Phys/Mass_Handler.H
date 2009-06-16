#ifndef ATOOLS_Phys_Mass_Handler_H
#define ATOOLS_Phys_Mass_Handler_H

#include "ATOOLS/Org/My_Limits.H"

#define MYINF std::numeric_limits<double>::max()

namespace ATOOLS {

  class Flavour;
  class Mass_Base;

  class Mass_Handler {
    Mass_Base * p_mass;
  public:
    Mass_Handler(Flavour flav);
    ~Mass_Handler();
    double GetMass(double min, double max);
  };

  class Mass_Base {
  protected:
    double m_peak, m_width;
  public:
    Mass_Base() {}
    virtual ~Mass_Base() {}
    virtual double GetMass(double min, double max)=0;
    double Peak() { return m_peak; }
    double Width() { return m_width; }
  };

  class Breit_Wigner: public Mass_Base {
  public:
    Breit_Wigner(double peak, double width);
    double GetMass(double min=0.0, double max=MYINF);
  };

  class Relativistic_Breit_Wigner: public Mass_Base  {
  public:
    Relativistic_Breit_Wigner(double peak, double width);
    double GetMass(double min=0.0, double max=MYINF);
  };
}

#endif
