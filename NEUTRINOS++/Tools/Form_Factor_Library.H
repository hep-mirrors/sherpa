#ifndef NEUTRINOS_Tools_Form_Factor_Library_H
#define NEUTRINOS_Tools_Form_Factor_Library_H

namespace NEUTRINOS {
  class Form_Factor_Base {
  protected:
    double m_norm, m_arg;
  public:
    Form_Factor_Base(const double & norm=1.,const double & arg=0.) :
      m_norm(norm), m_arg(arg) {}
    ~Form_Factor_Base() {}
    virtual double Calc(const double & q2) = 0;
  };

  class Dipole_Form_Factor : public Form_Factor_Base {
  public:
    Dipole_Form_Factor(const double & norm=1.,const double & arg=0.) :
      Form_Factor_Base(norm, arg) {}
    ~Dipole_Form_Factor() {}
    double Calc(const double & q2);
  };

  class Exponential_Form_Factor : public Form_Factor_Base {
  public:
    Exponential_Form_Factor(const double & norm=1.,const double & arg=0.) :
      Form_Factor_Base(norm, arg) {}
    ~Exponential_Form_Factor() {}
    double Calc(const double & q2);
  };

  class Gaussian_Form_Factor : public Form_Factor_Base {
  public:
    Gaussian_Form_Factor(const double & norm=1.,const double & arg=0.) :
      Form_Factor_Base(norm, arg) {}
    ~Gaussian_Form_Factor() {}
    double Calc(const double & q2);
  };
}
#endif
