#ifndef QCD_Splitting_Kernels_H
#define QCD_Splitting_Kernels_H

#include "Splitting_Kernel.H"

namespace PDF {

  class Q_QG: public Splitting_Kernel {
  public:
    
    // constructor
    Q_QG(const ATOOLS::Flavour &q);

    // member functions
    double operator()(const double &z) const;
    double Integral(const double &zmin,const double &zmax) const;

  };// end of class Q_QG

  class Q_GQ: public Splitting_Kernel {
  private:

    Q_QG m_swaped;
    
  public:
    
    // constructor
    Q_GQ(const ATOOLS::Flavour &q);

    // member functions
    double operator()(const double &z) const;
    double Integral(const double &zmin,const double &zmax) const;

  };// end of class Q_GQ

  class G_QQ: public Splitting_Kernel {
  public:
    
    // constructor
    G_QQ(const ATOOLS::Flavour &q);

    // member functions
    double operator()(const double &z) const;
    double Integral(const double &zmin,const double &zmax) const;

  };// end of class G_QQ

  class G_GG: public Splitting_Kernel {
  public:
    
    // constructor
    G_GG();

    // member functions
    double operator()(const double &z) const;
    double Value(const double &z) const;
    double Integral(const double &zmin,const double &zmax) const;

  };// end of class G_GG

}// end of namespace PDF

#endif
