#ifndef HD_ME_Base_H
#define HD_ME_Base_H

#include "Flavour.H"
#include "Vector.H"
#include <string>

namespace HADRONS {
  class HD_ME_Base {
  protected:
    std::string       m_metype;
    int               m_nout;
    ATOOLS::Flavour * p_flavs;
    double          * p_masses,* p_masses2;
  public:
    HD_ME_Base(int,ATOOLS::Flavour *);
    virtual ~HD_ME_Base();
    virtual double operator()(const ATOOLS::Vec4D *) { return 1.; }

    virtual void AddVector(double,double) {}
  };

  class Isotropic : public HD_ME_Base {
  public:
    Isotropic(int,ATOOLS::Flavour *,
	      std::string _met=std::string(""));
  };
}


#endif
