#ifndef Ordering_H
#define Ordering_H

#include "Particle.H"
#include "Getter_Function.H"

namespace ATOOLS {

  class Order_Base {
  public :

    virtual ~Order_Base();

    virtual bool operator()(const Vec4D &a,const Vec4D &b) const;
    virtual bool operator()(const Particle &a,const Particle &b) const;
    virtual bool operator()(Particle * const &a,Particle * const &b) const;

    void ShowOrders(const int mode=1);

  };

  typedef Getter_Function<Order_Base,std::string> Order_Getter;

}
#endif
