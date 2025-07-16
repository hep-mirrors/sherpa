#ifndef ATOOLS_Phys_Ordering_H
#define ATOOLS_Phys_Ordering_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  typedef bool (*Order_Vec_Func)(const Vec4D &a,const Vec4D &b);
  typedef bool (*Order_Part_Func)(const Particle &a,const Particle &b);
  typedef bool (*Order_PPart_Func)(Particle *const &a,Particle *const &b);


  class Order_Base {
  protected:

    Order_Vec_Func   p_vec;
    Order_Part_Func  p_part;
    Order_PPart_Func p_ppart;

    Order_Base(Order_Vec_Func vec,Order_Part_Func part,Order_PPart_Func ppart):
      p_vec(vec), p_part(part), p_ppart(ppart) {}

  private:

    Order_Base(): p_vec(NULL), p_part(NULL), p_ppart(NULL) {}

  public:

    virtual ~Order_Base();

    inline bool operator()(const Vec4D &a,const Vec4D &b) const
    { return (*p_vec)(a,b); }
    inline bool operator()(const Particle &a,const Particle &b) const
    { return (*p_part)(a,b); }
    inline bool operator()(Particle *const &a,Particle *const &b) const
    { return (*p_ppart)(a,b); }

    void ShowOrders(const int mode=1);

  };

  typedef Getter_Function<Order_Base,std::string> Order_Getter;

}
#endif
