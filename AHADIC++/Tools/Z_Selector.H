#ifndef AHADIC_Tools_Z_Selector_H
#define AHADIC_Tools_Z_Selector_H

namespace AHADIC {
  class Splitter_Base;
  
  class Z_Selector {
  private:
    Splitter_Base * p_splitterbase;
  public:
    Z_Selector();
    ~Z_Selector();
    void Init(Splitter_Base * splitterbase=0);

    double operator()(const double & zmin=1.e-6,const double & zmax=0.999999);
  };
}

#endif
