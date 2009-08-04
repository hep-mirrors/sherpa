#ifndef METOOLS_Main_Partial_Amplitude_Base_H
#define METOOLS_Main_Partial_Amplitude_Base_H

#include "METOOLS/Main/Spin_Structure.H"

namespace METOOLS {

  std::string GetName(ATOOLS::Flavour* flavs, int n);

  class  Partial_Amplitude_Base : public Spin_Structure<Complex> {
  protected:
    ATOOLS::Flavour         * p_flavs;
    int             * p_i; // external indices for flavs and moms
    bool            * p_out; // incoming/outgoing

    void AssertSpins(int, ...);
    void AssertIn(int in);
  public:
    Partial_Amplitude_Base(ATOOLS::Flavour* flavs,size_t size,int* i,bool* out);
    virtual ~Partial_Amplitude_Base();

    virtual void    operator()(const ATOOLS::Vec4D * moms,bool anti=false)=0;
    static Partial_Amplitude_Base* Select(ATOOLS::Flavour* flavs, int n);
  };

  class Isotropic : public Partial_Amplitude_Base {
  public:
    Isotropic(ATOOLS::Flavour *fl,int n,int *i,bool *out);
    ~Isotropic() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

}


#endif
