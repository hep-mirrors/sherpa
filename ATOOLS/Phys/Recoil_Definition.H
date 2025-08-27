#ifndef ATOOLS__Phys__Recoil_Definition_H
#define ATOOLS__Phys__Recoil_Definition_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  struct RecoilDefinition_Key { };

  class Recoil_Definition {
  public:

    static void ShowSyntax(const int mode);

    virtual Vec4D Recoil(const Cluster_Amplitude *ampl,
			 size_t splits=0,size_t specs=0,int mode=-1) = 0;

    virtual std::vector<int> RecoilTags
    (const Cluster_Amplitude *ampl,
     size_t splits=0,size_t specs=0,int mode=-1) = 0;

    virtual int Mode(const Cluster_Amplitude *ampl,int split) const = 0;
    
    double m_prefac = 1;
    
  };// end of class Recoil_Definition

  typedef ATOOLS::Getter_Function
  <Recoil_Definition,RecoilDefinition_Key>
  RecoilDefinition_Getter;

}// end of namespace ATOOLS

#endif
