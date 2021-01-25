#ifndef ATOOLS_Phys_Fragmentation_Base_H
#define ATOOLS_Phys_Fragmentation_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  class Fragmentation_Base {
  protected:
    std::string m_name;
  public:
    Fragmentation_Base();
    virtual ~Fragmentation_Base();

    virtual Return_Value::code Hadronize(ATOOLS::Blob_List *) = 0;
    
    const std::string & Name() const { return m_name; }
    void SetName(const std::string& name) { m_name=name; }
  };


  struct Fragmentation_Getter_Parameters {
    std::string shower;
    Fragmentation_Getter_Parameters(const std::string& _shower) :
      shower(_shower)
    {}
  };
  typedef Getter_Function<Fragmentation_Base,Fragmentation_Getter_Parameters> Fragmentation_Getter;
}

#define DEFINE_FRAGMENTATION_GETTER(CLASS,TAG)				        \
  DECLARE_GETTER(CLASS,TAG,Fragmentation_Base,Fragmentation_Getter_Parameters); \
  Fragmentation_Base *ATOOLS::Getter<Fragmentation_Base,Fragmentation_Getter_Parameters,CLASS>:: \
  operator()(const Fragmentation_Getter_Parameters &parameters) const \
  { Fragmentation_Base* ret = new CLASS(parameters.shower); ret->SetName(TAG); return ret;  } \
  void ATOOLS::Getter<Fragmentation_Base,Fragmentation_Getter_Parameters, CLASS>:: \
  PrintInfo(std::ostream &str,const size_t width) const \
  { str<<TAG<<" fragmentation"; }

#endif

