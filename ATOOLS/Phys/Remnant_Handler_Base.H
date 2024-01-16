#ifndef ATOOLS_Phys_Remnant_Handler_Base_H
#define ATOOLS_Phys_Remnant_Handler_Base_H

#include "ATOOLS/Phys/Remnant_Base.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"

namespace ATOOLS {    
  class Basic_Remnant_Handler_Base {
  protected:
    std::string             m_name;
  public:
    Basic_Remnant_Handler_Base(PDF::ISR_Handler * isr);
    virtual ~Basic_Remnant_Handler_Base();

    virtual void Reset();
    virtual ATOOLS::Return_Value::code MakeBeamBlobs(ATOOLS::Blob_List * blobs,
						     const double & B=0.) = 0;
    
    const std::string & Name()      const { return m_name; }
    void SetName(const std::string& name) { m_name=name; }
  };

}

#endif
