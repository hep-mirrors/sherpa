#ifndef PHASIC__Main__Enhance_Observable_Base_H
#define PHASIC__Main__Enhance_Observable_Base_H

#include "ATOOLS/Org/Getter_Function.H"

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"

namespace PHASIC {

  class Process_Base;

  struct Enhance_Arguments {
    Process_Base *p_proc;
    std::string m_enhance;
    Enhance_Arguments(Process_Base *const proc=NULL,
		      const std::string &enhance=""):
      p_proc(proc), m_enhance(enhance) {}
  };// end of class Enhance_Arguments

  class Enhance_Observable_Base {
  public:

    typedef ATOOLS::Getter_Function
    <Enhance_Observable_Base,Enhance_Arguments> Getter_Function;

  protected:

    Process_Base *const p_proc;

  public:

    Enhance_Observable_Base(const Enhance_Arguments &args);

    virtual ~Enhance_Observable_Base() {};
    virtual double operator()(const ATOOLS::Vec4D *p,
			      const ATOOLS::Flavour *fl,
			      const size_t n) = 0;

  };// end of class Enhance_Observable_Base

}// end of namespace PHASIC

#endif
