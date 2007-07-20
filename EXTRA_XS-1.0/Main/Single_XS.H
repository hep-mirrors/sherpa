#ifndef Single_XS_H
#define Single_XS_H

#include "XS_Base.H"

namespace ATOOLS { class Blob_Data_Base; }

namespace EXTRAXS {

  class Single_XS: public XS_Base {
  public:

    typedef Single_XS *(*Getter_Function)(const size_t nin,const size_t nout,
					  const ATOOLS::Flavour *flavours,
					  const size_t nqed,const size_t nqcd);

  public:

    // constructors
    Single_XS(const size_t nin,const size_t nout,
	      const ATOOLS::Flavour *flavours,XS_Model_Base *const model=NULL);
    Single_XS(const size_t nin,const size_t nout,
	      const ATOOLS::Flavour *flavours,
	      const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
	      BEAM::Beam_Spectra_Handler *const beamhandler,
	      PDF::ISR_Handler *const isrhandler,
	      ATOOLS::Selector_Data *const selectordata,
	      XS_Model_Base *const model=NULL);

    // member functions
    size_t Size() const;

    XS_Base *const operator[](const size_t i) const;

    double Differential2();

    void SetISR(PDF::ISR_Handler *const isrhandler);
    void SetTotal(const int mode=1);

    double Differential(const double s,const double t,const double u);
    double operator()(const double s,const double t,const double u);
    void   ResetMax(int);
    void   OptimizeResult();

    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create=false);
    void WriteOutXSecs(std::ofstream &outfile);

    bool Tests();

    ATOOLS::Blob_Data_Base *OneEvent();
    ATOOLS::Blob_Data_Base *WeightedEvent(const int mode=0);

    template <class Process_Type>
    static Single_XS *GetProcess(const size_t nin,const size_t nout,
				 const ATOOLS::Flavour *flavours,
				 const size_t nqed,const size_t nqcd);

    template <class Process_Type>
    static Single_XS *GetProcess(const size_t nin,const size_t nout,
				 const ATOOLS::Flavour *flavours,
				 XS_Model_Base *const model,
				 const size_t nqed,const size_t nqcd);

  };// end of class Single_XS

}// end of namespace EXTRAXS

#endif



