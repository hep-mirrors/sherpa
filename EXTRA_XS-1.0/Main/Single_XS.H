#ifndef Single_XS_H
#define Single_XS_H

#include "XS_Base.H"

namespace ATOOLS { class Blob_Data_Base; }

namespace EXTRAXS {

  class Single_XS: public XS_Base {
  public:

    // constructors
    Single_XS(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    Single_XS(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
	      const int scalescheme,const int kfactorscheme,const double scalefactor,
	      BEAM::Beam_Spectra_Handler *const beamhandler,
	      PDF::ISR_Handler *const isrhandler,
	      ATOOLS::Selector_Data *const selectordata);

    // member functions
    size_t Size() const;

    XS_Base *const operator[](const size_t i) const;

    double Differential2();

    void SetISR(PDF::ISR_Handler *const isrhandler);
    void SetTotal();

    double Differential(const double s,const double t,const double u);
    double operator()(const double s,const double t,const double u);
    void   AddPoint(const double xs);
    void   ResetMax(int);
    void   OptimizeResult();

    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create=false);
    void WriteOutXSecs(std::ofstream &outfile);

    bool OneEvent();

    ATOOLS::Blob_Data_Base *WeightedEvent(const int mode=0);

    template <class Process_Type>
    static Single_XS *GetProcess(const size_t nin,const size_t nout,
				 const ATOOLS::Flavour *flavours);

  };// end of class Single_XS

}// end of namespace EXTRAXS

#endif



