#ifndef Single_XS_H
#define Single_XS_H

#include "XS_Base.H"

namespace ATOOLS { class Blob_Data_Base; }

namespace EXTRAXS {

  class Single_XS : public XS_Base {
  public:

    // constructors
    Single_XS(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    Single_XS(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
	      const int scalescheme,const int kfactorscheme,const double scalefactor,
	      BEAM::Beam_Spectra_Handler *const beamhandler,
	      PDF::ISR_Handler *const isrhandler,
	      ATOOLS::Selector_Data *const selectordata);

    // member functions
    size_t Size();

    double Differential2();

    void SetISR(PDF::ISR_Handler *const isrhandler);
    void SetTotalXS();

    double Differential(const double s,const double t,const double u);
    double operator()(const double s,const double t,const double u);
    void   AddPoint(const double xs);

    bool CalculateTotalXSec(const std::string &resultpath);
    void WriteOutXSecs(std::ofstream &outfile);

    bool   OneEvent();

    ATOOLS::Blob_Data_Base *WeightedEvent();

  };// end of class Single_XS

}// end of namespace EXTRAXS

#endif



