#ifndef PDF_MRST_H
#define PDF_MRST_H

#include <iostream>
#include "PDF_Base.H"
#include "mrst99.h"

namespace PDF {
  class PDF_MRST : public PDF_Base {
    // only one MRST Library
    static c_mrst * proton;
    // but local content store
    s_partoncontent content;
    int anti;
    std::vector<APHYTOOLS::Flavour> partons;
  public:
    // initializes
    PDF_MRST(int isanti=0); 

    // calculates distributions
    void Calculate(const double x, const double Q2) {
      //std::cout<<" in PDF_MRST::Calculate("<<x<<','<<Q2<<")"<<std::endl;
      proton->mrst99(x,Q2,1);
      content=proton->cont;
    }
    // returns x times the parton distribution corresponding to the given flavour
    double GetXPDF(const APHYTOOLS::Flavour & fl) {
      //    kf::code kfc=fl.kfcode(); does not hold anti information
      int kfc=anti*int(fl);
      switch (kfc) {
      case  APHYTOOLS::kf::d : return (content.dnv + content.dsea);
      case -APHYTOOLS::kf::d : return content.dsea; 
      case  APHYTOOLS::kf::u : return (content.upv + content.usea);
      case -APHYTOOLS::kf::u : return content.usea; 
      case  APHYTOOLS::kf::s :
      case -APHYTOOLS::kf::s : return content.str;
      case  APHYTOOLS::kf::c : 
      case -APHYTOOLS::kf::c : return content.chm;
      case  APHYTOOLS::kf::b : 
      case -APHYTOOLS::kf::b : return content.bot;
      case APHYTOOLS::kf::gluon : 
      case -APHYTOOLS::kf::gluon :return content.glu; // pseudo anti gluon for anti-proton
      default: return 0.;
      }
    };
    PDF_Base * GetCopy() {
      if (anti==-1) return new PDF_MRST(1);  // anti-proton
      return new PDF_MRST; // proton
    }

    // return domain for x and Q of the PDF
    double             GetXMin()  { return MRST99::xmin; };
    double             GetXMax()  { return MRST99::xmax; };
    double             GetQ2Min() { return MRST99::qsqmin; };
    double             GetQ2Max() { return MRST99::qsqmax; };
    APHYTOOLS::Flavour Beam()     { 
      if (anti==-1) return APHYTOOLS::Flavour(APHYTOOLS::kf::p_plus).bar();
      return               APHYTOOLS::Flavour(APHYTOOLS::kf::p_plus);
    }
    std::vector<APHYTOOLS::Flavour> Partons() { return partons; }
  };
}  // namespace PDF


#endif // PDF_MRST_H
