#ifndef PDF_MRST99_H
#define PDF_MRST99_H

#include <iostream>
#include "PDF_Base.H"
#include "mrst99.h"

namespace PDF {
  class PDF_MRST99 : public PDF_Base {
  private:
    static c_mrst * p_proton;
    s_partoncontent m_content;
    std::string     m_path;
    int             m_set, m_anti;
  public:
    PDF_MRST99(const APHYTOOLS::Flavour,const int,const std::string); 
    void               Output();
    void               Calculate(const double, const double);
    double             GetXPDF(const APHYTOOLS::Flavour & fl); 
    PDF_Base         * GetCopy();
    double             GetXMin()  { return MRST99::xmin; };
    double             GetXMax()  { return MRST99::xmax; };
    double             GetQ2Min() { return MRST99::qsqmin; };
    double             GetQ2Max() { return MRST99::qsqmax; };
    APHYTOOLS::Flavour Bunch()    { return m_bunch; } 
    std::vector<APHYTOOLS::Flavour> Partons() { return m_partons; }
  };
  /*!
    \class PDF_MRST99
    \brief An interface to the C++-version of the MRST99 sets.
    
    This is an interface of the SHERPA package to the 
    <A HREF="http://durpdg.dur.ac.uk/hepdata/mrstc++.html>C++-version of the MRST99<\A>
    sets of parton distribution functions by Jeppe Andersen. 
  */
  /*!
    \var static c_mrst * PDF_MRST99::proton
    A pointer to the C++ object of the MRST99 set.
  */
}  // namespace PDF


#endif // PDF_MRST99_H
