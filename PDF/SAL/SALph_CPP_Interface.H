/*
 * SALph_CPP_Interface.H
 *
 *  Created on: 1 Jul 2021
 *      Author: pmeinzinger
 */

#ifndef PDF_SAL_SALPH_CPP_INTERFACE_H_
#define PDF_SAL_SALPH_CPP_INTERFACE_H_

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

namespace PDF {
  class SALph_CPP_Interface : public PDF_Base {
  private:
    //std::string m_set,m_path;
    float   m_d,m_u,m_s,m_c,m_b,m_t,m_g;
    int m_nf;
  public:
    SALph_CPP_Interface(const ATOOLS::Flavour);
    PDF_Base * GetCopy();

    void   CalculateSpec(const double&, const double&);
    double GetXPDF(const ATOOLS::Flavour&);
    double GetXPDF(const kf_code&, bool);

    double GetXMin()   { return 1.e-5; } // TODO: x_min is actually dependent on the quark flavour
    double GetXMax()   { return 0.9999; }
    double GetQ2Min()  { return 2.; }
    double GetQ2Max()  { return 8.e4; }
  };
}

#endif /* PDF_SAL_SALPH_CPP_INTERFACE_H_ */
