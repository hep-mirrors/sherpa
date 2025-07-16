#ifndef SHRIMPS_Ladders_Ladder_Generator_LDC_H
#define SHRIMPS_Ladders_Ladder_Generator_LDC_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Ladder_Generator_LDC : public Ladder_Generator_Base {
  private:
    Continued_PDF  * p_pdf[2];
    Sigma_Partonic   m_partonic;
    double           m_kt2min;
    double           m_Ebeam[2], m_y[2], m_kt2[2], m_zp[2], m_zm[2];
    ATOOLS::Vec4D    m_Pbeam[2], m_k[2];
    ATOOLS::Flavour  m_flavs[2];
    std::set<double> m_zps, m_zms;
    
    bool   FixInitialPartons();
    bool   FillZs();
    size_t SelectN(const double & a,const double & b);
    void   CalculateWeight();
    void   Output();
  public:
    Ladder_Generator_LDC();
    ~Ladder_Generator_LDC() {}

    void Initialise() { m_partonic.Initialise(); }
    Ladder * operator()(const ATOOLS::Vec4D & pos);

    const  double Weight() const { return m_weight; }
    inline void SetRemnantHandler(Remnant_Handler * remnants) {
      for (size_t beam=0;beam<2;beam++) {
	p_pdf[beam] = remnants->GetPDF(beam);
	m_partonic.SetPDF(beam,p_pdf[beam]);
      }
    }
  };
}
#endif
