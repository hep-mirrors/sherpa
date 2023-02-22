#ifndef SHRIMPS_Ladders_Ladder_Generator_QE_H
#define SHRIMPS_Ladders_Ladder_Generator_QE_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Histogram.H"
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_SD.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"


namespace SHRIMPS {
  class Ladder_Generator_QE : public Ladder_Generator_Base {
  private:
    ATOOLS::Vec4D   m_lastk;
    void   FillGluons();
    void   SelectPropagatorColours();
    void   CalculateWeight();
    void   FixEmissionsKinematics_elastic();
    Sigma_Elastic * p_sigma_el;
    Sigma_SD * p_sigma_sd;
    Sigma_Partonic   m_partonic;
    Continued_PDF  * p_pdf[2];
    double m_abs_t;
    ATOOLS::Flavour m_beam1, m_beam2, m_ladder_ends[2];                
    ATOOLS::Vec4D   m_p1, m_p2, m_Pbeam1, m_Pbeam2;
    double          m_pl12, m_pl22, m_pl1, m_pl2;
    int             m_sign1;
    double          m_fraction;
    int m_interactionType;
    std::map <ATOOLS::Flavour,double> m_PfractionMap;


  public:
    Ladder_Generator_QE();
    ~Ladder_Generator_QE();
    int GetIntType() {return 2;};
    void FixEmissionsKinematics_SD(int mode);
    
    Ladder * operator()(const ATOOLS::Vec4D & pos,Sigma_Elastic * sigma_el,Sigma_SD * sigma_sd);  
  };
}
#endif
