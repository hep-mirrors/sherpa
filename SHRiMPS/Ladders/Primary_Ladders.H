#ifndef SHRIMPS_Ladders_Primary_Ladders_H
#define SHRIMPS_Ladders_Primary_Ladders_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "ATOOLS/Math/Histogram.H"
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_SD.H"

namespace SHRIMPS {
  class Primary_Ladders {
  private:
    Ladder_Generator_Base * p_laddergenerator;
    Colour_Generator        m_colourgenerator;
    double                  m_Ecms, m_E[2];
    std::list<Ladder *>     m_ladders;

    bool m_test;
    long unsigned int n_calls, n_start, n_gen;
    std::map<std::string,ATOOLS::Histogram * > m_histos;

    bool   IsAllowed(Ladder * ladder);
    void   Add(Ladder * ladder);
    void   FillAnalysis(Ladder * ladder,const std::string & tag);

    Sigma_Elastic * p_sigma_el;
    Sigma_SD * p_sigma_sd;
  public:
    Primary_Ladders(Sigma_Elastic * sigma_el,Sigma_SD * sigma_sd);
    ~Primary_Ladders();
    
    bool operator()(Omega_ik * eikonal,const double & B,const size_t & N);
    void Initialise(Remnant_Handler * remnants);
    void Test();
    void Reset();

    inline void SetMaxEnergies(const double & E1,const double & E2) {
      m_E[0] = E1; m_E[1] = E2;
    }
    inline Colour_Generator    * GetColourGenerator()  { return &m_colourgenerator; }
    inline std::list<Ladder *> * GetLadders()          { return &m_ladders; }
    inline const double Yhat()   const { return p_laddergenerator->Yhat(); }
    inline const double Ymax()   const { return p_laddergenerator->Ymax(); }
    inline const double KT2min() const { return p_laddergenerator->KT2min(); }
  };
}
#endif
