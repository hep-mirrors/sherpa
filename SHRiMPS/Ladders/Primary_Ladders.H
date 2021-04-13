#ifndef SHRIMPS_Ladders_Primary_Ladders_H
#define SHRIMPS_Ladders_Primary_Ladders_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Primary_Ladders {
  private:
    Ladder_Generator_Base * p_laddergenerator;
    Colour_Generator        m_colourgenerator;
    double                  m_Ecms, m_E[2];
    std::list<Ladder *>     m_ladders;

    bool m_test;
    std::map<std::string,ATOOLS::Histogram * > m_histos;

    bool   IsAllowed(Ladder * ladder);
    void   Add(Ladder * ladder);
    void   FillAnalysis(Ladder * ladder,const std::string & tag);
  public:
    Primary_Ladders();
    ~Primary_Ladders();
    
    bool operator()(Omega_ik * eikonal,const double & B,const size_t & N);
    void Initialise();
    void Test();
    void Reset();

    inline Colour_Generator    * GetColourGenerator()  { return &m_colourgenerator; }
    inline std::list<Ladder *> * GetLadders()          { return &m_ladders; }
    inline void                  SetRemnantHandler(Remnant_Handler * remnants) {
      p_laddergenerator->SetRemnantHandler(remnants);
    }
    inline const double Ymax()   const { return p_laddergenerator->Ymax(); }
    inline const double KT2min() const { return p_laddergenerator->KT2min(); }
  };
}
#endif
