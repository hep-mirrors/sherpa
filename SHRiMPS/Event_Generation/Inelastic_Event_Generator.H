#ifndef SHRIMPS_Main_Inelastic_Event_Generator_H
#define SHRIMPS_Main_Inelastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Event_Generation/Ladder_Generator.H"
#include "SHRiMPS/Event_Generation/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Beam_Remnants/Rescatter_Handler.H"
#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace SHRIMPS {
  class Inelastic_Event_Generator : public Event_Generator_Base {
  private:
    Sigma_Inelastic         m_sigma;
    Parton_Luminosity       m_luminosity;
    Ladder_Generator        m_laddergenerator;
    Beam_Remnant_Handler  * p_beams;
    Rescatter_Handler       m_rescatterhandler;

    double     m_Bmin, m_Bmax, m_B;
    int        m_Bsteps;
    double     m_deltaB;
    bool       m_isUE, m_weighted, m_first, m_done;
    int        m_Nladders_fix,m_Nladders, m_Ngen, m_Nprim, m_Nsec;
    double     m_kt2fac, m_difffac;
    Ladder *   p_ladder;

    int	 m_test, m_output;

    long int m_connectblobs, m_laddercols, m_updatecols;
    bool m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;


    void FillCrossSections();
    void InitInelasticEvent();
    void FixEikonalAndImpact();
    int  AddScatter(ATOOLS::Blob_List * blobs,const double & xsec=1.);
    bool CreateBlob(ATOOLS::Blob_List * blobs,const double & xsec);

    void TestNumberOfLadders(Omega_ik * eikonal,const double & B);
  public:
    Inelastic_Event_Generator(Beam_Remnant_Handler * beams);
    ~Inelastic_Event_Generator();

    bool DressShowerBlob(ATOOLS::Blob * blob);
    int  GenerateEvent(ATOOLS::Blob_List * blobs,const bool & isUE=false);

    inline double XSec(Omega_ik * eikonal=NULL) {
      return m_sigma.XSec(eikonal);
    }

    double Smin() const;
    bool   IsLastRescatter() const;

    void Reset();

    Omega_ik * GetEikonal() const { return p_eikonal; }
    inline double TMax() const {return p_ladder->That(); }
    inline int    NLadders() const {return (m_Nladders+m_Nsec); }
  };
}
#endif
