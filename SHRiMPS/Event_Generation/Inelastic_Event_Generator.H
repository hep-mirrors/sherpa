#ifndef SHRIMPS_Main_Inelastic_Event_Generator_H
#define SHRIMPS_Main_Inelastic_Event_Generator_H

#include "SHRiMPS/Event_Generation/Event_Generator_Base.H"
#include "SHRiMPS/Event_Generation/Ladder_Generator.H"
#include "SHRiMPS/Event_Generation/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Cross_Sections/Parton_Luminosity.H"
#include "SHRiMPS/Beam_Remnants/Rescatter_Handler.H"
#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace SHRIMPS {
  class Inelastic_Event_Generator : public Event_Generator_Base {
  private:
    Sigma_Inelastic       * p_sigma;
    std::list<Omega_ik *> * p_eikonals;
    Parton_Luminosity       m_luminosity;
    Ladder_Generator        m_laddergenerator;
    Beam_Remnant_Handler  * p_beams;
    Rescatter_Handler       m_rescatterhandler;

    double     m_Bmin, m_Bmax, m_B;
    int        m_Bsteps;
    double     m_deltaB;
    bool       m_isUE, m_weighted, m_first, m_done;
    int        m_Nladders, m_Ngen, m_Nprim, m_Nsec;
    double     m_kt2fac, m_difffac;
    Ladder *   p_ladder;

    int	 m_test, m_output;

    long int m_connectblobs, m_laddercols, m_updatecols;
    bool m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;


    void FillGrids();
    void InitInelasticEvent(const bool & isUE=false,
			    const bool & weighted=false);
    int  AddScatter(ATOOLS::Blob_List * blobs,const double & xsec=1.);
    bool CreateBlob(ATOOLS::Blob_List * blobs,const double & xsec);

    void TestNumberOfLadders(Omega_ik * eikonal,const double & B);
  public:
    Inelastic_Event_Generator() : m_test(0), m_output(0) {}
    Inelastic_Event_Generator(Sigma_Inelastic * sigma,
			      std::list<Omega_ik *> * eikonals,
			      Beam_Remnant_Handler * beams,
			      const int & test);
    ~Inelastic_Event_Generator();

    bool DressShowerBlob(ATOOLS::Blob * blob);
    int  InelasticEvent(ATOOLS::Blob_List * blobs,const double & xsec=1.,
			const bool & isUE=false,const bool & weighted=false);

    double Smin() const;
    bool   IsLastRescatter() const;

    void Reset();

    Omega_ik * GetEikonal() const { return p_eikonal; }
    inline const double & XSec()   const { return p_sigma->XSec(); }
  };
}
#endif
