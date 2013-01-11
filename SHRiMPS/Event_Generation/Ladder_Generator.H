#ifndef SHRIMPS_Main_Ladder_Generator_H
#define SHRIMPS_Main_Ladder_Generator_H

#include "SHRiMPS/Event_Generation/Ladder.H"
#include "SHRiMPS/Event_Generation/Simple_MEs.H"
#include "SHRiMPS/Event_Generation/Initial_State.H"
#include "SHRiMPS/Event_Generation/Final_State.H"
#include "SHRiMPS/Cross_Sections/Parton_Luminosity.H"

#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Histogram_2D.H"

namespace SHRIMPS {
  class Ladder_Generator {
  private:
    Initial_State     m_IS;
    Final_State       m_FS;
    Simple_MEs        m_xsec;

    double            m_originalY, m_cutoffY;
    double            m_pplus, m_pminus;

    Omega_ik        * p_eikonal;
    double            m_B,m_b1,m_b2;
    double            m_Ycms, m_Ymin, m_Ymax, m_deltaY, m_Yrange;
    double            m_maxkt2;
    double            m_prodkt2;
    int               m_Nprops;
    ATOOLS::Vec4D     m_pos;

    Ladder          * p_ladder;

    long int          m_N1h,m_N1,m_N8,m_resc0,m_resc1,m_wtover1,m_wtover2,m_wt;
    bool              m_output;
    
    int               m_Nprim,m_Nsec,m_Ndd_p,m_Ndd_s,m_Nsd_p,m_Nsd_s,m_Ncep_p,m_Ncep_s;

    std::map<std::string, ATOOLS::Histogram * > m_histograms;
    std::map<std::string, ATOOLS::Histogram_2D * > m_histogram2ds;

    double InitialiseLadder(ATOOLS::Particle * part1,ATOOLS::Particle * part2,
			    const bool & rescatter);
    bool   FixFirstOutgoings();
    bool   Fix2To2Outgoing(const ATOOLS::Vec4D & inmom1,
			   const ATOOLS::Vec4D & inmom2,
			   ATOOLS::Vec4D & outmom1,
			   ATOOLS::Vec4D & outmom2,bool keep);
    double Weight(const double & isweight); 
    void   Analyse(const bool & isprimary);
  public:
    Ladder_Generator() : p_ladder(0), m_output(false) {}
    Ladder_Generator(Parton_Luminosity * lumi,const int & test=false);
    ~Ladder_Generator();

    void     InitCollision(Omega_ik * eikonal,const double & B);
    Ladder * operator()(ATOOLS::Particle * part1,ATOOLS::Particle * part2,
			const bool & rescatter,
			const bool & first=false,const bool & weighted=false);
    //     Ladder * MakeLadder(ATOOLS::Blob * blob);

//    inline void SetNPrim(const int & N) { m_FS.SetQ02eff(N); }
    inline void SetNPrim(const int & N) { m_FS.SetNLad(N); }
    inline void SetLadderGeneration(const int & N) { m_FS.SetQ02eff(N); }
    //     inline const double & B1()   const { return m_b1; }
    //     inline const double & B2()   const { return m_b2; }
    inline const double & B1()   const { return m_IS.B1(); }
    inline const double & B2()   const { return m_IS.B2(); }
    inline const double & B()    const { return m_B; }
    inline const double & Smin() const { return m_IS.Smin(); }
  };
}
#endif
