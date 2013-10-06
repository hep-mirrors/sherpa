#ifndef SHRIMPS_Main_Final_State_H
#define SHRIMPS_Main_Final_State_H

#include "SHRiMPS/Event_Generation/Ladder.H"
#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  class Final_State {
  private:
    ktform::code             m_ktform;
    ordering::code           m_ordering;
    resc_ktmin::code         m_resc_ktmin;
    resc_nosing::code        m_resc_nosing;
    int                      m_ktmin_mode;
    MODEL::Strong_Coupling * p_alphaS;
    double                   m_Q02,m_Q12,m_Q02eff,m_QN2,m_q02min;
    double                   m_kt1cut2,m_kt1min2,m_kt1max2;
    double                   m_d2, m_kdiff, m_Ylimit;
    int                      m_test;
    int                      m_nprimlad;

    Ladder                 * p_ladder;
    Omega_ik               * p_eikonal;
    double                   m_B,m_b1,m_b2,m_shat,m_Deltay;

    double                   m_lastwt,m_recombwt;

    LadderMap              * p_emissions;
    LadderMap::iterator      m_plusiter;
    LadderMap::iterator      m_minusiter;
    TPropList              * p_props;
    TPropList::iterator      m_propiter;
    ATOOLS::Vec4D            m_pprop, m_mprop;
    ATOOLS::Vec4D            m_k0,m_k1,m_k2,m_q0,m_q2;
    double                   m_q01_2,m_q12_2;

    bool     m_output;
    std::map<std::string, ATOOLS::Histogram * > m_histomap;
    long int m_rej_negkt1,m_rej_nokt2,m_rej_nokt0,m_rej_nophi0,m_rej_noy0;
    long int m_rej_offshell,m_rej_nofit,m_rej_order,m_trials,m_ys,m_alphaS;
    long int m_dir0,m_dir1,m_rej_order0,m_rej_order1,m_cols;
    long int m_singexit,m_firstsing;

    double GenerateEmissions();
    bool   FirstSinglet(const double & y1,const double & y2,
			const double & sup,const int & nbeam);
    bool   OneEmission(LadderMap::iterator & split,LadderMap::iterator & spect,
		       double & kt2);
    bool   SelectSplitterAndSpectator(LadderMap::iterator & split,
				      LadderMap::iterator & spect);
    void   SwapSplitterAndSpectator(LadderMap::iterator & split,
				    LadderMap::iterator & spect);
    bool   FixPropColours(const LadderMap::iterator & split,
			  const LadderMap::iterator & spect);
    bool   TryEmission(double & kt2,const bool & dir);
    bool   MomViolation(ATOOLS::Vec4D & k_0, ATOOLS::Vec4D & k_1,
			ATOOLS::Vec4D & k_2,bool dir);
    double KT2integral(const double & kt2max,const double & kt2min,
		       const double & q02,const double & ktexpo);
    double FKT2(const double &, const double &, const double &);
    bool IsOrdered(const bool & dir,ATOOLS::Vec4D & k_0,
			 ATOOLS::Vec4D & k_1,ATOOLS::Vec4D & k_2,
			 const double & t);
    double SuppressionTerm(const double & q02,const double & q12);
    double Q02MinEstimate(const double, const double);
    void   InitHistograms();
    void   OutputHistograms();
  public:
    Final_State() : m_output(false) {}
    Final_State(const int & test);
    inline ~Final_State() { if (m_output) OutputHistograms(); }

    double operator()(Ladder * ladder,const double & Deltay,
		      const bool & firstladder,const bool & firstattempt);
    inline MODEL::Strong_Coupling * GetAS(){ return p_alphaS; }
    inline void Init(Omega_ik * eikonal,const double & b1,const double & b2,const double & B) { 
      p_eikonal=eikonal;m_b1=b1;m_b2=b2;m_B=B;
    }
    inline void   SetQ02eff(const int & N) { 
      m_Q02eff = m_Q02+(N-1)*m_QN2;
    }
    inline void SetNLad(const int & N) { m_nprimlad = N; }
    double AlphaS(const double & kt2) { return (*p_alphaS)(kt2); }
    inline double AlphaSMax() const { return p_alphaS->MaxValue(); }
    double SelectKT2(const double & kt2max,const double & kt2min,
		     const double & mu2,const double & expo=1.);
    double Q02(const double & y=0.);
    double Saturation(const double & y=0.);
    const double & KTMin2() const { return m_kt1min2; }
  };
}
#endif
