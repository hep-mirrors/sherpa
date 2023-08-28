#ifndef AHADIC_Decays_Cluster_Splitter_H
#define AHADIC_Decays_Cluster_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"

namespace AHADIC {
  class Cluster_Splitter : public Splitter_Base {
  private:
    int    m_mode, m_defmode, m_beammode;

    // we have four cluster types, L, H, D, B each is getting a std::vector of
    // alpha/beta/gamme values
    std::vector<double> m_alpha[4], m_beta[4], m_gamma[4];
    int m_a[2], m_b[2], m_c[2];
    std::vector<double> m_kt02;
    double threshold_fac;
    double m_R2[2], m_mean[2], m_sigma[2], m_m2min[2], m_masses;
    std::vector<double> variation_weights;

    Cluster *       p_out[2];
    ATOOLS::Vec4D   m_cms, m_mom[2], m_newmom[2];
    ATOOLS::Flavour m_fl[2];
    double          m_mass2[2];

    bool   m_output;

    bool   MakeLongitudinalMomenta();
    bool   MakeLongitudinalMomentaZ();
    bool   MakeLongitudinalMomentaZSimple();
    bool   MakeLongitudinalMomentaMass();
    bool   MakeLongitudinalMomentaMassSimple();
    double DeltaM(const size_t & cl);
    void   FixCoefficients();
    void   CalculateLimits();
    bool   FillParticlesInLists();
    size_t MakeAndCheckClusters();
    void   MakeNewMomenta(size_t shuffle);
    void   FillHadronAndDeleteCluster(size_t i);
    void   UpdateAndFillCluster(size_t i);
    bool   RecalculateZs();
    bool   CheckKinematics() { return true; }
    Cluster * MakeCluster(size_t i);
    double FragmentationFunction(double z, double zmin, double zmax,
				 double alpha, double beta, double gamma,
				 double kt02);
  public:
    Cluster_Splitter(std::list<Cluster *> * cluster_list,
		     Soft_Cluster_Handler * softclusters,
		     Flavour_Selector     * flavourselector,
		     KT_Selector          * ktselector);
    void   Init(const bool & isgluon=false);
    double WeightFunction(const double & z,
			  const double & zmin,const double & zmax,
			  const unsigned int & cnt);
    void   SetOutput(const bool & out) { m_output = out; }

    std::vector<double>& get_variationweights() {
      return variation_weights;
    };

    void reset_variationweights(int nvars) {
      variation_weights.resize(nvars);
      // reset weights and Splitter_base weight
      std::fill(variation_weights.begin(), variation_weights.end(), 1);
    };

    void z_rejected(const double wgt, const double & z,
		    const double & zmin,const double & zmax,
		    const unsigned int & cnt) override;

    void z_accepted(const double wgt, const double & z,
		    const double & zmin,const double & zmax,
		    const unsigned int & cnt) override;
  };
}

#endif
