#ifndef AHADIC_Formation_Gluon_Splitter_H
#define AHADIC_Formation_Gluon_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"
#include "ATOOLS/Org/Message.H"


namespace AHADIC {
  class Gluon_Splitter : public Splitter_Base {
  private:
    int    m_mode;
    std::vector<double> m_alpha;
    double m_arg, m_x, m_y, m_lastmass;
    bool   m_lastB, m_lastC;

    long int  m_kin_fails{ 0 };

    bool      MakeLongitudinalMomenta();
    void      CalculateLimits();
    bool      CalculateXY();

    bool      FillParticlesInLists();
    bool      CheckKinematics();
    bool      CheckConstituentKinematics(const ATOOLS::Vec4D & newmom11,
					 const ATOOLS::Vec4D & newmom12);
    Cluster * MakeCluster();
    void      ReplaceClusterWithHadron(const ATOOLS::Flavour & fl,ATOOLS::Vec4D & mom);
    void      UpdateSpectator(const ATOOLS::Vec4D & clumom);
    std::vector<double> variation_weights;
  public:
    Gluon_Splitter(std::list<Cluster *> * cluster_list,
		   Soft_Cluster_Handler * softclusters) :
      Splitter_Base(cluster_list,softclusters) {}
    ~Gluon_Splitter();

    void   Init(const bool & isgluon=true);
    double FragmentationFunction(double z, double zmin, double zmax, double alpha);
    double WeightFunction(const double & z,
			  const double & zmin=0.,const double & zmax=1.,
			  const unsigned int & cnt=0);

    inline void GetLast(double & lastmass, bool & isB, bool & isC) {
      lastmass = m_lastmass; isB = m_lastB; isC = m_lastC;
    }
    std::vector<double>& get_variationweights() {
      const auto wgts = m_ktselector.get_variationweights();
      for(int i{0}; i<wgts.size(); i++)
	variation_weights[i] *= wgts[i];

      return variation_weights;
    };

    void reset_variationweights() {
      // reset weights and Splitter_base weight
      std::fill(variation_weights.begin(), variation_weights.end(), 1);
      m_ktselector.reset_variationweights();
    };

    void z_rejected(const double, const double&, const double&, const double&,
		    const unsigned int&) override;

    void z_accepted(const double, const double&, const double&, const double&,
		    const unsigned int& ) override;

  };
}

#endif
