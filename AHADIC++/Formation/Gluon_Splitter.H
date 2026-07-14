#ifndef AHADIC_Formation_Gluon_Splitter_H
#define AHADIC_Formation_Gluon_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"


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

  public:
    Gluon_Splitter(std::list<Cluster *> * cluster_list,
		   Soft_Cluster_Handler * softclusters,
		   Flavour_Selector     * flavourselector,
		   KT_Selector          * ktselector,
		   Hadronisation_Reweighting   * reweighting) :
      Splitter_Base(cluster_list,softclusters,flavourselector,ktselector,
		    reweighting) {}
    ~Gluon_Splitter();

    void   Init();
    double FragmentationFunction(double z, double zmin, double zmax, double a);
    double FragmentationFunctionProb(double z, double zmin, double zmax, double a);
    double WeightFunction(const double & z,
			  const double & zmin=0.,const double & zmax=1.,
			  const unsigned int & cnt=0);
    inline void GetLast(double & lastmass, bool & isB, bool & isC) {
      lastmass = m_lastmass; isB = m_lastB; isC = m_lastC;
    }
    void AcceptSplittingWeights() override {
      p_reweighting->AcceptGluonSplittingWeights();
    }

    void ResetSplittingWeights() override {
      p_reweighting->ResetGluonSplittingWeights();
    }

    int RecordSite() const override { return 0; } // OUTPUT

    void ZRejected(const double, const double&, const double&, const double&,
		   const unsigned int&) override;

    void ZAccepted(const double, const double&, const double&, const double&,
		   const unsigned int& ) override;

  };
}

#endif
