// HADRON_RESCATTERING/Main/Deuteron_Production.H
#ifndef DEUTERON_COALESCENCE_H
#define DEUTERON_COALESCENCE_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Flavour.H"
#include "HADRON_RESCATTERING/XSecs/BaryonBaryon.H"
#include <string>
#include <vector>

namespace HADRON_RESCATTERING {

  // HardDisk, Exponential: coalescence in k = |p_p - p_n| (pair frame) with
  // coalescence momentum p_coal, and in the impact parameter b of the pair
  // with radius b_coal (no spatial condition for b_coal <= 0).
  // CrossSectionModel: arXiv:1504.07242.
  enum class CoalescenceModel { HardDisk, Exponential, CrossSectionModel };

  // Deuteron formation channels of Table I in arXiv:1504.07242, written for
  // nucleons; for antinucleon pairs all final state flavours are conjugated.
  enum class Deuteron_Channel {
    pn_dgamma, pn_dpi0, pp_dpiplus, nn_dpiminus,
    pn_dpi0pi0, pn_dpipluspiminus, pp_dpipluspi0, nn_dpiminuspi0,
    none
  };

  struct Channel_Probability {
    Deuteron_Channel channel;
    double           sigma;  // in microbarn
    double           prob;   // P = sigma/sigma_0
  };

  class Deuteron_Coalescence {
  private:
    double           m_b_coal;
    double           m_p_coal;
    double           m_inv_sigma0;   // 1/sigma_0 in 1/barn
    CoalescenceModel m_model;
    BaryonBaryon     m_bb;           // for cross section model

    double Sigma(const Deuteron_Channel & channel,
		 ATOOLS::Particle * A,ATOOLS::Particle * B);
    double Threshold(const Deuteron_Channel & channel) const;
    void   Isotropic(const double & p,ATOOLS::Vec3D & dir) const;
    void   TwoBody(const ATOOLS::Vec4D & P,const double & m1,const double & m2,
		   ATOOLS::Vec4D & p1,ATOOLS::Vec4D & p2) const;
  public:
    Deuteron_Coalescence(double b_coal, double p_coal, double inv_sigma0,
			 CoalescenceModel model);

    // k = |p_A - p_B| in the pair c.m. frame [GeV]
    static double K(const ATOOLS::Vec4D & pA,const ATOOLS::Vec4D & pB);
    // squared impact parameter [fm^2] of the straight trajectories of A and B
    // in their c.m. frame: separation transverse to the relative momentum
    static double ImpactParameter2(ATOOLS::Particle * A,ATOOLS::Particle * B);

    static std::string Name(const Deuteron_Channel & channel);
    static std::vector<Deuteron_Channel> Channels();
    // Flavours of the final state, deuteron first, conjugated if anti.
    static std::vector<ATOOLS::Flavour> Products(const Deuteron_Channel & channel,
						 const bool & anti);
    // Returns 0 if A and B are not a (anti)nucleon pair, else 1 for nucleons
    // and -1 for antinucleons.
    static int PairType(ATOOLS::Particle * A,ATOOLS::Particle * B);
    static bool Matches(const Deuteron_Channel & channel,
			ATOOLS::Particle * A,ATOOLS::Particle * B);

    // Cross section model: sigma and P for all channels open for this pair.
    // For the HardDisk/Exponential models only pn -> d gamma is open.
    std::vector<Channel_Probability> Probabilities(ATOOLS::Particle * A,
						   ATOOLS::Particle * B,
						   const double & b2_fm2=0.);
    // Paper's prescription: one random number per channel; if more than one
    // channel succeeds, select by relative cross section.
    Deuteron_Channel SelectChannel(ATOOLS::Particle * A,ATOOLS::Particle * B,
				   const double & b2_fm2=0.);
    // Final state momenta (lab frame) in the order of Products(); false if
    // the channel is kinematically closed.
    bool FinalState(const Deuteron_Channel & channel,
		    const ATOOLS::Vec4D & pA,const ATOOLS::Vec4D & pB,
		    std::vector<ATOOLS::Vec4D> & moms) const;
    ATOOLS::Blob * BuildBlob(const Deuteron_Channel & channel,
			     ATOOLS::Particle * A,ATOOLS::Particle * B);
    ATOOLS::Blob * operator()(ATOOLS::Particle* A,
			      ATOOLS::Particle* B,
			      const double& b2_fm2=0.);

    void Tests(const std::string & dir);

    inline const double & InvSigma0() const { return m_inv_sigma0; }
    inline const CoalescenceModel & Model() const { return m_model; }
  };
}
#endif
