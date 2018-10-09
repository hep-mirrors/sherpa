// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_Ztautau : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Ztautau);
    double mT(const FourMomentum& p) {
      return sqrt(p.E2() - p.pz2());
    }

    double mT(const FourMomentum& p1, const FourMomentum& p2) {
      return sqrt((p1.Et()+p2.Et())*(p1.Et()+p2.Et()) - (p1+p2).perp()*(p1+p2).perp());
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs;
      IdentifiedFinalState Bare_Leptons(fs);
      IdentifiedFinalState Photons(Cuts::E > 0.1*GeV); // replace by E, run for 0.1, 0.01
      IdentifiedFinalState Neutrinos(fs);
      Bare_Leptons.acceptIdPair(15);
      Photons.acceptId(22);
      Neutrinos.acceptIdPair(16);
      DressedLeptons Dressed_1(Photons, Bare_Leptons, 0.1, Cuts::open(),1,1);//, cut, 1);
      DressedLeptons Dressed_2(Photons, Bare_Leptons, 0.2, Cuts::open(),1,1);//, cut, 1);
      declare(Photons, "Photons");
      declare(Bare_Leptons, "BareLeptons");
      declare(Dressed_1, "Dressed01");
      declare(Dressed_2, "Dressed02");
      declare(Neutrinos, "Neutrinos");

      // Book histograms - if 3 different histograms, index 0 - bare electrons
      // 1 - electrons dressed with dR 0.1, 2 - electron dressed with dR 0.2
      // Photon energy/pT diagrams
      _h_sum_E_gamma = bookHisto1D("Sum_E_Gamma_bare_t", 100, 0., 200.);
      _h_sum_pT_gamma = bookHisto1D("Sum_pT_Gamma_bare_t", 100, 0., 200.);
      _h_sum_E_gamma_rest_frame = bookHisto1D("Sum_E_Gamma_rest_frame_bare_t", 100, 0., 100.);
      _h_sum_pT_gamma_rest_frame = bookHisto1D("Sum_pT_Gamma_rest_frame_bare_t", 100, 0., 100.);
      _h_sum_E_gamma_coarse_rest_frame = bookHisto1D("Sum_E_Gamma_coarse_rest_frame_bare_t", 25, 0., 100.);
      _h_sum_pT_gamma_coarse_rest_frame = bookHisto1D("Sum_pT_Gamma_coarse_rest_frame_bare_t", 25, 0., 100.);

      // Photon angle diagrams
      // h_theta_phot - angle to reference direction, all photons
      // h_theta_hardest_phot - angle to reference direction, hardest photon
      // h_theta_closest_phot - smallest photon lepton angle in units 2ml/mZ
      _h_theta_phot = bookHisto1D("Theta_photon_bare_t", 50, 0., pi);
      _h_theta_closest_phot = bookHisto1D("Theta_photon_closest_bare_t", 50, 0., 10.);
      _h_theta_hardest_phot = bookHisto1D("Theta_hardest_photon_bare_t", 50, 0., pi);

      // invariant mass lepton antilepton
      _h_mll = bookHisto1D("mtt_bare", 100, 0., 150.);
      // invariant mass lepton antilepton photon - index 0 - hardest photon
      // index 1 = closest photon
      _h_mllg[0] = bookHisto1D("mttg_hardest_g_bare", 100, 0., 150.);
      _h_mllg[1] = bookHisto1D("mttg_closest_g_bare", 100, 0., 150.);
      // lepton pT and mT distributions (charged and all leptons)
      _h_pt_l = bookHisto1D("pT_t_charged_bare", 100, 0., 75.);
      _h_pt_ln = bookHisto1D("pT_t_all_bare", 100, 0., 75.);
      _h_mt_l = bookHisto1D("mT_t_charged_bare", 100, 0., 150.);
      _h_mt_ln = bookHisto1D("mT_t_all_bare", 100, 0., 150.);
      // mT lepton/antilepton system
      _h_mt_system_l = bookHisto1D("mT_t_system_bare", 100, 0., 150.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const Particles& Photons = 
	apply<IdentifiedFinalState>(event, "Photons").particles(Cuts::E > 0.1*GeV,cmpMomByE);
      const Particles& BareLeptons = apply<IdentifiedFinalState>(event, "BareLeptons").particles();
      const Particles& Neutrinos = apply<FinalState>(event, "Neutrinos").particles();
 
      /// @todo Do the event by event analysis here
      std::vector<DressedLepton> UseLeptons;
      Particles UsePhotons;
      for (size_t j = 0; j < BareLeptons.size(); j++) UseLeptons.push_back(BareLeptons[j]);
      UsePhotons = Photons;
      for (size_t j = 0; j < Neutrinos.size(); j++) UseLeptons.push_back(Neutrinos[j]);
      // Find Lorentz transform to go into back-to-back lepton frame
      // Store boosted lepton directions: Lepton_COM[0] is the reference
      // for h_theta_phot - either positively charged lepton or negatively
      // charged lepton if the other lepton is antineutrino
      FourMomentum UseLeptons_COM[2];
      LorentzTransform lep_lt;
      if (UseLeptons.size() > 1) {
	FourMomentum lep_sum = UseLeptons[0].momentum() + UseLeptons[1].momentum();
	lep_lt.setBetaVec(-lep_sum.betaVec());
	lep_lt.mkFrameTransformFromBeta(-lep_sum.betaVec());
	if (UseLeptons[0].charge() > 0.) {
	  UseLeptons_COM[0] = lep_lt.transform(UseLeptons[0].momentum());
	  UseLeptons_COM[1] = lep_lt.transform(UseLeptons[1].momentum());
	}
	else if (UseLeptons[0].charge() <= 0.) {
	  UseLeptons_COM[0] = lep_lt.transform(UseLeptons[1].momentum());
	  UseLeptons_COM[1] = lep_lt.transform(UseLeptons[0].momentum());
	}
      }
      // containers for hardest photon, closest photon to one of the leptons
      FourMomentum Hardest_Photon, Closest_Photon;
      // containers for sum of photon energies/pTs, smallest lepton-photon angle
      double sumPtgamma(0), sumEgamma(0),sumEgammaBoost(0),sumPtgammaBoost(0), small_angle(pi);
      // define boost into Z rest frame (Z momentum from IS momentat)
      LorentzTransform lt;
      const GenEvent* _genEvt = event.genEvent();
      const GenVertex* _sigVtx = _genEvt->signal_process_vertex();
      Particles in;
      foreach (const GenParticle* gp, particles(_sigVtx, HepMC::parents)) {
	in.push_back(Particle(gp));
      }
      Vector3 beta = (in[0].momentum()+in[1].momentum()).betaVec();

      // Need to do this first, otherwise beta is not set!
      lt.setBetaVec(-beta);
      lt.mkFrameTransformFromBeta(-beta);
      for (size_t i(0); i < UsePhotons.size(); ++i) {
	const Particle p(UsePhotons[i]);
	// Fill photon histograms
	sumPtgamma += p.pT();
	sumEgamma += p.E();
	FourMomentum phot_Z_boost = lt.transform(p.momentum());
	sumEgammaBoost += phot_Z_boost.E();
	sumPtgammaBoost += phot_Z_boost.pT();

	// define photon momentum in lept antilepton rest frame
	FourMomentum phot_l_boost = lep_lt.transform(p.momentum());
	// container for closest photon in lepton antilepton rest frame
	FourMomentum Closest_Photon_Boost;
	if (i == 0) {
	  // define initial closest photon, hardest photon
	  Hardest_Photon = p.momentum();
	  Closest_Photon = p.momentum();
	  Closest_Photon_Boost = phot_l_boost;
	  small_angle = min(UseLeptons_COM[0].angle(Closest_Photon_Boost.vector3()),UseLeptons_COM[1].angle(Closest_Photon_Boost.vector3()));
	}
	else {
	  if (p.momentum().E() > Hardest_Photon.E()) std::cout << "Not ordered by hardness" << std::endl;
	}
	// Fill (all, hardest photon) angle plots
	_h_theta_phot->fill(UseLeptons_COM[0].angle(phot_l_boost.vector3()),weight);
	if (i == 0) _h_theta_hardest_phot->fill(UseLeptons_COM[0].angle(phot_l_boost.vector3()),weight);
	// Update smallest photon-lepton angle
	if (i != 0) {
	  double angle = min(UseLeptons_COM[0].angle(phot_l_boost.vector3()),UseLeptons_COM[1].angle(phot_l_boost.vector3()));
	  if (angle < small_angle) small_angle = angle;
	}
      }
      // Fill photon energy/pT distributions
      _h_sum_pT_gamma->fill(sumPtgamma, weight);
      _h_sum_E_gamma->fill(sumEgamma, weight);
      _h_sum_E_gamma_rest_frame->fill(sumEgammaBoost, weight);
      _h_sum_pT_gamma_rest_frame->fill(sumPtgammaBoost, weight);
      _h_sum_E_gamma_coarse_rest_frame->fill(sumEgammaBoost, weight);
      _h_sum_pT_gamma_coarse_rest_frame->fill(sumPtgammaBoost, weight);
      if (UseLeptons.size() > 0) {
	double mZ = 91.1876; // check how to get from Rivet
	double ml = UseLeptons[0].mass(); // Lepton[0] is guaranteed to be charged lepton
	if (ml == 0.) { std::cout << "charged lepton massless" << std::endl; }
	if (UsePhotons.size() !=0 && ml != 0.) _h_theta_closest_phot->fill(small_angle*mZ/(2.*ml),weight);
	for (size_t i = 0; i < UseLeptons.size(); ++i) {
	  if (UseLeptons[i].charge() != 0.) {
	    _h_pt_l->fill(UseLeptons[i].momentum().pT(), weight);
	    _h_mt_l->fill(mT(UseLeptons[i].momentum()), weight);
	  }
	  _h_pt_ln->fill(UseLeptons[i].momentum().pT(), weight);
	  _h_mt_ln->fill(mT(UseLeptons[i].momentum()), weight);
	}
	if (UseLeptons.size() > 1) {
	  _h_mt_system_l->fill(mT(UseLeptons[0].momentum(),UseLeptons[1].momentum()),weight);
	  FourMomentum mom_ll = UseLeptons[0].momentum() + UseLeptons[1].momentum();
	  _h_mll->fill(sqrt(mom_ll*mom_ll), weight);
	  _h_mllg[0]->fill(sqrt((mom_ll+Hardest_Photon)*(mom_ll+Hardest_Photon)),weight);
	  _h_mllg[1]->fill(sqrt((mom_ll+Closest_Photon)*(mom_ll+Closest_Photon)),weight);
	}
      }
    
	

      /// @todo Do the event by event analysis here

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double incxs = crossSection();
      double sumweight = sumOfWeights();
      scale(_h_sum_E_gamma,1./sumweight);	       
      scale(_h_sum_pT_gamma,1./sumweight);	 
      scale(_h_sum_E_gamma_rest_frame,1./sumweight);
      scale(_h_sum_pT_gamma_rest_frame,1./sumweight);
      scale(_h_sum_E_gamma_coarse_rest_frame,1./sumweight);
      scale(_h_sum_pT_gamma_coarse_rest_frame,1./sumweight);
      scale(_h_mll,1./sumweight);	       
      scale(_h_mllg[0],1./sumweight);
      scale(_h_mllg[1],1./sumweight);
      scale(_h_pt_l,1./sumweight);
      scale(_h_pt_ln,1./sumweight);
      scale(_h_mt_l,1./sumweight);
      scale(_h_mt_ln,1./sumweight);
      scale(_h_mt_system_l,1./sumweight);
      scale(_h_theta_phot,1./sumweight);
      scale(_h_theta_closest_phot,1./sumweight);
      scale(_h_theta_hardest_phot,1./sumweight);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_sum_E_gamma;
    Histo1DPtr _h_sum_pT_gamma;
    Histo1DPtr _h_sum_E_gamma_rest_frame;
    Histo1DPtr _h_sum_pT_gamma_rest_frame;
    Histo1DPtr _h_sum_E_gamma_coarse_rest_frame;
    Histo1DPtr _h_sum_pT_gamma_coarse_rest_frame;
    Histo1DPtr _h_mll;
    Histo1DPtr _h_mllg[2];
    Histo1DPtr _h_pt_l;
    Histo1DPtr _h_pt_ln;
    Histo1DPtr _h_mt_l;
    Histo1DPtr _h_mt_ln;
    Histo1DPtr _h_mt_system_l;
    Histo1DPtr _h_theta_phot;
    Histo1DPtr _h_theta_closest_phot;
    Histo1DPtr _h_theta_hardest_phot;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Ztautau);


}
