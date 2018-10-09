// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_Zee : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Zee);
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

      FinalState fs;
      IdentifiedFinalState Bare_Leptons(fs);
      IdentifiedFinalState Photons(Cuts::E > 0.1*GeV);
      IdentifiedFinalState Neutrinos(fs);
      Bare_Leptons.acceptIdPair(11);
      Photons.acceptId(22);
      Neutrinos.acceptIdPair(12);
      Cut cut = (Cuts::abseta < 2.4 && Cuts::E > 0.5*GeV);
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
      _h_sum_E_gamma[0] = bookHisto1D("Sum_E_Gamma_bare_e", 100, 0., 200.);
      _h_sum_E_gamma[1] = bookHisto1D("Sum_E_Gamma_dressed_1", 100, 0., 200.);
      _h_sum_E_gamma[2] = bookHisto1D("Sum_E_Gamma_dressed_2", 100, 0., 200.);
      _h_sum_pT_gamma[0] = bookHisto1D("Sum_pT_Gamma_bare_e", 100, 0., 200.);
      _h_sum_pT_gamma[1] = bookHisto1D("Sum_pT_Gamma_dressed_1", 100, 0., 200.);
      _h_sum_pT_gamma[2] = bookHisto1D("Sum_pT_Gamma_dressed_2", 100, 0., 200.);
      _h_sum_E_gamma_rest_frame[0] = bookHisto1D("Sum_E_Gamma_rest_frame_bare_e", 100, 0., 100.);
      _h_sum_E_gamma_rest_frame[1] = bookHisto1D("Sum_E_Gamma_rest_frame_dressed_1", 100, 0., 100.);
      _h_sum_E_gamma_rest_frame[2] = bookHisto1D("Sum_E_Gamma_rest_frame_dressed_2", 100, 0., 100.);
      _h_sum_pT_gamma_rest_frame[0] = bookHisto1D("Sum_pT_Gamma_rest_frame_bare_e", 100, 0., 100.);
      _h_sum_pT_gamma_rest_frame[1] = bookHisto1D("Sum_pT_Gamma_rest_frame_dressed_1", 100, 0., 100.);
      _h_sum_pT_gamma_rest_frame[2] = bookHisto1D("Sum_pT_Gamma_rest_frame_dressed_2", 100, 0., 100.);
      _h_sum_E_gamma_coarse_rest_frame[0] = bookHisto1D("Sum_E_Gamma_coarse_rest_frame_bare_e", 25, 0., 100.);
      _h_sum_E_gamma_coarse_rest_frame[1] = bookHisto1D("Sum_E_Gamma_coarse_rest_frame_dressed_1", 25, 0., 100.);
      _h_sum_E_gamma_coarse_rest_frame[2] = bookHisto1D("Sum_E_Gamma_coarse_rest_frame_dressed_2", 25, 0., 100.);
      _h_sum_pT_gamma_coarse_rest_frame[0] = bookHisto1D("Sum_pT_Gamma_coarse_rest_frame_bare_e", 25, 0., 100.);
      _h_sum_pT_gamma_coarse_rest_frame[1] = bookHisto1D("Sum_pT_Gamma_coarse_rest_frame_dressed_1", 25, 0., 100.);
      _h_sum_pT_gamma_coarse_rest_frame[2] = bookHisto1D("Sum_pT_Gamma_coarse_rest_frame_dressed_2", 25, 0., 100.);

      // Photon angle diagrams
      // h_theta_phot - angle to reference direction, all photons
      // h_theta_hardest_phot - angle to reference direction, hardest photon
      // h_theta_closest_phot - smallest photon lepton angle in units 2ml/mZ
      _h_theta_phot[0] = bookHisto1D("Theta_photon_bare_e", 50, 0., pi);
      _h_theta_phot[1] = bookHisto1D("Theta_photon_dressed_1", 50, 0., pi);
      _h_theta_phot[2] = bookHisto1D("Theta_photon_dressed_2", 50, 0., pi);
      _h_theta_closest_phot[0] = bookHisto1D("Theta_photon_closest_bare_e", 50, 0., 10.);
      _h_theta_closest_phot[1] = bookHisto1D("Theta_photon_closest_dressed_1", 50, 0., 10.);
      _h_theta_closest_phot[2] = bookHisto1D("Theta_photon_closest_dressed_2", 50, 0., 10.);
      _h_theta_hardest_phot[0] = bookHisto1D("Theta_hardest_photon_bare_e", 50, 0., pi);
      _h_theta_hardest_phot[1] = bookHisto1D("Theta_hardest_photon_dressed_1", 50, 0., pi);
      _h_theta_hardest_phot[2] = bookHisto1D("Theta_hardest_photon_dressed_2", 50, 0., pi);

      // invariant mass lepton antilepton
      _h_mll[0] = bookHisto1D("mee_bare", 100, 0., 150.);
      _h_mll[1] = bookHisto1D("mee_dressed_1", 100, 0., 150.);
      _h_mll[2] = bookHisto1D("mee_dressed_2", 100, 0., 150.);
      // invariant mass lepton antilepton photon - index 0 - hardest photon
      // index 1 = closest photon
      _h_mllg[0][0] = bookHisto1D("meeg_hardest_g_bare", 100, 0., 150.);
      _h_mllg[1][0] = bookHisto1D("meeg_hardest_g_dressed_1", 100, 0., 150.);
      _h_mllg[2][0] = bookHisto1D("meeg_hardest_g_dressed_2", 100, 0., 150.);
      _h_mllg[0][1] = bookHisto1D("meeg_closest_g_bare", 100, 0., 150.);
      _h_mllg[1][1] = bookHisto1D("meeg_closest_g_dressed_1", 100, 0., 150.);
      _h_mllg[2][1] = bookHisto1D("meeg_closest_g_dressed_2", 100, 0., 150.);
      // lepton pT and mT distributions (charged and all leptons)
      _h_pt_l[0] = bookHisto1D("pT_e_charged_bare", 100, 0., 75.);
      _h_pt_l[1] = bookHisto1D("pT_e_charged_dressed_1", 100, 0., 75.);
      _h_pt_l[2] = bookHisto1D("pT_e_charged_dressed_2", 100, 0., 75.);
      _h_pt_ln[0] = bookHisto1D("pT_e_all_bare", 100, 0., 75.);
      _h_pt_ln[1] = bookHisto1D("pT_e_all_dressed_1", 100, 0., 75.);
      _h_pt_ln[2] = bookHisto1D("pT_e_all_dressed_2", 100, 0., 75.);
      _h_mt_l[0] = bookHisto1D("mT_e_charged_bare", 100, 0., 150.);
      _h_mt_l[1] = bookHisto1D("mT_e_charged_dressed_1", 100, 0., 150.);
      _h_mt_l[2] = bookHisto1D("mT_e_charged_dressed_2", 100, 0., 150.);
      _h_mt_ln[0] = bookHisto1D("mT_e_all_bare", 100, 0., 150.);
      _h_mt_ln[1] = bookHisto1D("mT_e_all_dressed_1", 100, 0., 150.);
      _h_mt_ln[2] = bookHisto1D("mT_e_all_dressed_2", 100, 0., 150.);
      // mT lepton/antilepton system
      _h_mt_system_l[0] = bookHisto1D("mT_e_system_bare", 100, 0., 150.);
      _h_mt_system_l[1] = bookHisto1D("mT_e_system_dressed_1", 100, 0., 150.);
      _h_mt_system_l[2] = bookHisto1D("mT_e_system_dressed_2", 100, 0., 150.);

      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const Particles& Photons = 
	apply<IdentifiedFinalState>(event, "Photons").particles(Cuts::E > 0.1*GeV,cmpMomByE);
      const Particles& BareLeptons = apply<IdentifiedFinalState>(event, "BareLeptons").particles();
      const Particles& Dressed_1 = apply<DressedLeptons>(event, "Dressed01").particles();
      const vector<DressedLepton>& Dressed_1_dressed = apply<DressedLeptons>(event, "Dressed01").dressedLeptons();
      const Particles& Dressed_2 = apply<DressedLeptons>(event, "Dressed02").particles();
      const vector<DressedLepton>& Dressed_2_dressed = apply<DressedLeptons>(event, "Dressed02").dressedLeptons();
      const Particles& Neutrinos = apply<FinalState>(event, "Neutrinos").particles();

      // std::cout << "\n\nPhotons:" << std::endl;
      // for (size_t j = 0; j < Photons.size(); j++) {
      // 	std::cout << "j " << j << " mom " << Photons[j].momentum() << std::endl;
      // }
      // std::cout << "Bare leptons:" << std::endl;
      // for (size_t j = 0; j < BareLeptons.size(); j++) {
      // 	std::cout << "j " << j << " mom " << BareLeptons[j].momentum() << std::endl;
      // }
      // std::cout << "Dressed leptons 0.1:" << std::endl;
      // for (size_t j = 0; j < Dressed_1.size(); j++) {
      // 	std::cout << "j " << j << " pid " << Dressed_1[j].pid() << " mom " << Dressed_1[j].momentum() << std::endl;
      // }
      // for (size_t j = 0; j < Dressed_1_dressed.size(); j++) {
      // 	std::cout << "j " << j << " pid " << Dressed_1_dressed[j].pid() << " mom " << Dressed_1_dressed[j].momentum() << std::endl;
      // }
      // std::cout << "Dressed leptons 0.2:" << std::endl;
      // for (size_t j = 0; j < Dressed_2.size(); j++) {
      // 	std::cout << "j " << j << " pid " << Dressed_2[j].pid() << " mom " << Dressed_2[j].momentum() << std::endl;
      // }
      // for (size_t j = 0; j < Dressed_2_dressed.size(); j++) {
      // 	std::cout << "j " << j << " pid " << Dressed_2_dressed[j].pid() << " mom " << Dressed_2_dressed[j].momentum() << std::endl;
      // }

      for (int d = 0; d < 3; d++) {
	// No projection is filled - need to use electrons in runcard ...

	// Dressed leptons: do not use photons which were clustered?
	std::vector<DressedLepton> UseLeptons;
	Particles UsePhotons;
	Particles ClusteredPhotons;
	

	if (d == 0) {
	  for (size_t j = 0; j < BareLeptons.size(); j++) UseLeptons.push_back(BareLeptons[j]);
	  UsePhotons = Photons;
	}
	else {
	  if (d == 1) {
	    for (size_t j = 0; j < Dressed_1_dressed.size(); j++) {
	      if (Dressed_1_dressed[j].isLepton()) UseLeptons.push_back(Dressed_1_dressed[j]);
	    }
	    for (size_t j = 0; j < Dressed_1.size(); j++) {
	      if (Dressed_1[j].pid() == 22) ClusteredPhotons.push_back(Dressed_1[j]);
	    }
	  }
	  else if (d == 2) {
	    for (size_t j = 0; j < Dressed_2_dressed.size(); j++) {
	      if (Dressed_2_dressed[j].isLepton()) UseLeptons.push_back(Dressed_2_dressed[j]);
	    }
	    for (size_t j = 0; j < Dressed_2.size(); j++) {
	      if (Dressed_2[j].pid() == 22) ClusteredPhotons.push_back(Dressed_2[j]);
	    }
	  }
	  for (size_t j = 0; j < Photons.size(); j++) {
	    bool clustered=false;
	    for (size_t k = 0; k < ClusteredPhotons.size(); k++) {
	      if (Photons[j] == ClusteredPhotons[k]) clustered=true;
	    }
	    if (!clustered) UsePhotons.push_back(Photons[j]);
	  }
	}

	for (size_t j = 0; j < Neutrinos.size(); j++) UseLeptons.push_back(Neutrinos[j]);

	// if (d == 1 && UseLeptons.size() <= 1) {
	//   std::cout << "\n\nPhotons:" << std::endl;
	//   for (size_t j = 0; j < Photons.size(); j++) {
	//     std::cout << "j " << j << " mom " << Photons[j].momentum() << std::endl;
	//   }
	//   std::cout << "\nBareLeptons:" << std::endl;
	//   for (size_t j = 0; j < BareLeptons.size(); j++) {
	//     std::cout << "j " << j << " mom " << BareLeptons[j].momentum() << std::endl;
	//   }
	//   std::cout << "\nDressed leptons (particles) 0.1:" << std::endl;
	//   for (size_t j = 0; j < Dressed_1.size(); j++) {
	//     std::cout << "j " << j << " pid " << Dressed_1[j].pid() << " mom " << Dressed_1[j].momentum() << std::endl;
	//   }
	//   std::cout << "\nDressed leptons (dressed) 0.1:" << std::endl;
	//   for (size_t j = 0; j < Dressed_1_dressed.size(); j++) {
	//     std::cout << "j " << j << " pid " << Dressed_1_dressed[j].pid() << " mom " << Dressed_1_dressed[j].momentum() << std::endl;
	//   }
	//   std::cout << "\n\nUseLeptons:" << std::endl;
	//   for (size_t j = 0; j < UseLeptons.size(); j++) {
	//     std::cout << "j " << j << " mom " << UseLeptons[j].momentum() << std::endl;
	//   }
	//   std::cout << "\nUsePhotons:" << std::endl;
	//   for (size_t j = 0; j < UsePhotons.size(); j++) {
	//     std::cout << "j " << j << " mom " << UsePhotons[j].momentum() << std::endl;
	//   }
	// }
	// if (d == 2 && UseLeptons.size() <= 1) {
	//   std::cout << "\n\nPhotons:" << std::endl;
	//   for (size_t j = 0; j < Photons.size(); j++) {
	//     std::cout << "j " << j << " mom " << Photons[j].momentum() << std::endl;
	//   }
	//   std::cout << "\nDressed leptons (particles) 0.2:" << std::endl;
	//   for (size_t j = 0; j < Dressed_2.size(); j++) {
	//     std::cout << "j " << j << " pid " << Dressed_2[j].pid() << " mom " << Dressed_2[j].momentum() << std::endl;
	//   }
	//   std::cout << "\nDressed leptons (dressed) 0.2:" << std::endl;
	//   for (size_t j = 0; j < Dressed_2_dressed.size(); j++) {
	//     std::cout << "j " << j << " pid " << Dressed_2_dressed[j].pid() << " mom " << Dressed_2_dressed[j].momentum() << std::endl;
	//   }
	//   std::cout << "\n\nUseLeptons:" << std::endl;
	//   for (size_t j = 0; j < UseLeptons.size(); j++) {
	//     std::cout << "j " << j << " mom " << UseLeptons[j].momentum() << std::endl;
	//   }
	//   std::cout << "\nUsePhotons:" << std::endl;
	//   for (size_t j = 0; j < UsePhotons.size(); j++) {
	//     std::cout << "j " << j << " mom " << UsePhotons[j].momentum() << std::endl;
	//   }
	// }
      

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
	// if (UseLeptons.size() <= 1) {
	//   std::cout << "Only " << UseLeptons.size() << " lepton" << std::endl;
	//   std::cout << UseLeptons_COM[0] << std::endl;
	// }
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
	  _h_theta_phot[d]->fill(UseLeptons_COM[0].angle(phot_l_boost.vector3()),weight);
	  if (i == 0) _h_theta_hardest_phot[d]->fill(UseLeptons_COM[0].angle(phot_l_boost.vector3()),weight);
	  // Update smallest photon-lepton angle
	  if (i != 0) {
	    double angle = min(UseLeptons_COM[0].angle(phot_l_boost.vector3()),UseLeptons_COM[1].angle(phot_l_boost.vector3()));
	    if (angle < small_angle) small_angle = angle;
	  }
	}
	// Fill photon energy/pT distributions
	_h_sum_pT_gamma[d]->fill(sumPtgamma, weight);
	_h_sum_E_gamma[d]->fill(sumEgamma, weight);
	_h_sum_E_gamma_rest_frame[d]->fill(sumEgammaBoost, weight);
	_h_sum_pT_gamma_rest_frame[d]->fill(sumPtgammaBoost, weight);
	_h_sum_E_gamma_coarse_rest_frame[d]->fill(sumEgammaBoost, weight);
	_h_sum_pT_gamma_coarse_rest_frame[d]->fill(sumPtgammaBoost, weight);
	if (UseLeptons.size() > 0) {
	  double mZ = 91.1876; // check how to get from Rivet
	  double ml = UseLeptons[0].mass(); // Lepton[0] is guaranteed to be charged lepton
	  if (ml == 0.) { std::cout << "charged lepton massless" << std::endl; }
	  if (UsePhotons.size() !=0 && ml != 0.) _h_theta_closest_phot[d]->fill(small_angle*mZ/(2.*ml),weight);
	  for (size_t i = 0; i < UseLeptons.size(); ++i) {
	    if (UseLeptons[i].charge() != 0.) {
	      _h_pt_l[d]->fill(UseLeptons[i].momentum().pT(), weight);
	      _h_mt_l[d]->fill(mT(UseLeptons[i].momentum()), weight);
	    }
	    _h_pt_ln[d]->fill(UseLeptons[i].momentum().pT(), weight);
	    _h_mt_ln[d]->fill(mT(UseLeptons[i].momentum()), weight);
	  }
	  if (UseLeptons.size() > 1) {
	    _h_mt_system_l[d]->fill(mT(UseLeptons[0].momentum(),UseLeptons[1].momentum()),weight);
	    FourMomentum mom_ll = UseLeptons[0].momentum() + UseLeptons[1].momentum();
	    _h_mll[d]->fill(sqrt(mom_ll*mom_ll), weight);
	    _h_mllg[d][0]->fill(sqrt((mom_ll+Hardest_Photon)*(mom_ll+Hardest_Photon)),weight);
	    _h_mllg[d][1]->fill(sqrt((mom_ll+Closest_Photon)*(mom_ll+Closest_Photon)),weight);
	  }
	}
      }
      /// @todo Do the event by event analysis here
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double incxs = crossSection();
      double sumweight = sumOfWeights();

      for (int d = 0; d < 3; d++) {
	scale(_h_sum_E_gamma[d],1./sumweight);	       
	scale(_h_sum_pT_gamma[d],1./sumweight);	 
	scale(_h_sum_E_gamma_rest_frame[d],1./sumweight);
	scale(_h_sum_pT_gamma_rest_frame[d],1./sumweight);
	scale(_h_sum_E_gamma_coarse_rest_frame[d],1./sumweight);
	scale(_h_sum_pT_gamma_coarse_rest_frame[d],1./sumweight);
	scale(_h_mll[d],1./sumweight);	       
	scale(_h_mllg[d][0],1./sumweight);
	scale(_h_mllg[d][1],1./sumweight);
	scale(_h_pt_l[d],1./sumweight);
	scale(_h_pt_ln[d],1./sumweight);
	scale(_h_mt_l[d],1./sumweight);
	scale(_h_mt_ln[d],1./sumweight);
	scale(_h_mt_system_l[d],1./sumweight);
	scale(_h_theta_phot[d],1./sumweight);
	scale(_h_theta_closest_phot[d],1./sumweight);
	scale(_h_theta_hardest_phot[d],1./sumweight);
      }
      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity

    }      
      
    //@}


  private:

    
    /// @name Histograms
    //@{
    Histo1DPtr _h_sum_E_gamma[3];
    Histo1DPtr _h_sum_pT_gamma[3];
    Histo1DPtr _h_sum_E_gamma_rest_frame[3];
    Histo1DPtr _h_sum_pT_gamma_rest_frame[3];
    Histo1DPtr _h_sum_E_gamma_coarse_rest_frame[3];
    Histo1DPtr _h_sum_pT_gamma_coarse_rest_frame[3];
    Histo1DPtr _h_mll[3];
    Histo1DPtr _h_mllg[3][2];
    Histo1DPtr _h_pt_l[3];
    Histo1DPtr _h_pt_ln[3];
    Histo1DPtr _h_mt_l[3];
    Histo1DPtr _h_mt_ln[3];
    Histo1DPtr _h_mt_system_l[3];
    Histo1DPtr _h_theta_phot[3];
    Histo1DPtr _h_theta_closest_phot[3];
    Histo1DPtr _h_theta_hardest_phot[3];
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Zee);


}
