#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  class MC_TTSpins : public Analysis {
  private:
    double _leta, _lpt, _missET, _iso_lj; 

    AIDA::IHistogram1D *_h_njets;
    AIDA::IHistogram1D *_h_jet_1_pT,*_h_jet_2_pT,*_h_jet_3_pT,*_h_jet_4_pT;
    AIDA::IHistogram1D *_h_jet_HT;
    AIDA::IHistogram1D *_h_bjet_1_pT,*_h_bjet_2_pT,*_h_ljet_1_pT,*_h_ljet_2_pT;
    AIDA::IHistogram1D *_h_bb_dR, *_h_bb_deta, *_h_bb_dphi;
    AIDA::IHistogram1D *_h_ll_dR, *_h_ll_deta, *_h_ll_dphi;
    AIDA::IHistogram1D *_h_bl_dR, *_h_bl_deta, *_h_bl_dphi;

    void inithistos() {
      _h_njets    = bookHistogram1D("jet_mult", 11, -0.5, 10.5);
      _h_jet_1_pT = bookHistogram1D("jet_1_pT", logspace(20.0, 500.0, 50));
      _h_jet_2_pT = bookHistogram1D("jet_2_pT", logspace(20.0, 400.0, 50));
      _h_jet_3_pT = bookHistogram1D("jet_3_pT", logspace(20.0, 300.0, 50));
      _h_jet_4_pT = bookHistogram1D("jet_4_pT", logspace(20.0, 200.0, 50));
      _h_jet_HT   = bookHistogram1D("jet_HT", logspace(100.0, 2000.0, 50));
      //
      _h_bjet_1_pT = bookHistogram1D("jetb_1_pT", logspace(20.0, 400.0, 50));
      _h_bjet_2_pT = bookHistogram1D("jetb_2_pT", logspace(20.0, 300.0, 50));
      _h_ljet_1_pT = bookHistogram1D("jetl_1_pT", logspace(20.0, 400.0, 50));
      _h_ljet_2_pT = bookHistogram1D("jetl_2_pT", logspace(20.0, 300.0, 50));
      //
      _h_bb_dR    = bookHistogram1D("bb_dR",28, 0.0, 7.0);
      _h_bb_deta  = bookHistogram1D("bb_deta",28, 0.0, 7.0);
      _h_bb_dphi  = bookHistogram1D("bb_dphi",32, 0.0, 6.2);
      _h_ll_dR    = bookHistogram1D("ll_dR",28, 0.0, 7.0);
      _h_ll_deta  = bookHistogram1D("ll_deta",28, 0.0, 7.0);
      _h_ll_dphi  = bookHistogram1D("ll_dphi",32, 0.0, M_PI);
      _h_bl_dR    = bookHistogram1D("bl_dR",28, 0.0, 7.0);
      _h_bl_deta  = bookHistogram1D("bl_deta",28, 0.0, 7.0);
      _h_bl_dphi  = bookHistogram1D("bl_dphi",32, 0.0, 6.2);
      std::cerr<<"out\n";
    }
  public:
    MC_TTSpins() : 
      Analysis("MC_TTSpins"),
      _leta(4.2), _lpt(30.), _missET(30.), _iso_lj(0.3)
    {}
    
    void init() {
      ChargedLeptons lfs(FinalState(-_leta, _leta, _lpt*GeV));
      addProjection(lfs, "LFS");
      VetoedFinalState fs(FinalState(-_leta, _leta, 0*GeV));
      fs.addVetoOnThisFinalState(lfs);
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.6), "Jets");
      addProjection(MissingMomentum(fs), "MissingET");
      
      inithistos();
    }
    void analyze(const Event& event) {
      const double weight = event.weight();
      // Use the "LFS" projection to require at least one hard charged
      // lepton. This is an experimental signature for the leptonically 
      // decaying W. This helps to reduce pure QCD backgrounds.
      const ChargedLeptons& lfs = 
	applyProjection<ChargedLeptons>(event, "LFS");
      MSG_DEBUG("Charged lepton multiplicity = " << 
		lfs.chargedLeptons().size());
      foreach (const Particle& lepton, lfs.chargedLeptons()) {
        MSG_DEBUG("Lepton pT = " << lepton.momentum().pT());
      }
      if (lfs.chargedLeptons().empty()) {
        MSG_DEBUG("Event failed lepton multiplicity cut");
        vetoEvent;
      }
      
      // Use a missing ET cut to bias toward events with a hard neutrino 
      // from the leptonically decaying W. This helps to reduce pure QCD 
      // backgrounds.
      const MissingMomentum& met = 
	applyProjection<MissingMomentum>(event, "MissingET");
      MSG_DEBUG("Vector ET = " << met.vectorEt().mod() << " GeV");
      if (met.vectorEt().mod() < _missET*GeV) {
        MSG_DEBUG("Event failed missing ET cut");
        vetoEvent;
      }
      
      // Use the "Jets" projection to check that there are at least 4 
      // jets of any pT.  Getting the jets sorted by pT ensures that the 
      // first jet is the hardest, and so on. We apply no pT cut here only 
      // because we want to plot all jet pTs to help optimise our jet pT cut.
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets alljets = jetpro.jetsByPt();
      //if (alljets.size() < 2) {
      //  MSG_DEBUG("Event failed jet multiplicity cut");
      //  vetoEvent;
      //}
      
      // Update passed-cuts counter and fill all-jets histograms
      int njets(alljets.size());
      if (njets>0)
	_h_jet_1_pT->fill(alljets[0].momentum().pT()/GeV, weight);
      if (njets>1)
	_h_jet_2_pT->fill(alljets[1].momentum().pT()/GeV, weight);
      if (njets>2)
	_h_jet_3_pT->fill(alljets[2].momentum().pT()/GeV, weight);
      if (njets>3)
	_h_jet_4_pT->fill(alljets[3].momentum().pT()/GeV, weight);
      
      
      // Sort the jets into b-jets and light jets. We expect one hard b-jet 
      // from each top decay, so our 4 hardest jets should include two b-jets. 
      // The Jet::containsBottom() method is equivalent to perfect experimental
      // b-tagging, in a generator-independent way.
      Jets bjets, ljets;
      foreach (const Jet& jet, alljets) {
        // Veto events with jets that overlap with the hard leptons
        bool isolated = true;
        foreach (const Particle& lepton, lfs.chargedLeptons()) {
          if (deltaR(jet.momentum(), lepton.momentum()) < _iso_lj) {
            isolated = false;
            break;
          }
        }
        if (!isolated) {
          MSG_DEBUG("Jet failed lepton isolation cut");
	  vetoEvent;
        }
        if (jet.containsBottom()) bjets.push_back(jet);
        else ljets.push_back(jet);
      }
      MSG_DEBUG("Number of b-jets = " << bjets.size());
      MSG_DEBUG("Number of l-jets = " << ljets.size());

      // Plot the pTs of the identified jets.
      int nbjets(bjets.size()), nljets(ljets.size());
      if (nbjets>0) 
	_h_bjet_1_pT->fill(bjets[0].momentum().pT(), weight);
      if (nbjets>1) 
	_h_bjet_2_pT->fill(bjets[1].momentum().pT(), weight);
      if (nljets>0) 
	_h_ljet_1_pT->fill(ljets[0].momentum().pT(), weight);
      if (nljets>1) 
	_h_ljet_2_pT->fill(ljets[1].momentum().pT(), weight);

      if (nbjets==2) {
	_h_bb_dR->fill(deltaR(bjets[0].momentum(), 
			      bjets[1].momentum()),weight);
        _h_bb_deta->fill(fabs(bjets[0].momentum().eta()-
			      bjets[1].momentum().eta()),weight);
        _h_bb_dphi->fill(deltaPhi(bjets[0].momentum(),
				  bjets[1].momentum()),weight);
      }


      if (lfs.chargedLeptons().size()>0 && nbjets>0) {
	FourMomentum l=lfs.chargedLeptons()[0].momentum();
	_h_bl_dR->fill(deltaR(bjets[0].momentum(), l),weight);
	_h_bl_deta->fill(fabs(bjets[0].momentum().eta()-l.eta()),weight);
	_h_bl_dphi->fill(deltaPhi(bjets[0].momentum(),l),weight);
      }
      if (lfs.chargedLeptons().size()>1 && nbjets>0) {
	FourMomentum l1=lfs.chargedLeptons()[0].momentum();
	FourMomentum l2=lfs.chargedLeptons()[1].momentum();
	_h_ll_dR->fill(deltaR(l1, l2),weight);
	_h_ll_deta->fill(fabs(l1.eta()-l2.eta()),weight);
	_h_ll_dphi->fill(deltaPhi(l1,l2),weight);
      }
    }


    void finalize() {
      normalize(_h_njets);
      normalize(_h_jet_1_pT);
      normalize(_h_jet_2_pT);
      normalize(_h_jet_3_pT);
      normalize(_h_jet_4_pT);
      normalize(_h_jet_HT);
      normalize(_h_bjet_1_pT);
      normalize(_h_bjet_2_pT);
      normalize(_h_ljet_1_pT);
      normalize(_h_ljet_2_pT);
      normalize(_h_bb_dR);
      normalize(_h_bb_deta);
      normalize(_h_bb_dphi);
      normalize(_h_bl_dR);
      normalize(_h_bl_deta);
      normalize(_h_bl_dphi);
      normalize(_h_ll_dR);
      normalize(_h_ll_deta);
      normalize(_h_ll_dphi);
    }
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTSpins);

}
