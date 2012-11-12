#include "Rivet/Analysis.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  class MC_TTW : public Analysis {
  private:
    double _leta, _mupt, _ept, _lpt, _missET;
    double _isoE_lj, _isoR_lj, _jeta, _jR, _jpt; 

    AIDA::IHistogram1D *_h_njets; //HT?
    AIDA::IHistogram1D *_h_mass_ll_ss,*_h_mass_ll_os,*_h_mass_lll;
    AIDA::IHistogram1D *_h_pT_l1,*_h_pT_l2,*_h_pT_l3,*_h_missET;
    AIDA::IHistogram1D *_h_eta_l1,*_h_eta_l2,*_h_eta_l3; 
    AIDA::IHistogram1D *_h_pT_t1,*_h_pT_t2,*_h_pT_W,*_h_pT_tt,*_h_pT_ttW;
    AIDA::IHistogram1D *_h_mass_tt,*_h_mass_ttW;
    AIDA::IHistogram1D *_h_rap_t1,*_h_rap_t2,*_h_rap_W,*_h_rap_tt,*_h_rap_ttW;
    AIDA::IHistogram1D *_h_pT_b1,*_h_pT_b2,*_h_pT_lj;
    AIDA::IHistogram1D *_h_eta_b1,*_h_eta_b2,*_h_eta_lj; 

    void inithistos() {
      _h_njets      = bookHistogram1D("jet_mult", 11, -0.5, 10.5);
      _h_mass_ll_ss = bookHistogram1D("mass_ll_ss", 50, 0.0, 150.0);
      _h_mass_ll_os = bookHistogram1D("mass_ll_os", 50, 0.0, 150.0);
      _h_mass_lll   = bookHistogram1D("mass_lll", 50, 0.0, 150.0);
      _h_pT_l1      = bookHistogram1D("pT_l1", logspace(10.0, 500.0, 50));
      _h_pT_l2      = bookHistogram1D("pT_l2", logspace(10.0, 500.0, 50));
      _h_pT_l3      = bookHistogram1D("pT_l3", logspace(10.0, 500.0, 50));
      _h_eta_l1     = bookHistogram1D("eta_l1", 10, -_leta, _leta);
      _h_eta_l2     = bookHistogram1D("eta_l2", 10, -_leta, _leta); 
      _h_eta_l3     = bookHistogram1D("eta_l3", 10, -_leta, _leta); 
      _h_missET     = bookHistogram1D("missET", logspace(10.0, 500.0, 50));
      _h_pT_t1      = bookHistogram1D("pT_t1",logspace(1.0, 500.0, 50));
      _h_pT_t2      = bookHistogram1D("pT_t2", logspace(1.0, 500.0, 50));
      _h_pT_W       = bookHistogram1D("pT_W", logspace(1.0, 500.0, 50));
      _h_pT_tt      = bookHistogram1D("pT_tt", logspace(1.0, 500.0, 50));
      _h_pT_ttW     = bookHistogram1D("pT_ttW", logspace(1.0, 500.0, 50));
      _h_mass_tt    = bookHistogram1D("mass_tt", 50, 350.0, 1350.0);
      _h_mass_ttW   = bookHistogram1D("mass_ttW", 50, 450.0, 1450.0);
      _h_rap_t1     = bookHistogram1D("rap_t1", 20, -5., 5.);
      _h_rap_t2     = bookHistogram1D("rap_t2", 20, -5., 5.);
      _h_rap_W      = bookHistogram1D("rap_W", 20, -5., 5.); 
      _h_rap_tt     = bookHistogram1D("rap_tt", 20, -5., 5.);
      _h_rap_ttW    = bookHistogram1D("rap_ttW", 20, -5., 5.);
      _h_pT_b1      = bookHistogram1D("pT_b1", logspace(10.0, 500.0, 50));
      _h_pT_b2      = bookHistogram1D("pT_b2", logspace(10.0, 500.0, 50));
      _h_pT_lj      = bookHistogram1D("pT_b2", logspace(10.0, 500.0, 50));
      _h_eta_b1     = bookHistogram1D("eta_b1", 20, -_jeta, _jeta);
      _h_eta_b2     = bookHistogram1D("eta_b2", 20, -_jeta, _jeta); 
      _h_eta_lj     = bookHistogram1D("eta_b2", 20, -_jeta, _jeta); 
    }
  public:
    MC_TTW() : 
      Analysis("MC_TTW"),
      _leta(2.5), _mupt(10.), _ept(10.), _lpt(min(_ept,_mupt)),
      _missET(50.), _isoE_lj(0.05), _isoR_lj(0.2), 
      _jeta(5.), _jR(0.4), _jpt(20.) 
    {}
    
    void init() {
      ChargedLeptons lfs(FinalState(-_leta, _leta, _lpt*GeV));
      addProjection(lfs, "LFS");
      VisibleFinalState vfs(VisibleFinalState(-_jeta, _jeta));
      addProjection(vfs, "VFS");
      VetoedFinalState fs(FinalState(-_jeta, _jeta, 0*GeV));
      fs.addVetoOnThisFinalState(lfs);
      addProjection(FastJets(fs, FastJets::ANTIKT, _jR), "Jets");
      addProjection(MissingMomentum(fs), "MissingET");
      IdentifiedFinalState nfs(-50., 50., 0.0*GeV);
      nfs.acceptNeutrinos();
      addProjection(nfs, "NFS");
      
      inithistos();
    }
    void analyze(const Event& event) {
      const double weight = event.weight();
      // Use the "LFS" projection to require at least one hard charged
      // lepton. This is an experimental signature for the leptonically 
      // decaying W. This helps to reduce pure QCD backgrounds.
      const ParticleVector& leptons = 
	applyProjection<ChargedLeptons>(event, "LFS").particlesByPt();
      const ParticleVector& tracks = 
	applyProjection<VisibleFinalState>(event, "VFS").particles();
      const ParticleVector& neutrinos = 
	applyProjection<IdentifiedFinalState>(event, "NFS").particlesByPt();
      const MissingMomentum& met = 
	applyProjection<MissingMomentum>(event, "MissingET");
      const Jets alljets = 
	applyProjection<FastJets>(event, "Jets").jetsByPt();

      double misset = met.vectorEt().mod();
      if (leptons.size()<3 || neutrinos.size()<3 ||
	  misset < _missET*GeV) {
	vetoEvent;
      }
      
      ParticleVector isolatedLeptons;
      foreach (const Particle& lepton, leptons) {
	double coneET(0.);
	foreach (const Particle& track, tracks) {
	  if (deltaR(track.momentum(), lepton.momentum()) < _isoR_lj) {
	    coneET += track.momentum().pT();
	  }
	}
	if (coneET < (1.+_isoE_lj)*lepton.momentum().pT())
	  isolatedLeptons.push_back(lepton);
      }
      if (isolatedLeptons.size()!=3) vetoEvent;

      Jets bjets, ljets;
      foreach (const Jet& jet, alljets) {
        if (jet.containsBottom()) bjets.push_back(jet);
        else ljets.push_back(jet);
      }
      if (bjets.size()<2) vetoEvent;

      // reconstruct 3 W's.
      FourMomentum W1, W2, W3;
      FourMomentum lepmoms[3], neumoms[3];
      for (size_t i=0;i<3;i++) {
	lepmoms[i] = leptons[i].momentum();
	neumoms[i] = neutrinos[i].momentum();
      }
      std::set<int> lepIDs, neuIDs;
      lepIDs.insert(0);lepIDs.insert(1);lepIDs.insert(2);
      neuIDs.insert(0);neuIDs.insert(1);neuIDs.insert(2);
    }


    void finalize() {
      normalize(_h_njets); 
      normalize(_h_mass_ll_ss); 
      normalize(_h_mass_ll_os); 
      normalize(_h_mass_lll); 
      normalize(_h_pT_l1); 
      normalize(_h_pT_l2); 
      normalize(_h_pT_l3); 
      normalize(_h_eta_l1); 
      normalize(_h_eta_l2); 
      normalize(_h_eta_l3); 
      normalize(_h_missET); 
      normalize(_h_pT_t1); 
      normalize(_h_pT_t2); 
      normalize(_h_pT_W); 
      normalize(_h_pT_tt); 
      normalize(_h_pT_ttW); 
      normalize(_h_mass_tt); 
      normalize(_h_mass_ttW); 
      normalize(_h_rap_t1); 
      normalize(_h_rap_t2); 
      normalize(_h_rap_W); 
      normalize(_h_rap_tt); 
      normalize(_h_rap_ttW); 
      normalize(_h_pT_b1); 
      normalize(_h_pT_b2); 
      normalize(_h_pT_lj); 
      normalize(_h_eta_b1); 
      normalize(_h_eta_b2); 
      normalize(_h_eta_lj); 
    }
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTW);

}
