//  rivetinstall/rivetenv.sh
//  rivet-buildplugin RivetMC_WplusWplus.so MC_WplusWplus.cc
//  export RIVET_ANALYSIS_PATH=$(rivet-config --libdir):$PWD;  (also script)
//

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <map>

namespace Rivet {
  class MC_WplusWplus : public Analysis {
  private:
    double _pt_lepton1, _pt_lepton2, _eta_lepton, _ETmiss;
    double _pt_J, _eta_J, _Delta_R;
    double _Delta_Iso, _ptfrac_Iso, _pt_softlepton;
    double _m_JJ, _Delta_eta_JJ;
    double _pt_Jveto, _eta_offset_Jveto;
    double _m_ll_min, _m_ll_max;
    std::vector<double> _weights;
    
  public:
    MC_WplusWplus() : 
      Analysis("MC_WplusWplus"),
      _pt_lepton1(20.), _pt_lepton2(20.), _eta_lepton(2.5), _ETmiss(20.),
      _pt_J(30.), _eta_J(4.7), _Delta_R(0.4),
      _Delta_Iso(0.3), _ptfrac_Iso(0.15), _pt_softlepton(1.),
      _m_JJ(600.), _Delta_eta_JJ(4.), 
      _pt_Jveto(30.), _eta_offset_Jveto(1.),
      _m_ll_min(0.), _m_ll_max(1.e6)
    { }    

    void init() {
      // Set up projections
      // all fs particles in eta range
      FinalState incl_fs(-6.0, 6.0);
      addProjection(incl_fs, "FS");

      // for pTmiss
      // all visible fs particle in eta range
      VisibleFinalState vfs(-5.0,5.0);
      addProjection(vfs,"vfs");

      // all leptons in eta range with minimal pt
      // accept 11, -11, 13, -13
      IdentifiedFinalState fs_leptons(-_eta_lepton,_eta_lepton,0.*GeV);
      fs_leptons.acceptIdPair(11);
      fs_leptons.acceptIdPair(13);
      addProjection(fs_leptons, "Leptons");

      // jets
      FastJets fs_jets(incl_fs, FastJets::ANTIKT, _Delta_R);
      addProjection(fs_jets, "Jets");

      init_histos();

      // weights to count how many events survive successive cuts
      _weights.push_back(0.);
      _weights.push_back(0.);
      _weights.push_back(0.);
      _weights.push_back(0.);
      _weights.push_back(0.);
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& evt) {
      const double weight = evt.weight(); 
      _weights[0] += weight;

      const FinalState incl_fs = applyProjection<FinalState>(evt, "FS");
      const Jets          jets = 
	applyProjection<FastJets>(evt, "Jets").
	jetsByPt(_pt_J*GeV,14000.*GeV,-_eta_J,_eta_J); 

      // two leptons at least.
      const ParticleVector lepton_candidates = 
	applyProjection<FinalState>(evt, "Leptons").particles();
      if (lepton_candidates.size()<2) vetoEvent;

      // define missing ET and demand it to be above the cut
      ParticleVector vfs_particles = 
	applyProjection<VisibleFinalState>(evt, "vfs").particles();
      FourMomentum pTmiss;
      foreach ( const Particle & p, vfs_particles ) pTmiss -= p.momentum();
      double missET = pTmiss.pT();
      if (missET<_ETmiss) vetoEvent;

      // ensure leptons are in eta/pt range and sufficiently isolated,
      // no soft isolated lepton allowed.
      std::vector<FourMomentum> p_leptons;
      if (!ExtractLeptons(incl_fs,lepton_candidates,p_leptons)) vetoEvent;
      
      
      // filter leptons out of jets, must have at least two jets
      std::vector<FourMomentum> p_jets;
      if (!FilterJets(jets,p_leptons,p_jets)) vetoEvent;
      
      // search for two tagging jets.  Will have to check for proper cuts ....
      _weights[1] += weight; 
      double ptjp(0.), ptjm(0.),yjp(0.),yjm(0.);
      int jp(-1),jm(-1);
      const std::string mode(std::string("trivial"));
      if (!FindTwoTaggingJets(p_jets,jp,ptjp,yjp,jm,ptjm,yjm,mode)) vetoEvent;

      // apply jet veto
      _weights[2] += weight; 

      bool veto(false);
      double massjj((p_jets[jp]+p_jets[jm]).mass());
      double HT(p_jets[jp].pT()+p_jets[jm].pT());
      double yjj(abs(yjp-yjm)), ymean((yjp+yjm)/2.);
      _normed_histos["M_jetpm_bc"]->fill(massjj,weight);
      _normed_histos["dy_jetpm_bc"]->fill(yjj,weight);
      if (massjj<_m_JJ || yjj<_Delta_eta_JJ) veto = true;
      double yj(0), ptj, yjstar;
      for (int j=0;j<int(p_jets.size());j++) {
	if (j==jp || j==jm) continue;
	yj  = p_jets[j].rapidity();
	HT += ptj = p_jets[j].pT();
	if (ptj>_pt_Jveto &&
	    yj>yjm+_eta_offset_Jveto && 
	    yj<yjp-_eta_offset_Jveto) {
	  veto = true;
	  break;
	}
      }

      int n(0);
      for (int j=0;j<int(p_jets.size());j++) {
	yj     = p_jets[j].rapidity();
	yjstar = yj-ymean;
	ptj    = p_jets[j].pT();
	if (j==jp) {
	  _normed_histos["pT_jetp_bc"]->fill(ptj,weight);
	  _normed_histos["y_jetp_bc"]->fill(yj,weight);
	  _normed_histos["y*_jetp_bc"]->fill(yjstar,weight);
	  if (!veto) {
	    _normed_histos["pT_jetp"]->fill(ptj,weight);
	    _normed_histos["y_jetp"]->fill(yj,weight);
	    _normed_histos["y*_jetp"]->fill(yjstar,weight);
	  }
	}
	else if (j==jm) {
	  _normed_histos["pT_jetm_bc"]->fill(ptj,weight);
	  _normed_histos["y_jetm_bc"]->fill(yj,weight);
	  _normed_histos["y*_jetm_bc"]->fill(yjstar,weight);
	  if (!veto) {
	    _normed_histos["pT_jetm"]->fill(ptj,weight);
	    _normed_histos["y_jetm"]->fill(yj,weight);
	    _normed_histos["y*_jetm"]->fill(yjstar,weight);
	  }
	}
	else {
	  n++;
	  switch (n) {
	  case 1: 
	    _normed_histos["pT_jet1_bc"]->fill(ptj,weight);
	    _normed_histos["y_jet1_bc"]->fill(yj,weight);
	    _normed_histos["y*_jet1_bc"]->fill(yjstar,weight);
	    if (!veto) {
	      _normed_histos["pT_jet1"]->fill(ptj,weight);
	      _normed_histos["y_jet1"]->fill(yj,weight);
	      _normed_histos["y*_jet1"]->fill(yjstar,weight);
	    }
	    break;
	  case 2: 
	    _normed_histos["pT_jet2_bc"]->fill(ptj,weight);
	    _normed_histos["y_jet2_bc"]->fill(yj,weight);
	    _normed_histos["y*_jet2_bc"]->fill(yjstar,weight);
	    if (!veto) {
	      _normed_histos["pT_jet2"]->fill(ptj,weight);
	      _normed_histos["y_jet1"]->fill(yj,weight);
	      _normed_histos["y*_jet2"]->fill(yjstar,weight);
	    }
	    break;
	  }
	}
      }
      if (veto) vetoEvent;
      _weights[3] += weight;
      _normed_histos["M_jetpm"]->fill(massjj,weight);
      _normed_histos["dy_jetpm"]->fill(yjj,weight);
      _normed_histos["HT_jets"]->fill(HT,weight);
      _normed_histos["missET"]->fill(missET,weight);
      HT += p_leptons[0].pT()+p_leptons[1].pT()+missET;
      _normed_histos["HT"]->fill(HT,weight);
      double phijj = abs(p_jets[jp].phi()-p_jets[jm].phi());
      if (phijj>PI) phijj = abs(phijj-PI);
      _normed_histos["dphi_jetpm"]->fill(phijj,weight);
    }
  private:
    std::map<std::string,AIDA::IHistogram1D*> _normed_histos;
    
    void init_histos() { 
      _normed_histos["pT_jetp"]     = bookHistogram1D("pT_jetp",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["pT_jetm"]     = bookHistogram1D("pT_jetm",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["pT_jet1"]     = bookHistogram1D("pT_jet1",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["pT_jet2"]     = bookHistogram1D("pT_jet2",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["y_jetp"]      = bookHistogram1D("y_jetp",20,-5.,5.0);
      _normed_histos["y_jetm"]      = bookHistogram1D("y_jetm",20,-5.,5.0);
      _normed_histos["y_jet1"]      = bookHistogram1D("y_jet1",20,-5.,5.0);
      _normed_histos["y_jet2"]      = bookHistogram1D("y_jet2",20,-5.,5.0);
      _normed_histos["y*_jetp"]     = bookHistogram1D("y*_jetp",20,-5.,5.0);
      _normed_histos["y*_jetm"]     = bookHistogram1D("y*_jetm",20,-5.,5.0);
      _normed_histos["y*_jet1"]     = bookHistogram1D("y*_jet1",20,-5.,5.0);
      _normed_histos["y*_jet2"]     = bookHistogram1D("y*_jet2",20,-5.,5.0);
      _normed_histos["HT_jets"]     = bookHistogram1D("HT_jets",100,100.,1100.);
      _normed_histos["HT"]          = bookHistogram1D("HT",100,100.,1100.);
      _normed_histos["missET"]      = bookHistogram1D("missET",100,_ETmiss,
						      _ETmiss+100+10.);
      _normed_histos["M_jetpm"]     = bookHistogram1D("M_jetpm",50,2.*_pt_J,
						      1000.0+2.*_pt_J);
      _normed_histos["dy_jetpm"]    = bookHistogram1D("dy_jetpm",20,0.,10.0);
      _normed_histos["dphi_jetpm"]  = bookHistogram1D("dphi_jetpm",63,0.,PI);
      
      
      _normed_histos["pT_jetp_bc"]  = bookHistogram1D("pT_jetp_bc",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["pT_jetm_bc"]  = bookHistogram1D("pT_jetm_bc",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["pT_jet1_bc"]  = bookHistogram1D("pT_jet1_bc",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["pT_jet2_bc"]  = bookHistogram1D("pT_jet2_bc",100,_pt_J,
						      5.*100.+_pt_J);
      _normed_histos["y_jetp_bc"]   = bookHistogram1D("y_jetp_bc",20,-5.,5.0);
      _normed_histos["y_jetm_bc"]   = bookHistogram1D("y_jetm_bc",20,-5.,5.0);
      _normed_histos["y_jet1_bc"]   = bookHistogram1D("y_jet1_bc",20,-5.,5.0);
      _normed_histos["y_jet2_bc"]   = bookHistogram1D("y_jet2_bc",20,-5.,5.0);
      _normed_histos["y*_jetp_bc"]  = bookHistogram1D("y*_jetp_bc",20,-5.,5.0);
      _normed_histos["y*_jetm_bc"]  = bookHistogram1D("y*_jetm_bc",20,-5.,5.0);
      _normed_histos["y*_jet1_bc"]  = bookHistogram1D("y*_jet1_bc",20,-5.,5.0);
      _normed_histos["y*_jet2_bc"]  = bookHistogram1D("y*_jet2_bc",20,-5.,5.0);
      _normed_histos["M_jetpm_bc"]  = bookHistogram1D("M_jetpm_bc",50,2.*_pt_J,
						      1000.0+2.*_pt_J);
      _normed_histos["dy_jetpm_bc"] = bookHistogram1D("dy_jetpm_bc",20,0.,10.0);
    }  

    bool ExtractLeptons(const FinalState& incl_fs,
			const ParticleVector &lepton_candidates,
			std::vector<FourMomentum> &p_leptons) {
      std::list<FourMomentum> p_softleptons;
      foreach (const Particle& cand,lepton_candidates) {
	FourMomentum candmom = cand.momentum();
	double cand_pT = candmom.pT(), cone_pT = 0.0;
	if (candmom.eta()>_eta_lepton ||
	    candmom.eta()<-_eta_lepton) continue;
	foreach (const Particle& p, incl_fs.particles()) {
	  if (deltaR(candmom, p.momentum()) < _Delta_Iso &&
	      candmom.E() != p.momentum().E()) {
	    cone_pT += p.momentum().pT();
	  }
	}
	if (cone_pT/cand_pT < _ptfrac_Iso) {
	  if (cand_pT>_pt_lepton2) 
	    p_leptons.push_back(candmom);
	  else if (cand_pT>_pt_softlepton)
	    p_softleptons.push_back(candmom);
	}
      }
      if (p_leptons.size()!=2 || p_softleptons.size()!=0) return false;
      if (p_leptons[0].pT()<p_leptons[1].pT()) {
	FourMomentum swap0(p_leptons[1]), swap1(p_leptons[0]);
	p_leptons.clear();
	p_leptons.push_back(swap0);
	p_leptons.push_back(swap1);
      }
      if (p_leptons[0].pT()<_pt_lepton1 ||
	  p_leptons[1].pT()<_pt_lepton2)  return false;
      return true;
    }
    bool FilterJets(const Jets& jets,
		    std::vector<FourMomentum> &p_leptons,
		    std::vector<FourMomentum> &p_jets)
    { 
      FourMomentum p_J;
      foreach (const Jet& jet,jets) {
	p_J = jet.momentum();
	bool vetoJ(false);
	foreach (FourMomentum & mom,p_leptons) {
	  if (deltaR(mom,p_J)<_Delta_R) {
	    vetoJ = true;
	    break;
	  }
	}
	if (vetoJ) continue;
	p_jets.push_back(p_J);
      }
      if (p_jets.size()<2) return false;
      return true;
    }
    bool FindTwoTaggingJets(std::vector<FourMomentum> &p_jets,
			    int & jp,double & ptjp,double & yjp,
			    int & jm,double & ptjm,double & yjm,
			    const std::string & mode)
    {
      jp = jm = -1;
      double ptj, yj;
      FourMomentum pjet;
      for (size_t i(0);i<p_jets.size();i++) {
	pjet = p_jets[i];
	yj   = pjet.rapidity();
	ptj  = pjet.pT();
	if (mode==std::string("trivial")) {
	  if (yj>0. && ptj>ptjp) {
	    jp = i; ptjp = ptj; yjp = yj;
	  }
	  else if (yj<0. && ptj>ptjm) {
	    jm = i; ptjm = ptj; yjm = yj;
	  }
	}
      }
      if (jp<0||jm<0) return false;
      return true;
    }
  public:    
    /// Normalise histograms etc., after the run
    void finalize() {
      for (std::map<std::string,AIDA::IHistogram1D*>::iterator 
	     hit=_normed_histos.begin();
	   hit!=_normed_histos.end();hit++) {
	scale(hit->second, crossSection()/sumOfWeights());
      }
    }
  };

  // The hook for the plugin system
  AnalysisBuilder<MC_WplusWplus> plugin_MC_WplusWplus;
}
