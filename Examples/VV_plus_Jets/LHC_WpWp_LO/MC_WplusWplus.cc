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
    double _pt_Jlow, _pt_Jhigh, _eta_J, _Delta_R;
    double _Delta_Iso, _ptfrac_Iso, _pt_softlepton;
    double _m_JJ, _Delta_eta_JJ;
    double _pt_Jveto, _eta_offset_Jveto;
    double _m_ll_min, _m_ll_max;
    std::vector<double> _weights;
    
  public:
    MC_WplusWplus() : 
      Analysis("MC_WplusWplus"),
      _pt_lepton1(20.), _pt_lepton2(20.), _eta_lepton(2.5), _ETmiss(20.),
      _pt_Jlow(25.), _pt_Jhigh(40.), _eta_J(4.7), _Delta_R(0.4),
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
      _weights.push_back(0.);
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& evt) {
      const double weight = evt.weight(); 
      _weights[0] += weight;
      _histos["weights"]->fill(0,weight);

      const FinalState incl_fs = applyProjection<FinalState>(evt, "FS");
      // two leptons at least.
      const ParticleVector lepton_candidates = 
	applyProjection<FinalState>(evt, "Leptons").particles();
      if (lepton_candidates.size()<2) vetoEvent;
      _weights[1] += weight;
      _histos["weights"]->fill(1,weight);

      // ensure leptons are in eta/pt range and sufficiently isolated,
      // no soft isolated lepton allowed.
      std::vector<FourMomentum> p_leptons;
      if (!ExtractLeptons(incl_fs,lepton_candidates,p_leptons)) vetoEvent;
      _weights[2] += weight;
      _histos["weights"]->fill(2,weight);

      // define missing ET and demand it to be above the cut
      ParticleVector vfs_particles = 
	applyProjection<VisibleFinalState>(evt, "vfs").particles();
      FourMomentum pTmiss;
      foreach ( const Particle & p, vfs_particles ) pTmiss -= p.momentum();
      double missET = pTmiss.pT();
      if (missET<_ETmiss) vetoEvent;
      _weights[3] += weight;
      _histos["weights"]->fill(3,weight);

      const Jets      jets_low = 
	applyProjection<FastJets>(evt, "Jets").
	jetsByPt(_pt_Jlow*GeV,8000.*GeV,-_eta_J,_eta_J);       
      const Jets     jets_high = 
	applyProjection<FastJets>(evt, "Jets").
	jetsByPt(_pt_Jlow*GeV,8000.*GeV,-_eta_J,_eta_J); 

      // filter leptons out of jets
      std::vector<FourMomentum> p_jets_low;
      int njets_low = FilterJets(jets_low,p_leptons,p_jets_low);
      std::vector<FourMomentum> p_jets;
      int njets_high = FilterJets(jets_high,p_leptons,p_jets);
      _histos["Njets_low"]->fill(njets_low,weight);
      _histos["Njets_high"]->fill(njets_high,weight);
      
      _histos["missET_all"]->fill(missET,weight);
      _histos["ptL1_all"]->fill(p_leptons[0].pT(),weight);
      _histos["ptL2_all"]->fill(p_leptons[1].pT(),weight);
      _histos["dphiLL_all"]->fill(deltaPhi(p_leptons[0],
					   p_leptons[1]),weight);
      _histos["detaLL_all"]->fill(fabs(p_leptons[0].eta()-
				       p_leptons[1].eta()),weight);

      if (njets_high==0) {
	_histos["missET_0_high"]->fill(missET,weight);
	_histos["ptL1_0_high"]->fill(p_leptons[0].pT(),weight);
	_histos["ptL2_0_high"]->fill(p_leptons[1].pT(),weight);
	_histos["dphiLL_0_high"]->fill(deltaPhi(p_leptons[0],
						       p_leptons[1]),weight);
	_histos["detaLL_0_high"]->fill(fabs(p_leptons[0].eta()-
						   p_leptons[1].eta()),weight);
	_weights[4] += weight; 
	_histos["weights"]->fill(4,weight);
      }
      if (njets_low==0) {
	_histos["missET_0_low"]->fill(missET,weight);
	_histos["ptL1_0_low"]->fill(p_leptons[0].pT(),weight);
	_histos["ptL2_0_low"]->fill(p_leptons[1].pT(),weight);
	_histos["dphiLL_0_low"]->fill(deltaPhi(p_leptons[0],
						      p_leptons[1]),weight);
	_histos["detaLL_0_low"]->fill(fabs(p_leptons[0].eta()-
						  p_leptons[1].eta()),weight);
	_weights[5] += weight; 
	_histos["weights"]->fill(5,weight);
      }
      else if (njets_low==1) {
	_histos["missET_1_low"]->fill(missET,weight);
	_histos["ptL1_1_low"]->fill(p_leptons[0].pT(),weight);
	_histos["ptL2_1_low"]->fill(p_leptons[1].pT(),weight);
	_histos["dphiLL_1_low"]->fill(deltaPhi(p_leptons[0],
						      p_leptons[1]),weight);
	_histos["detaLL_1_low"]->fill(fabs(p_leptons[0].eta()-
						  p_leptons[1].eta()),weight);
      }
	

      // Now go for VBF-like signatures
      // the two tagging jets.  Will have to check for proper cuts ....
      double ptjp(0.), ptjm(0.),yjp(0.),yjm(0.);
      int jp(-1),jm(-1);
      const std::string mode(std::string("trivial"));
      if (!FindTwoTaggingJets(p_jets,jp,ptjp,yjp,jm,ptjm,yjm,mode)) vetoEvent;
      _weights[6] += weight; 
      _histos["weights"]->fill(6,weight);

      // apply mass and jet vetoes 
      double massjj((p_jets[jp]+p_jets[jm]).mass());
      double HT(p_jets[jp].pT()+p_jets[jm].pT());
      double yjj(abs(yjp-yjm)), ymean((yjp+yjm)/2.);
      _histos["M_jetpm"]->fill(massjj,weight);
      _histos["dy_jetpm"]->fill(yjj,weight);

      if (massjj<_m_JJ || yjj<_Delta_eta_JJ)  vetoEvent;
      _weights[7] += weight; 
      _histos["weights"]->fill(7,weight);

      bool veto(false);
      double yj(0), ptj, yjstar;
      for (int j=0;j<int(p_jets.size());j++) {
	if (j==jp || j==jm) continue;
	yj  = p_jets[j].rapidity();
	HT += ptj = p_jets[j].pT();
	if (ptj>_pt_Jveto &&
	    yj>yjm+_eta_offset_Jveto && 
	    yj<yjp-_eta_offset_Jveto) {
	  _histos["pT_jet_veto"]->fill(ptj,weight);
	  _histos["y_jet_veto"]->fill(fabs(yj),weight);
	  _histos["y*_jet_veto"]->fill(fabs(yj-ymean),weight);
	  veto = true;
	}
      }
      if (veto) vetoEvent;
      _weights[8] += weight; 
      _histos["weights"]->fill(8,weight);

      int n(0);
      for (int j=0;j<int(p_jets.size());j++) {
	yj     = p_jets[j].rapidity();
	yjstar = yj-ymean;
	ptj    = p_jets[j].pT();
	if (j==jp) {
	  _histos["pT_jetp"]->fill(ptj,weight);
	  _histos["y_jetp"]->fill(fabs(yj),weight);
	  _histos["y*_jetp"]->fill(fabs(yjstar),weight);
	}
	else if (j==jm) {
	  _histos["pT_jetm"]->fill(ptj,weight);
	  _histos["y_jetm"]->fill(fabs(yj),weight);
	  _histos["y*_jetm"]->fill(fabs(yjstar),weight);
	}
	else {
	  n++;
	  switch (n) {
	  case 1: 
	    _histos["pT_jet1"]->fill(ptj,weight);
	    _histos["y_jet1"]->fill(fabs(yj),weight);
	    _histos["y*_jet1"]->fill(fabs(yjstar),weight);
	    break;
	  case 2: 
	    _histos["pT_jet2"]->fill(ptj,weight);
	    _histos["y_jet1"]->fill(fabs(yj),weight);
	    _histos["y*_jet2"]->fill(fabs(yjstar),weight);
	    break;
	  }
	}
      }
      _histos["HT_jets"]->fill(HT,weight);
      HT += p_leptons[0].pT()+p_leptons[1].pT()+missET;
      _histos["HT"]->fill(HT,weight);
      double phijj = abs(p_jets[jp].phi()-p_jets[jm].phi());
      if (phijj>PI) phijj = abs(phijj-PI);
      _histos["dphi_jetpm"]->fill(phijj,weight);
      _histos["missET_VBF"]->fill(missET,weight);
      _histos["ptL1_VBF"]->fill(p_leptons[0].pT(),weight);
      _histos["ptL2_VBF"]->fill(p_leptons[1].pT(),weight);
      _histos["dphiLL_VBF"]->fill(deltaPhi(p_leptons[0],
					   p_leptons[1]),weight);
      _histos["detaLL_VBF"]->fill(fabs(p_leptons[0].eta()-
				       p_leptons[1].eta()),weight);
    }
  private:
    std::map<std::string,AIDA::IHistogram1D*> _histos;
    
    void init_histos() { 
      _histos["Njets_low"]      = bookHistogram1D("Njets_low", 11, -0.5, 10.5);
      _histos["Njets_high"]     = bookHistogram1D("Njets_high", 11, -0.5, 10.5);
      
      _histos["missET_all"]     = 
	bookHistogram1D("missET_all",logspace(20.,500.0,10));
      _histos["ptL1_all"]       = 
	bookHistogram1D("ptL1_all",logspace(20.,500.0,10));
      _histos["ptL2_all"]       = 
	bookHistogram1D("ptL2_all",logspace(20.,500.0,10));
      _histos["dphiLL_all"]     = bookHistogram1D("dphiLL_all",8,0.,PI);
      _histos["detaLL_all"]     = bookHistogram1D("detaLL_all",6,0.,6.);
      _histos["missET_0_low"]   = 
	bookHistogram1D("missET_0_low",logspace(20.,500.0,10));
      _histos["ptL1_0_low"]     = 
	bookHistogram1D("ptL1_0_low",logspace(20.,500.0,10));
      _histos["ptL2_0_low"]     = 
	bookHistogram1D("ptL2_0_low",logspace(20.,500.0,10));
      _histos["dphiLL_0_low"]   = bookHistogram1D("dphiLL_0_low",8,0.,PI);
      _histos["detaLL_0_low"]   = bookHistogram1D("detaLL_0_low",6,0.,6.);
      _histos["missET_1_low"]   = 
	bookHistogram1D("missET_1_low",logspace(20.,500.0,10));
      _histos["ptL1_1_low"]     = 
	bookHistogram1D("ptL1_1_low",logspace(20.,500.0,10));
      _histos["ptL2_1_low"]     = 
	bookHistogram1D("ptL2_1_low",logspace(20.,500.0,10));
      _histos["dphiLL_1_low"]   = bookHistogram1D("dphiLL_1_low",8,0.,PI);
      _histos["detaLL_1_low"]   = bookHistogram1D("detaLL_1_low",6,0.,6.);
      _histos["missET_0_high"]  = 
	bookHistogram1D("missET_0_high",logspace(20.,500.0,10));
      _histos["ptL1_0_high"]    = 
	bookHistogram1D("ptL1_0_high",logspace(20.,500.0,10));
      _histos["ptL2_0_high"]    = 
	bookHistogram1D("ptL2_0_high",logspace(20.,500.0,10));
      _histos["dphiLL_0_high"]  = bookHistogram1D("dphiLL_0_high",8,0.,PI);
      _histos["detaLL_0_high"]  = bookHistogram1D("detaLL_0_high",6,0.,6.);
      _histos["missET_VBF"]     = 
	bookHistogram1D("missET_VBF",logspace(20.,500.0,10));
      _histos["ptL1_VBF"]       = 
	bookHistogram1D("ptL1_VBF",logspace(20.,500.0,10));
      _histos["ptL2_VBF"]       = 
	bookHistogram1D("ptL2_VBF",logspace(20.,500.0,10));
      _histos["dphiLL_VBF"]     = bookHistogram1D("dphiLL_VBF",8,0.,PI);
      _histos["detaLL_VBF"]     = bookHistogram1D("detaLL_VBF",6,0.,6.);

      _histos["pT_jetp"]        = 
	bookHistogram1D("pT_jetp",logspace(20.,500.,10));
      _histos["pT_jetm"]        = 
	bookHistogram1D("pT_jetm",logspace(20.,500.,10));
      _histos["pT_jet1"]        = 
	bookHistogram1D("pT_jet1",logspace(20.,500.,10));
      _histos["pT_jet2"]        = 
	bookHistogram1D("pT_jet2",logspace(20.,500.,10));
      _histos["y_jetp"]         = bookHistogram1D("y_jetp",5,0.,5.0);
      _histos["y_jetm"]         = bookHistogram1D("y_jetm",5,0.,5.0);
      _histos["y_jet1"]         = bookHistogram1D("y_jet1",5,0.,5.0);
      _histos["y_jet2"]         = bookHistogram1D("y_jet2",5,0.,5.0);
      _histos["y*_jetp"]        = bookHistogram1D("y*_jetp",5,0.,5.0);
      _histos["y*_jetm"]        = bookHistogram1D("y*_jetm",5,0.,5.0);
      _histos["y*_jet1"]        = bookHistogram1D("y*_jet1",5,0.,5.0);
      _histos["y*_jet2"]        = bookHistogram1D("y*_jet2",5,0.,5.0);
      _histos["pT_jet_veto"]    = 
	bookHistogram1D("pT_jet_veto",logspace(20.,500.,10));
      _histos["y_jet_veto"]     = bookHistogram1D("y_jet_veto",5,0.,5.0);
      _histos["y*_jet_veto"]    = bookHistogram1D("y*_jet_veto",5,0.,5.0);


      _histos["HT_jets"]        = 
	bookHistogram1D("HT_jets",logspace(50,2000.,20));
      _histos["HT"]             = 
	bookHistogram1D("HT",logspace(50,2000.,20));
      _histos["M_jetpm"]        = 
	bookHistogram1D("M_jetpm",logspace(100.,2000.,50));
      _histos["dy_jetpm"]       = bookHistogram1D("dy_jetpm",10,0.,10.0);
      _histos["dphi_jetpm"]     = bookHistogram1D("dphi_jetpm",8,0.,PI);
      _histos["weights"]        = bookHistogram1D("weights",10,-0.5,9.5);
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

    int FilterJets(const Jets& jets,
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
      return p_jets.size();
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
	     hit=_histos.begin();
	   hit!=_histos.end();hit++) {
	scale(hit->second, crossSection()/sumOfWeights());
      }
    }
  };

  // The hook for the plugin system
  AnalysisBuilder<MC_WplusWplus> plugin_MC_WplusWplus;
}
