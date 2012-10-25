// An analysis for anomalous gauge couplings in W gamma

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {
  class MC_Anomalous_Wgamma : public Analysis {
  private:
    int m_nbins_pt, m_nbins_eta;
    double m_ptlep_min, m_ptgam_min, m_ptmiss_min, m_eta_min, m_ptjet_min,m_etajet_min;
    double m_deltaR, m_deltaR_gamlep, m_deltaR_lepjet, m_deltaR_gamjet;
    AIDA::IHistogram1D * h_ptlep, * h_ptlepl, * h_ptgam, * h_ptgaml, * h_ptmiss;
    AIDA::IHistogram1D * h_etalep, * h_etagam, * h_signed_delta_eta_gamlep;
    AIDA::IHistogram1D * h_ptjet, * h_etajet, * h_multjet;
    AIDA::IHistogram1D * h_massgamlep, * h_masseffgamlepptmiss;
    AIDA::IHistogram1D * h_cross_lepgam, * h_cross_lepgam_h;
    AIDA::IHistogram1D * h_cross_lepgam_long, * h_cross_lepgam_trans; 
    AIDA::IHistogram1D * h_cross_lepgam_long_h, * h_cross_lepgam_trans_h; 

    void init_histos() {     
      h_ptlep  = bookHistogram1D("pt_lepton",m_nbins_pt,m_ptlep_min,m_ptlep_min+100*5.);
      h_ptlepl = bookHistogram1D("pt_lepton_low",2*m_nbins_pt,m_ptlep_min,50.);
      h_ptgam  = bookHistogram1D("pt_photon",m_nbins_pt,m_ptgam_min,m_ptgam_min+100*5.);
      h_ptgaml = bookHistogram1D("pt_photon_low",2*m_nbins_pt,m_ptgam_min,50);
      h_ptmiss = bookHistogram1D("miss_pt",m_nbins_pt,m_ptmiss_min,m_ptmiss_min+100*5.);
      h_etalep = bookHistogram1D("eta_lepton",m_nbins_eta,-m_eta_min,m_eta_min);
      h_etagam = bookHistogram1D("eta_photon",m_nbins_eta,-m_eta_min,m_eta_min);
      h_signed_delta_eta_gamlep = bookHistogram1D("delta_eta",m_nbins_eta,-2.*m_eta_min,2.*m_eta_min);
      h_ptjet   = bookHistogram1D("pt_jet",m_nbins_pt,m_ptjet_min,m_ptjet_min+100*5.);
      h_etajet  = bookHistogram1D("eta_jet",m_nbins_eta,-m_etajet_min,m_etajet_min);
      h_multjet = bookHistogram1D("mult_jet",5,0,5);
      h_massgamlep           = bookHistogram1D("mass_gamma_lepton",m_nbins_pt,40.,540.);
      h_masseffgamlepptmiss  = bookHistogram1D("eff_mass",m_nbins_pt,40.,540.);
      h_cross_lepgam         = bookHistogram1D("cross_gamma_lepton",m_nbins_eta,0.,1.);
      h_cross_lepgam_long    = bookHistogram1D("cross_long_gamma_lepton",m_nbins_eta/2,-1.,1.);
      h_cross_lepgam_trans   = bookHistogram1D("cross_trans_gamma_lepton",m_nbins_eta/2,-1.,1.);
      h_cross_lepgam_h       = bookHistogram1D("cross_gamma_lepton_highmass",m_nbins_eta,0.,1.);
      h_cross_lepgam_long_h  = bookHistogram1D("cross_long_gamma_lepton_highmass",m_nbins_eta/2,-1.,1.);
      h_cross_lepgam_trans_h = bookHistogram1D("cross_trans_gamma_lepton_highmass",m_nbins_eta/2,-1.,1.);
    }
  public:
    MC_Anomalous_Wgamma() :
      Analysis("MC_Anomalous_Wgamma"),
      m_nbins_pt(20), m_nbins_eta(m_nbins_pt/2), 
      m_ptlep_min(20.), m_ptgam_min(20.), m_ptmiss_min(20.), m_eta_min(2.5), 
      m_ptjet_min(30.), m_etajet_min(5.), m_deltaR(0.4), 
      m_deltaR_gamlep(0.2), m_deltaR_lepjet(0.2+m_deltaR), m_deltaR_gamjet(0.2+m_deltaR)
      {}
    void init() {
      // projection to find the electrons and photons
      std::vector<std::pair<double, double> > etas;
      etas.push_back(make_pair(-m_eta_min,m_eta_min));
      IdentifiedFinalState elecs(etas, m_ptlep_min*GeV);
      elecs.acceptIdPair(ELECTRON);
      addProjection(elecs, "elecs");
      IdentifiedFinalState photons(etas, m_ptgam_min*GeV);
      photons.acceptId(PHOTON);
      addProjection(photons, "photons");

      // Jet finder
      FinalState fs;
      addProjection(fs, "FS");
      addProjection(FastJets(fs, FastJets::ANTIKT, m_deltaR),
		    "AntiKtJets");

      VisibleFinalState vfs;
      addProjection(vfs, "VFS");
      // All tracks (to do deltaR with leptons and photons)
      addProjection(ChargedFinalState(-m_eta_min,m_eta_min,0.5*GeV),"CFS");

      init_histos();
    }
    void analyze(const Event& event) {
      const double weight = event.weight(); 

      // candidate jets
      Jets cand_jets;
      foreach (const Jet& jet,
	       applyProjection<FastJets>(event, "AntiKtJets").jetsByPt(m_ptjet_min*GeV) ) {
        if (fabs(jet.momentum().eta())<m_etajet_min) cand_jets.push_back(jet);
      }
      // candidate electrons
      ParticleVector cand_electrons =
	applyProjection<IdentifiedFinalState>(event, "elecs").particlesByPt();
      // candidate electrons
      ParticleVector cand_photons =
	applyProjection<IdentifiedFinalState>(event, "photons").particlesByPt();
      // charged tracks for later isolation
      ParticleVector chg_tracks =
	applyProjection<ChargedFinalState>(event, "CFS").particles();
      // missing pt
      ParticleVector visibles =
	applyProjection<VisibleFinalState>(event, "VFS").particles();
      FourMomentum missmom(0.,0.,0.,0.);
      foreach (Particle p,visibles) {
	missmom -= p.momentum();
      }
      double ptmiss = missmom.pT();
      if (ptmiss<m_ptmiss_min) vetoEvent;

      if (cand_electrons.size()<1 || cand_photons.size()<1) vetoEvent;
      
      Particle & leading_electron = cand_electrons[0];
      Particle & leading_photon   = cand_photons[0];
      /*
	std::cout<<"Check electron: "<<cand_electrons[0].momentum()<<".\n";
	std::cout<<"Check photon:   "<<cand_photons[0].momentum()<<".\n";      
	std::cout<<"Check iso: egamma    :"
	<<deltaR(leading_photon,leading_electron)<<" < "
	<<m_deltaR_gamlep<<".\n";
      */
      if (deltaR(leading_photon,leading_electron)<m_deltaR_gamlep) vetoEvent;
      
      FourMomentum plep(leading_electron.momentum());
      FourMomentum pgam(leading_photon.momentum());
      double ptlep(plep.pT()), etalep(plep.eta());
      double ptgam(pgam.pT()), etagam(pgam.eta());
      double ptjet;

      foreach(Jet jet,cand_jets) {
	ptjet = jet.momentum().pT();
	/*
	  std::cout<<"Jet momentum :"<<jet.momentum()<<", "
	  <<jet.momentum().pT()<<".\n";
	  std::cout<<"Check iso: e jet     :"
	  <<deltaR(jet,leading_electron)<<" < "
	  <<m_deltaR_lepjet<<".\n";
	  std::cout<<"Check iso: gamma jet :"
	  <<deltaR(jet,leading_photon)<<" < "
	  <<m_deltaR_gamjet<<".\n";
	*/
	if (deltaR(jet,leading_electron)<m_deltaR_lepjet) {
	  if (fabs(1.-ptlep/ptjet)>0.05) vetoEvent;
	  //else std::cout<<"electron jet!\n";
	}
	if (deltaR(jet,leading_photon)<m_deltaR_gamjet) {
	  if (fabs(1.-ptgam/ptjet)>0.05) vetoEvent;
	  //else std::cout<<"photonic jet!\n";
	}
	if (cand_electrons.size()>1) {
	  for (size_t i=1;i<cand_electrons.size();i++) {
	    if (deltaR(jet,cand_electrons[i])>m_deltaR_lepjet) vetoEvent;
	  }
	}
	if (cand_photons.size()>1) {
	  for (size_t i=1;i<cand_photons.size();i++) {
	    if (deltaR(jet,cand_photons[i])>m_deltaR_gamjet) vetoEvent;
	  }
	}
      }
      Vector3 vlep(plep.vector3()), vgam(pgam.vector3());
      // need to further isolate photon and electron vs. gunk tracks.
      
      int charge(PID::threeCharge(leading_electron.pdgId())/3);
      double mass((plep+pgam).mass()/GeV), mass2(mass*mass);
      Vector3 vcross_gamlep(charge*cross(vlep.unit(),vgam.unit()));
      double  cross_gamlep(vcross_gamlep.rho());
      double  cross_gamleplong(dot(vcross_gamlep,Vector3(0.,0.,1.)));

      h_ptlep->fill(ptlep,weight);
      h_ptlepl->fill(ptlep,weight);
      h_etalep->fill(etalep,weight);
      h_ptmiss->fill(ptmiss,weight);
      h_ptgam->fill(ptgam,weight);
      h_ptgaml->fill(ptgam,weight);
      h_etagam->fill(etagam,weight);
      h_signed_delta_eta_gamlep->fill(charge*etalep-etagam,weight);
      h_multjet->fill(cand_jets.size(),weight);
      h_massgamlep->fill(mass,weight);
      foreach(Jet jet,cand_jets) {
	h_ptjet->fill(jet.momentum().pT(),weight);
	h_etajet->fill(jet.momentum().eta(),weight);
      }
      h_cross_lepgam->fill(cross_gamlep,weight);
      h_cross_lepgam_long->fill(cross_gamleplong,weight);
      if (mass>100.) {
	h_cross_lepgam_h->fill(cross_gamlep,weight);
	h_cross_lepgam_long_h->fill(cross_gamleplong,weight);
      }
    }
    void finalize() {
      scale(h_ptlep,crossSection()/sumOfWeights());
      scale(h_ptlepl,crossSection()/sumOfWeights());
      scale(h_ptgam,crossSection()/sumOfWeights());
      scale(h_ptgaml,crossSection()/sumOfWeights());
      scale(h_ptmiss,crossSection()/sumOfWeights());
      scale(h_etalep,crossSection()/sumOfWeights());
      scale(h_etagam,crossSection()/sumOfWeights());
      scale(h_signed_delta_eta_gamlep,crossSection()/sumOfWeights());
      scale(h_multjet,crossSection()/sumOfWeights());
      scale(h_ptjet,crossSection()/sumOfWeights());
      scale(h_etajet,crossSection()/sumOfWeights());
      scale(h_massgamlep,crossSection()/sumOfWeights());
      scale(h_masseffgamlepptmiss,crossSection()/sumOfWeights());
      scale(h_cross_lepgam,crossSection()/sumOfWeights());
      scale(h_cross_lepgam_long,crossSection()/sumOfWeights());
      scale(h_cross_lepgam_trans,crossSection()/sumOfWeights());
      scale(h_cross_lepgam_h,crossSection()/sumOfWeights());
      scale(h_cross_lepgam_long_h,crossSection()/sumOfWeights());
      scale(h_cross_lepgam_trans_h,crossSection()/sumOfWeights());
    }
  };

  // The hook for the plugin system
  AnalysisBuilder<MC_Anomalous_Wgamma> plugin_MC_Anomalous_Wgamma;
}
