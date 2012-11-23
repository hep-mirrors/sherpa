// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include <map>

namespace Rivet {


  /// @brief MC validation analysis for higgs [-> tau tau] + jets events
  class MC_HJETS_STABLE_CUTS : public Analysis {
  private:
    double _jeta, _jR, _jpt, _jjdy, _jjmass;
    std::map<std::string,AIDA::IHistogram1D *> histos;
  public:
    MC_HJETS_STABLE_CUTS() : 
      Analysis("MC_HJETS_STABLE_CUTS"),
      _jeta(5.), _jR(0.4), _jpt(25.),
      _jjdy(2.8), _jjmass(400.)
    {}
    
    void init() {
      FinalState fs;
      IdentifiedFinalState hfinder(-10.,10.,0.*GeV);
      hfinder.acceptId(HIGGS);
      VetoedFinalState hexcl(FinalState(-_jeta,_jeta,0.*GeV));
      hexcl.addVetoDetail(HIGGS,0.*GeV,14000.*GeV);
      addProjection(hfinder, "Hfinder");
      addProjection(FastJets(hexcl, FastJets::ANTIKT, _jR), "Jets");

      inithistos();
    }


    void inithistos() {
      histos["H_pt_peak"]      = bookHistogram1D("H_pT_peak", 20, 1.0, 21.0);
      histos["H_pt_low"]       = bookHistogram1D("H_pT_low", 20, 0.0, 50.0);
      histos["H_pt_long"]      = bookHistogram1D("H_pT_long", 50, 0.0, 500.0);
      histos["H_pt"]           = bookHistogram1D("H_pT", 20, 0.0, 100.0);
      histos["jet1_pt"]        = bookHistogram1D("jet1_pt",7,25.,200.);
      histos["jet2_pt"]        = bookHistogram1D("jet2_pt",7,25.,200.);
      histos["jet3_pt"]        = bookHistogram1D("jet3_pt",7,25.,200.);
      histos["jet1_y"]         = bookHistogram1D("jet1_y",10,-5.,5.);
      histos["jet2_y"]         = bookHistogram1D("jet2_y",10,-5.,5.);
      histos["jet3_y"]         = bookHistogram1D("jet3_y",10,-5.,5.);
      histos["jet12_dy"]       = bookHistogram1D("jet12_dy",8,0.,8.);
      histos["jet12_dphi"]     = bookHistogram1D("jet12_dphi",10,0.,PI);
      histos["jet12_mass"]     = bookHistogram1D("jet12_mass",20,0.,800.);

      histos["H_pt_2jet"]      = bookHistogram1D("H_pT_2j", 20, 0.0, 100.0);
      histos["jet1_pt_2jet"]   = bookHistogram1D("jet1_pt_2j",7,25.,200.);
      histos["jet1_y_2jet"]    = bookHistogram1D("jet1_y_2j",10,-5.,5.);

      histos["H_pt_WBF"]       = bookHistogram1D("H_pT_WBF", 20, 0.0, 100.0);
      histos["jet1_pt_WBF"]    = bookHistogram1D("jet1_pt_WBF",7,25.,200.);
      histos["jet2_pt_WBF"]    = bookHistogram1D("jet2_pt_WBF",7,25.,200.);
      histos["jet3_pt_WBF"]    = bookHistogram1D("jet3_pt_WBF",7,25.,200.);
      histos["jet1_y_WBF"]     = bookHistogram1D("jet1_y_WBF",10,-5.,5.);
      histos["jet2_y_WBF"]     = bookHistogram1D("jet2_y_WBF",10,-5.,5.);
      histos["jet3_y_WBF"]     = bookHistogram1D("jet3_y_WBF",10,-5.,5.);
      histos["jet3_y*_WBF"]    = bookHistogram1D("jet3_y*_WBF",10,-5.,5.);
      histos["jet12_dy_WBF"]   = bookHistogram1D("jet12_dy_WBF",8,0.,8.);
      histos["jet12_dphi_WBF"] = bookHistogram1D("jet12_dphi_WBF",10,0.,PI);
      histos["jet12_mass_WBF"] = bookHistogram1D("jet12_mass_WBF",20,0.,800.);
    }

    /// Do the analysis
    void analyze(const Event & e) {
      const IdentifiedFinalState hfinder = 
	applyProjection<IdentifiedFinalState>(e, "Hfinder");
      if (hfinder.size()!=1) vetoEvent;
      const double weight = e.weight();

      FourMomentum hmom(hfinder.particles()[0].momentum());
      const FastJets& jetpro = applyProjection<FastJets>(e, "Jets");
      const Jets& jets = jetpro.jetsByPt(_jpt*GeV);

      histos["H_pt"]->fill(hmom.pT(),weight);
      histos["H_pt_low"]->fill(hmom.pT(),weight);
      histos["H_pt_peak"]->fill(hmom.pT(),weight);
      histos["H_pt_long"]->fill(hmom.pT(),weight);
      if (jets.size()>0) {
	histos["jet1_pt"]->fill(jets[0].momentum().pT(),weight);
	histos["jet1_y"]->fill(jets[0].momentum().rapidity(),weight);
	if (jets.size()>1) {
	  histos["H_pt_2jet"]->fill(hmom.pT(),weight);
	  histos["jet1_pt_2jet"]->fill(jets[0].momentum().pT(),weight);
	  histos["jet1_y_2jet"]->fill(jets[0].momentum().rapidity(),weight);
	  histos["jet2_pt"]->fill(jets[1].momentum().pT(),weight);
	  histos["jet2_y"]->fill(jets[1].momentum().rapidity(),weight);
	  double yjj(fabs(jets[0].momentum().rapidity()-
			  jets[1].momentum().rapidity()));
	  double massjj((jets[0].momentum()+jets[1].momentum()).mass());
	  histos["jet12_dy"]->fill(yjj,weight);
	  histos["jet12_dphi"]->fill(deltaPhi(jets[0].momentum(),
					      jets[1].momentum()),weight);
	  histos["jet12_mass"]->fill((jets[0].momentum()+
				      jets[1].momentum()).mass(),weight);
	  if (jets.size()>2) {
	    histos["jet3_pt"]->fill(jets[2].momentum().pT(),weight);
	    histos["jet3_y"]->fill(jets[2].momentum().rapidity(),weight);
	  }
	  if (yjj>_jjdy && massjj>_jjmass) {
	    histos["H_pt_WBF"]->fill(hmom.pT(),weight);
	    histos["jet1_pt_WBF"]->fill(jets[0].momentum().pT(),weight);
	    histos["jet1_y_WBF"]->fill(jets[0].momentum().rapidity(),weight);
	    histos["jet2_pt_WBF"]->fill(jets[1].momentum().pT(),weight);
	    histos["jet2_y_WBF"]->fill(jets[1].momentum().rapidity(),weight);
	    histos["jet12_dy_WBF"]->fill(yjj,weight);
	    histos["jet12_dphi_WBF"]->fill(deltaPhi(jets[0].momentum(),
						    jets[1].momentum()),weight);
	    histos["jet12_mass_WBF"]->fill((jets[0].momentum()+
					    jets[1].momentum()).mass(),weight);
	    if (jets.size()>2) {
	      histos["jet3_pt_WBF"]->fill(jets[2].momentum().pT(),weight);
	      histos["jet3_y_WBF"]->fill(jets[2].momentum().rapidity(),weight);
	      double ystar((jets[0].momentum().rapidity()+
			    jets[1].momentum().rapidity())/2.-
			   jets[2].momentum().rapidity());
	      histos["jet3_y*_WBF"]->fill(ystar,weight);
	    }
	  }
	} 
      }
    }


    /// Finalize
    void finalize() {
      double scalefactor(crossSection()/sumOfWeights());
      for (std::map<std::string,AIDA::IHistogram1D *>::iterator 
	     hit=histos.begin(); hit!=histos.end();hit++) 
	scale(hit->second,scalefactor);
    }
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HJETS_STABLE_CUTS);
}
