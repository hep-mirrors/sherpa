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
    double _mW, _mt, _mZ, _sigMT, _sigMW, _sigMZ;
    
    AIDA::IHistogram1D *_h_njets; //HT?
    AIDA::IHistogram1D *_h_mass_ll_ss,*_h_mass_ll_os,*_h_mass_lll;
    AIDA::IHistogram1D *_h_pT_l1,*_h_pT_l2,*_h_pT_l3,*_h_missET;
    AIDA::IHistogram1D *_h_eta_l1,*_h_eta_l2,*_h_eta_l3; 
    AIDA::IHistogram1D *_h_pT_t1,*_h_pT_t2,*_h_pT_W,*_h_pT_tt,*_h_pT_ttW;
    AIDA::IHistogram1D *_h_mass_tt,*_h_mass_ttW;
    AIDA::IHistogram1D *_h_rap_t1,*_h_rap_t2,*_h_rap_W,*_h_rap_tt,*_h_rap_ttW;
    AIDA::IHistogram1D *_h_pT_b1,*_h_pT_b2; //,*_h_pT_lj;
    //   AIDA::IHistogram1D *_h_eta_b1,*_h_eta_b2,*_h_eta_lj; 
    
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
      _h_pT_lj      = bookHistogram1D("pT_lj", logspace(10.0, 500.0, 50));
      _h_eta_b1     = bookHistogram1D("eta_b1", 20, -_jeta, _jeta);
      _h_eta_b2     = bookHistogram1D("eta_b2", 20, -_jeta, _jeta); 
      _h_eta_lj     = bookHistogram1D("eta_lj", 20, -_jeta, _jeta); 
    }
  public:
    MC_TTW() : 
      Analysis("MC_TTW"),
      _leta(2.5), _mupt(10.), _ept(10.), _lpt(min(_ept,_mupt)),
      _missET(50.), _isoE_lj(0.05), _isoR_lj(0.2), 
      _jeta(5.), _jR(0.4), _jpt(20.), 
      _mW(80.419), _mt(172.5), _mZ(90.188),
      _sigMT(17.), _sigMW(10.), _sigMZ(10.)
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
      FourMomentum W1, W2, W3;
      FourMomentum t1, t2;
      FourMomentum W1_test, W2_test, W3_test;
      FourMomentum t1_test, t2_test;
      FourMomentum lepmoms[3], neumoms[3], bmoms[2];
      FourMomentum lepton_ss1, lepton_ss2, lepton_os;
      for (size_t i=0;i<3;i++) {
       	lepmoms[i] = leptons[i].momentum();
       	neumoms[i] = neutrinos[i].momentum();
      }
      for (size_t i=0;i<2;i++) {
       	bmoms[i] = bjets[i].momentum();
      }
      float chi_t1;
      float chi_t2;
      float chi_W1;
      float chi_W2;
      float chi_W3;
      
      float chisq=1.e12;
      float chisq_test;
      for (size_t ltag1=0;ltag1<3;ltag1++) {
      	for (size_t nutag1=0;nutag1<3;nutag1++) {
      	  for (size_t btag=0;btag<2;btag++) {
      	    W1_test = lepmoms[ltag1]+neumoms[nutag1];
      	    t1_test = lepmoms[ltag1]+neumoms[nutag1]+bmoms[btag];
      	    for (size_t ltag2=0;ltag2<2;ltag2++) {
      	      if (ltag1==ltag2) continue;
	      if (PID::charge(leptons[ltag1].pdgId()==PID::charge(leptons[ltag2].pdgId()))){
		lepton_ss1 = lepmoms[ltag1];
		lepton_ss2 = lepmoms[ltag2];
		lepton_os  = lepmoms[3-ltag1-ltag2];
		continue;
	      } else {
      		for (size_t nutag2=0;nutag2<3;nutag2++){
      		  if (nutag1==nutag2) continue;
      		  else {
      		    W2_test=lepmoms[ltag2]+neumoms[nutag2];
      		    t2_test=lepmoms[ltag2]+neumoms[nutag2]+bmoms[1-btag];
      		    W3_test=lepmoms[3-ltag1-ltag2]+neumoms[3-nutag1-nutag2];
      		    chi_t1 = sqr((t1_test.mass()-_mt)/_sigMT);
      		    chi_t2 = sqr((t2_test.mass()-_mt)/_sigMT);
      		    chi_W1 = sqr((W1_test.mass()-_mW)/_sigMW);
      		    chi_W2 = sqr((W2_test.mass()-_mW)/_sigMW);
      		    chi_W3 = sqr((W3_test.mass()-_mW)/_sigMW);
      		    chisq_test=chi_t1+chi_t2+chi_W1+chi_W2+chi_W3;
      		    if (chisq_test<chisq){
      		      chisq=chisq_test;
      		      W1=W1_test;
      		      W2=W2_test;
      		      W3=W3_test;
      		      t1=t1_test;
      		      t2=t2_test;
      		    }
      		  }
      		}
      	      }
      	    }
      	  }
      	}
      }
      
      // std::cout << "W1  "<<W1<<std::endl;
      //  std::cout << "Done W/t" << std::endl;
      _h_njets->fill(alljets.size(),weight);
      _h_mass_ll_ss->fill(lepton_ss1.mass()/GeV+lepton_ss2.mass()/GeV,weight);
      float Zmass2(pow(_mZ,2));
      float test;
      test=fabs((lepton_ss1+lepton_os)*(lepton_ss1+lepton_os))-Zmass2;
      if (fabs((lepton_ss2+lepton_os)*(lepton_ss2+lepton_os))-Zmass2<test)
	_h_mass_ll_os->fill((lepton_ss2+lepton_os).mass(),weight); 
      else
	_h_mass_ll_os->fill((lepton_ss1+lepton_os).mass(),weight); 
      _h_mass_lll->fill((lepmoms[0]+lepmoms[1]+lepmoms[2]).mass()/GeV,weight);
      _h_pT_l1->fill(lepmoms[0].pT()/GeV,weight);
      _h_pT_l2->fill(lepmoms[1].pT()/GeV,weight);
      _h_pT_l3->fill(lepmoms[2].pT()/GeV,weight);
      _h_eta_l1->fill(lepmoms[0].pseudorapidity(),weight);
      _h_eta_l2->fill(lepmoms[1].pseudorapidity(),weight);
      _h_eta_l3->fill(lepmoms[2].pseudorapidity(),weight);
      _h_missET->fill(neumoms[0].Et()/GeV,weight); // change

      //std::cout << "1"<<std::endl;
      _h_pT_t1->fill(t1.pT()/GeV,weight);
      _h_pT_t2->fill(t2.pT()/GeV,weight);
      _h_pT_W->fill(W3.pT()/GeV,weight);
      _h_pT_tt->fill((t1+t2).pT()/GeV,weight);
      _h_pT_ttW->fill((t1+t2+W3).pT()/GeV,weight);
      _h_mass_tt->fill((t1+t2).mass()/GeV,weight);
      _h_mass_ttW->fill((t1+t2+W3).mass()/GeV,weight);
      _h_rap_t1->fill(t1.y(),weight);
      _h_rap_t2->fill(t2.y(),weight);
      _h_rap_W->fill(W3.y(),weight);
      _h_rap_tt->fill((t1+t2).y(),weight);
      _h_rap_ttW->fill((t1+t2+W3).y(),weight);
      
      //     std::cout <<"2"<<std::endl;
      _h_pT_b1->fill(bjets[0].momentum().pT()/GeV,weight);
      //      if (bjets.size()>0){
      _h_pT_b2->fill(bjets[1].momentum().pT()/GeV,weight);
      //    }
      //    std::cout << "2.5"<<std::endl;
      //if (ljets.size()>0){
      //_h_pT_lj->fill(ljets[0].momentum().pT()/GeV,weight);
      // }
      //      std::cout <<"2.75"<<std::endl;
      //      if (bjets.size()>0){
      // _h_eta_b1->fill(bjets[0].momentum().pseudorapidity(),weight); 
      //      }
      //     if (bjets.size()>1){
      // _h_eta_b2->fill(bjets[1].momentum().pseudorapidity(),weight);  
      //     }
      //      std::cout << "3"<<std::endl;
      //if (ljets.size()>0){
      //	_h_eta_lj->fill(ljets[0].momentum().pT()/GeV,weight);
      //}  
      //     std::cout << "Filled!"<<std::endl;
      //   }
    }
    
    void finalize() {
      // std::cout << "Finalising"<<std::endl;
      normalize(_h_njets); 
      normalize(_h_mass_ll_ss); 
      normalize(_h_mass_ll_os); 
      normalize(_h_mass_lll); 
      normalize(_h_pT_l1); 
      normalize(_h_pT_l2); 
      normalize(_h_pT_l3); 
      // std::cout << "first set"<<std::endl;
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
      //std::cout << "second set"<<std::endl;
      normalize(_h_rap_t2); 
      normalize(_h_rap_W); 
      normalize(_h_rap_tt); 
      normalize(_h_rap_ttW); 
      normalize(_h_pT_b1); 
      normalize(_h_pT_b2); 
      // normalize(_h_pT_lj); 
      // std::cout << "b eta"<<std::endl;
      // normalize(_h_eta_b1); 
      //std::cout << "b1 eta"<<std::endl;
      //normalize(_h_eta_b2); 
      //std::cout << "b2 eta"<<std::endl;
      // normalize(_h_eta_lj); 
      // std::cout << "Done"<<std::endl;
    }
  };

  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TTW);

}
