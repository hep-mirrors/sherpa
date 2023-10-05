// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <iostream>
#include <fstream>

namespace Rivet {

  /// @brief Add a short analysis description here
  class DiffractiveAnalysis : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DiffractiveAnalysis);
    int cont = 0;
    int DDcont = 0;
    int SD1cont = 0;
    int SD2cont = 0;
    int ELcont = 0;
    //int elastic_events = 0;
    //int non_elastic_events = 0;

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs(Cuts::open());
      declare(fs,"allParticles");

      book(_h["DD"], "DD", 1, 0, 1);
      book(_h["SD1"], "SD1", 1, 0, 1);
      book(_h["SD2"], "SD2", 1, 0, 1);
      book(_h["elastic"], "elastic", 1, 0, 1);

      book(_h["rapidity-all"], "rapidity-all", 100, -10, 10);
      book(_h["DD-rapidity"], "DD-rapidity", 100, -10, 10);
      book(_h["DD-gap"], "DD-gap", 30, 0, 15);
      book(_h["gap"], "gap", 30, 0, 15);
      book(_h["pt"], "pt", 100, 0, 20);

      book(_h["SDdiffmass2"],"SDdiffmass2",150,0,150);
      book(_h2["DDdiffmass2"],"DDdiffmass2",200,0.,200000.,200,0.,200000.);
      //book(_h2["DDdiffmass2"],"DDdiffmass2",200,0.,20.,200,0.,20.);
      book(_h2["byhemDDdiffmass2"],"byhemDDdiffmass2",200,0.,200000.,200,0.,200000.);
      //book(_h2["byhemDDdiffmass2"],"byhemDDdiffmass2",200,0.,40.,200,0.,40.);
      book(_h["elastic_pt2"],"elastic_pt2",200,0,5);

      book(_c["DDevents"],"DDevents");
      book(_c["SD1events"],"SD1events");
      book(_c["SD2events"],"SD2events");
      book(_c["ELevents"],"ELevents");
      book(_c["ALLevents"],"ALLevents");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      cont++;
      _c["ALLevents"]->fill();
      int numcharged = 0;
      cout << "\n\n-----> Event " << cont << endl;
      std::map<double,Particle> particleRapidityMap;
      std::map<double,Particle> POSparticleAbsRapMap;
      std::map<double,Particle> NEGparticleAbsRapMap;
      Particles allparticles = apply<FinalState>(event, "allParticles").particles();

      int numProtonsPOS = 0;
      int numProtonsNEG = 0;
      FourMomentum pos_hem_mom, neg_hem_mom;
      multiply(pos_hem_mom,0.);
      multiply(neg_hem_mom,0.);
      for(Particle part:allparticles) {
        if (part.pid() != 22) {
          particleRapidityMap.insert(pair<double,Particle>(part.rapidity(),part));
          _h["rapidity-all"]->fill(part.eta());
          _h["pt"]->fill(sqrt(part.pT2()));
          if(part.isCharged()) numcharged++;
          if(part.rapidity() > 0) {
            pos_hem_mom += part.momentum();
            POSparticleAbsRapMap.insert(pair<double,Particle>(part.rapidity(),part));
            if(part.pid() == PID::PROTON) numProtonsPOS++;
          }
          else {
            neg_hem_mom += part.momentum();
            NEGparticleAbsRapMap.insert(pair<double,Particle>(-1*part.rapidity(),part));
            if(part.pid() == PID::PROTON) numProtonsNEG++;
          }
        }
      }
      std::map<double,Particle>::iterator firstPOS = POSparticleAbsRapMap.begin();
      std::map<double,Particle>::iterator lastPOS = POSparticleAbsRapMap.end();
      --lastPOS;
      std::map<double,Particle>::iterator firstNEG = NEGparticleAbsRapMap.begin();
      std::map<double,Particle>::iterator lastNEG = NEGparticleAbsRapMap.end();
      --lastNEG;

      int numPosParts = POSparticleAbsRapMap.size();
      int numNegParts = NEGparticleAbsRapMap.size();
      double gapPoles = firstPOS->first + firstNEG->first; //gap between most "central" particles of each pole
      _h["gap"] -> fill(gapPoles);

      if((numProtonsPOS == 1 && numProtonsNEG == 1) && particleRapidityMap.size() == 2) {
        //only protons out: EL
        _c["ELevents"]->fill();
        _h["elastic"] -> fill(0);
        _h["elastic_pt2"] -> fill(pos_hem_mom.pT2());
        ELcont++;
        
      }
      else if((numNegParts > 1 && numPosParts == 1) && numProtonsPOS == 1){
        //a proton in the positive side and more than one particle in the negative one: SD1 (1 or 2 is arbitrary)
        _c["SD1events"]->fill();
        _h["SD1"] -> fill(0);
        SD1cont++;
        _h["SDdiffmass2"] -> fill(neg_hem_mom.mass2());
      }
      else if((numNegParts == 1 && numPosParts > 1) && numProtonsNEG == 1){ 
        //a proton in the positive side and more than one particle in the negative one: SD2 (1 or 2 is arbitrary)
        _c["SD2events"]->fill();
        _h["SD2"] -> fill(0);
        SD2cont++;
        _h["SDdiffmass2"] -> fill(pos_hem_mom.mass2());
      }
      else {
        //DD
        _c["DDevents"]->fill();
        _h["DD"] -> fill(0);
        _h["DD-gap"] -> fill(gapPoles);
        double largermass(0.), smallermass(0.);
        _h2["byhemDDdiffmass2"]->fill(pos_hem_mom.mass2(),neg_hem_mom.mass2());
        if(pos_hem_mom.mass2() > neg_hem_mom.mass2()) {
          largermass = pos_hem_mom.mass2();
          smallermass = neg_hem_mom.mass2();
        }
        else {
          largermass = neg_hem_mom.mass2();
          smallermass = pos_hem_mom.mass2();
        }
        _h2["DDdiffmass2"]->fill(largermass,smallermass);
        for(std::map<double,Particle>::iterator it = particleRapidityMap.begin(); it != particleRapidityMap.end(); ++it) {
          _h["DD-rapidity"] -> fill(it->first);
        }
        DDcont++;
      }
      
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double sumDDwts = _c["DDevents"]->val();
      double sumSD1wts = _c["SD1events"]->val();
      double sumSD2wts = _c["SD2events"]->val();
      double sumELwts = _c["ELevents"]->val();
      double sumALLwts = _c["ALLevents"]->val();
      scale(_h["DD"], crossSection()/sumDDwts);
      scale(_h["SD1"], crossSection()/sumSD1wts);
      scale(_h["SD2"], crossSection()/sumSD2wts);
      scale(_h["elastic"], crossSection()/sumELwts);
      normalize(_h["rapidity-all"]);
      normalize(_h["DD-rapidity"]);
      normalize(_h["DD-gap"]);
      normalize(_h["gap"]);
      normalize(_h["pt"]);
      normalize(_h["SDdiffmass2"]);
      normalize(_h2["DDdiffmass2"]);
      normalize(_h2["byhemDDdiffmass2"]);
      normalize(_h["elastic_pt2"]);
      cout << "#EL events = " << ELcont << endl;
      cout << "#SD1 events = " << SD1cont << endl;
      cout << "#SD2 events = " << SD2cont << endl;
      cout << "#DD events = " << DDcont << endl;
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    ///@}


  };

  RIVET_DECLARE_PLUGIN(DiffractiveAnalysis);

}
