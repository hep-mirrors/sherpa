// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"


namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  class m_wbjet : public Analysis {

  #include "NLOHisto1D.cc"

  public:

    /// Constructor
    m_wbjet(string name="ee_wwbb_analysis")
      : Analysis(name)
    {
    }

    const double Rcut  = 0.4; // ->0.4 (not relevant in DURHAM-algorithm = eekt-algorithm in Sherpa)
    const double dcut25  = 25.;
    const double dcut100 = 100.;
    const double dcut400 = 400.;

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      using namespace Rivet::PID;

      // Init projections

      /* define jet without jet algorithm properties */
      FinalState jet_fs;
      VetoedFinalState jetfs(jet_fs);
      jetfs.addVetoPairId(TQUARK);
      jetfs.addVetoPairId(WPLUSBOSON);
      jetfs.addVetoPairId(MUON);
      jetfs.addVetoPairId(ELECTRON);
      jetfs.addVetoPairId(NU_MU);
      jetfs.addVetoPairId(NU_E);
      jetfs.addVetoPairId(ZBOSON);
      jetfs.addVetoPairId(HIGGS);
      jetfs.addVetoPairId(PHOTON);
      addProjection(jetfs, "JetFS");

      /* define jet algorithm properties */
      FastJets jet_DURHAM(jetfs, FastJets::DURHAM, Rcut); // Rcut irrelevant
      addProjection(jet_DURHAM, "DurhamJets");

      /* define W+ boson */
      IdentifiedFinalState wplus;
      wplus.acceptId(WPLUSBOSON);
      addProjection(wplus, "Wplus");

      /* define W- boson */
      IdentifiedFinalState wminus;
      wminus.acceptId(WMINUSBOSON);
      addProjection(wminus, "Wminus");

      /* define b-quark */
      IdentifiedFinalState bquark;
      bquark.acceptId(BQUARK);
      addProjection(bquark, "Bottom");

      /* define anti-b-quark */
      IdentifiedFinalState antibquark;
      antibquark.acceptId(-BQUARK);
      addProjection(antibquark, "AntiBottom");

      /* define gluon */
      IdentifiedFinalState gluon;
      gluon.acceptId(GLUON);
      addProjection(gluon, "Gluon");

      /* book Histograms: name, # bins, start, end */
      _histos["dmax25_m_Wplusb"]         = bookNLOHisto1D("dmax25_m_Wplusb", 100,50,300);
      _histos["dmax25_m_Wminusbb"]       = bookNLOHisto1D("dmax25_m_Wminusbb", 100,50,300);
      _histos["dmax25_pt-Wm"]            = bookNLOHisto1D("dmax25_pt-Wminus", 100,0,310);
      _histos["dmax25_angle-Wm-bjet"]    = bookNLOHisto1D("dmax25_angle-Wm-bjet", 100,0,M_PI);
      _histos["dmax25_angle-Wp-bjet"]    = bookNLOHisto1D("dmax25_angle-Wp-bjet", 100,0,M_PI);
      _histos["dmax25_angle-bjet-bbjet"] = bookNLOHisto1D("dmax25_angle-bjet-bbjet", 100,0,M_PI);
      _histos["dmax25_pt-bjet"]          = bookNLOHisto1D("dmax25_pt-bjet", 100,0,310);
      _histos["dmax100_m_Wplusb"]         = bookNLOHisto1D("dmax100_m_Wplusb", 100,50,300);
      _histos["dmax100_m_Wminusbb"]       = bookNLOHisto1D("dmax100_m_Wminusbb", 100,50,300);
      _histos["dmax100_pt-Wm"]            = bookNLOHisto1D("dmax100_pt-Wminus", 100,0,310);
      _histos["dmax100_angle-Wm-bjet"]    = bookNLOHisto1D("dmax100_angle-Wm-bjet", 100,0,M_PI);
      _histos["dmax100_angle-Wp-bjet"]    = bookNLOHisto1D("dmax100_angle-Wp-bjet", 100,0,M_PI);
      _histos["dmax100_angle-bjet-bbjet"] = bookNLOHisto1D("dmax100_angle-bjet-bbjet", 100,0,M_PI);
      _histos["dmax100_pt-bjet"]          = bookNLOHisto1D("dmax100_pt-bjet", 100,0,310);
      _histos["dmax400_m_Wplusb"]         = bookNLOHisto1D("dmax400_m_Wplusb", 100,50,300);
      _histos["dmax400_m_Wminusbb"]       = bookNLOHisto1D("dmax400_m_Wminusbb", 100,50,300);
      _histos["dmax400_pt-Wm"]            = bookNLOHisto1D("dmax400_pt-Wminus", 100,0,310);
      _histos["dmax400_angle-Wm-bjet"]    = bookNLOHisto1D("dmax400_angle-Wm-bjet", 100,0,M_PI);
      _histos["dmax400_angle-Wp-bjet"]    = bookNLOHisto1D("dmax400_angle-Wp-bjet", 100,0,M_PI);
      _histos["dmax400_angle-bjet-bbjet"] = bookNLOHisto1D("dmax400_angle-bjet-bbjet", 100,0,M_PI);
      _histos["dmax400_pt-bjet"]          = bookNLOHisto1D("dmax400_pt-bjet", 100,0,310);
      _histos["E-gluon"]                  = bookNLOHisto1D("E-gluon", 100,0,310);
      _histos["pt-gluon"]                 = bookNLOHisto1D("pt-gluon", 100,0,310);
      _histos["rapidity-gluon"]           = bookNLOHisto1D("rapidity-gluon", 100,-4,4);
      _histos["phi-gluon"]                = bookNLOHisto1D("phi-gluon", 100,0,M_PI);
      _histos["theta-b"]                  = bookNLOHisto1D("theta-b", 100,0,M_PI);
      _histos["theta-bb"]                 = bookNLOHisto1D("theta-bb", 100,0,M_PI);
      _histos["theta-gluon"]              = bookNLOHisto1D("theta-gluon", 100,0,M_PI);
      _histos["angle-b-gluon"]            = bookNLOHisto1D("angle-b-gluon", 100,0,M_PI);
      _histos["angle-Wm-gluon"]           = bookNLOHisto1D("angle-Wm-gluon", 100,0,M_PI);
      _histos["pt-b"]                     = bookNLOHisto1D("pt-bquark", 100,0,310);
      _histos["m_Wminusbb"]               = bookNLOHisto1D("m_Wminusbb", 100,50,300);
      _histos["m_Wplusb"]                 = bookNLOHisto1D("m_Wplusb", 100,50,300);
      _h_count_emissions                 = bookHisto1D("count_emissions", 2, 0.0, 1.0);
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState& Wplusbosons = applyProjection<IdentifiedFinalState>(event, "Wplus");
      const ParticleVector w_plus_par = Wplusbosons.particlesByPt();
      const FinalState& Wminusbosons = applyProjection<IdentifiedFinalState>(event, "Wminus");
      const ParticleVector w_minus_par = Wminusbosons.particlesByPt();

      if(w_plus_par.size() != 1 || w_minus_par.size() != 1 ){
       	cout << "wrong number of Ws -> veto"<<endl;
 	    vetoEvent;
      }
      FourMomentum Wplus_mom  = w_plus_par[0].momentum();
      FourMomentum Wminus_mom = w_minus_par[0].momentum();
 
      /* translate partons to raw pseudo-jets (no real jets yet) with user-index */
      const Particles& particles = applyProjection<FinalState>(event, "JetFS").particles();
      PseudoJets input;

      foreach(Particle p, particles){
        PseudoJet temp = p.pseudojet();
        if     (p.pdgId() == PID::BQUARK)   temp.set_user_index(5); 
        else if(p.pdgId() == -PID::BQUARK)  temp.set_user_index(-5); 
        input.push_back(temp);
      }

      /* define cluster-sequence according to jet-definition and build actual (pseudo-)jets
         from raw pseudo-jets */
      fastjet::JetDefinition ee_kt_def(fastjet::ee_kt_algorithm);
      fastjet::ClusterSequence cluster(input, ee_kt_def);
      PseudoJets jets25  = cluster.exclusive_jets(dcut25);
      PseudoJets jets100 = cluster.exclusive_jets(dcut100);
      PseudoJets jets400 = cluster.exclusive_jets(dcut400);

      /* Create b-(pseudo-)jets */
      PseudoJets b_jets25;
      PseudoJets bbar_jets25;
      PseudoJets bbbar_jets25;
      foreach(const PseudoJet & jet25_it, jets25){
        int btag = 0;
        int bbartag = 0;
        foreach(PseudoJet particle, jet25_it.constituents()){
          if(particle.user_index() == 5) btag++;
          if(particle.user_index() == -5) bbartag++;
        }
        if(btag==1 && bbartag==0) b_jets25.push_back(jet25_it);
        if(btag==0 && bbartag==1) bbar_jets25.push_back(jet25_it);
        if(btag==1 && bbartag==1) bbbar_jets25.push_back(jet25_it);
      }
      PseudoJets b_jets100;
      PseudoJets bbar_jets100;
      PseudoJets bbbar_jets100;
      foreach(const PseudoJet & jet100_it, jets100){
        int btag = 0;
        int bbartag = 0;
        foreach(PseudoJet particle, jet100_it.constituents()){
          if(particle.user_index() == 5) btag++;
          if(particle.user_index() == -5) bbartag++;
        }
        if(btag==1 && bbartag==0) b_jets100.push_back(jet100_it);
        if(btag==0 && bbartag==1) bbar_jets100.push_back(jet100_it);
        if(btag==1 && bbartag==1) bbbar_jets100.push_back(jet100_it);
      }
      PseudoJets b_jets400;
      PseudoJets bbar_jets400;
      PseudoJets bbbar_jets400;
      foreach(const PseudoJet & jet400_it, jets400){
        int btag = 0;
        int bbartag = 0;
        foreach(PseudoJet particle, jet400_it.constituents()){
          if(particle.user_index() == 5) btag++;
          if(particle.user_index() == -5) bbartag++;
        }
        if(btag==1 && bbartag==0) b_jets400.push_back(jet400_it);
        if(btag==0 && bbartag==1) bbar_jets400.push_back(jet400_it);
        if(btag==1 && bbartag==1) bbbar_jets400.push_back(jet400_it);
      }

      /* calculate observables */
      if ( b_jets25.size()==1 && bbar_jets25.size()==1 ){
        FourMomentum BottomJet_mom       =
        FourMomentum(b_jets25[0].E(),b_jets25[0].px(),b_jets25[0].py(),b_jets25[0].pz());
        FourMomentum AntiBottomJet_mom   =
        FourMomentum(bbar_jets25[0].E(),bbar_jets25[0].px(),
                     bbar_jets25[0].py(),bbar_jets25[0].pz());
        FourMomentum Top_mom             = BottomJet_mom + Wplus_mom;
        FourMomentum AntiTop_mom         = AntiBottomJet_mom + Wminus_mom;

        _histos["dmax25_m_Wplusb"]->fill( Top_mom.mass() ,event  );
        _histos["dmax25_m_Wminusbb"]->fill( AntiTop_mom.mass() ,event  );
        _histos["dmax25_pt-Wm"]->fill(Wplus_mom.pt(), event);
        _histos["dmax25_angle-Wm-bjet"]->fill(BottomJet_mom.angle(Wminus_mom) ,event);
        _histos["dmax25_angle-Wp-bjet"]->fill(BottomJet_mom.angle(Wplus_mom) ,event);
        _histos["dmax25_angle-bjet-bbjet"]->fill(BottomJet_mom.angle(AntiBottomJet_mom) ,event);
        _histos["dmax25_pt-bjet"]->fill(BottomJet_mom.pt(), event);
      }
      if ( b_jets100.size()==1 && bbar_jets100.size()==1 ){
        FourMomentum BottomJet_mom       =
        FourMomentum(b_jets100[0].E(),b_jets100[0].px(),b_jets100[0].py(),b_jets100[0].pz());
        FourMomentum AntiBottomJet_mom   =
        FourMomentum(bbar_jets100[0].E(),bbar_jets100[0].px(),
                     bbar_jets100[0].py(),bbar_jets100[0].pz());
        FourMomentum Top_mom             = BottomJet_mom + Wplus_mom;
        FourMomentum AntiTop_mom         = AntiBottomJet_mom + Wminus_mom;

        _histos["dmax100_m_Wplusb"]->fill( Top_mom.mass() ,event  );
        _histos["dmax100_m_Wminusbb"]->fill( AntiTop_mom.mass() ,event  );
        _histos["dmax100_pt-Wm"]->fill(Wplus_mom.pt(), event);
        _histos["dmax100_angle-Wm-bjet"]->fill(BottomJet_mom.angle(Wminus_mom) ,event);
        _histos["dmax100_angle-Wp-bjet"]->fill(BottomJet_mom.angle(Wplus_mom) ,event);
        _histos["dmax100_angle-bjet-bbjet"]->fill(BottomJet_mom.angle(AntiBottomJet_mom) ,event);
        _histos["dmax100_pt-bjet"]->fill(BottomJet_mom.pt(), event);
      }
      if ( b_jets400.size()==1 && bbar_jets400.size()==1 ){
        FourMomentum BottomJet_mom       =
        FourMomentum(b_jets400[0].E(),b_jets400[0].px(),b_jets400[0].py(),b_jets400[0].pz());
        FourMomentum AntiBottomJet_mom   =
        FourMomentum(bbar_jets400[0].E(),bbar_jets400[0].px(),
                     bbar_jets400[0].py(),bbar_jets400[0].pz());
        FourMomentum Top_mom             = BottomJet_mom + Wplus_mom;
        FourMomentum AntiTop_mom         = AntiBottomJet_mom + Wminus_mom;

        _histos["dmax400_m_Wplusb"]->fill( Top_mom.mass() ,event  );
        _histos["dmax400_m_Wminusbb"]->fill( AntiTop_mom.mass() ,event  );
        _histos["dmax400_pt-Wm"]->fill(Wplus_mom.pt(), event);
        _histos["dmax400_angle-Wm-bjet"]->fill(BottomJet_mom.angle(Wminus_mom) ,event);
        _histos["dmax400_angle-Wp-bjet"]->fill(BottomJet_mom.angle(Wplus_mom) ,event);
        _histos["dmax400_angle-bjet-bbjet"]->fill(BottomJet_mom.angle(AntiBottomJet_mom) ,event);
        _histos["dmax400_pt-bjet"]->fill(BottomJet_mom.pt(), event);
      }



      /* calculate not IR-safe "observables" */
      const FinalState& Bottoms     = applyProjection<IdentifiedFinalState>(event, "Bottom");
      const FinalState& AntiBottoms = applyProjection<IdentifiedFinalState>(event, "AntiBottom");
      const FinalState& Gluons = applyProjection<IdentifiedFinalState>(event, "Gluon");
      const ParticleVector bottoms_par     = Bottoms.particlesByPt();
      const ParticleVector antibottoms_par = AntiBottoms.particlesByPt();
      const ParticleVector gluons_par = Gluons.particlesByPt();
      FourMomentum Bottom_mom     = bottoms_par[0].momentum();
      FourMomentum AntiBottom_mom = antibottoms_par[0].momentum();
      FourMomentum Top_mom        = Bottom_mom + Wplus_mom;
      FourMomentum AntiTop_mom    = AntiBottom_mom + Wminus_mom;

      _histos["m_Wplusb"]->fill(Top_mom.mass(), event);
      _histos["m_Wminusbb"]->fill(AntiTop_mom.mass(), event);

      const int number_emissions = gluons_par.size();
      _histos["pt-b"]->fill(AntiBottom_mom.pt(), event);
      _h_count_emissions->fill( number_emissions==0 ? 0.25 : 0.75, 1. );

      if(number_emissions>=1){
        FourMomentum Gluon_mom = gluons_par[0].momentum();

        _histos["E-gluon"]->fill(Gluon_mom.E(), event);
        _histos["pt-gluon"]->fill(Gluon_mom.pt(), event);
        _histos["rapidity-gluon"]->fill(Gluon_mom.rapidity(), event);
        _histos["phi-gluon"]->fill(Gluon_mom.phi(), event);
        _histos["theta-gluon"]->fill(Gluon_mom.theta(), event);
        _histos["angle-b-gluon"]->fill(Gluon_mom.angle(Bottom_mom) ,event);
        _histos["angle-Wm-gluon"]->fill(Gluon_mom.angle(Wminus_mom) ,event);
      }
        _histos["theta-b"]->fill(Bottom_mom.theta(), event);
        _histos["theta-bb"]->fill(AntiBottom_mom.theta(), event);
    }


    /// Finalize
    void finalize() {

       double NORM=crossSection()/sumOfWeights();

       for (std::map<string,NLOHisto1DPtr>::iterator hit=_histos.begin(); hit!=_histos.end();hit++) {
         hit->second->finalize();
         scale(hit->second,NORM);
       }
    }

  private:
	  std::map<string,NLOHisto1DPtr> _histos;
      Histo1DPtr _h_count_emissions;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(m_wbjet);


}
