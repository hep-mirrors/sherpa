// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <map>

namespace Rivet {


  /// @brief A quick and dirty MC analysis to check instanton mass and rapidity
  // plus the number of outgoing gluons in different mass bins
  class MC_INSTANTON : public Analysis {
  private:
    std::map<std::string, Histo1DPtr> m_histos;
  public:
    MC_INSTANTON() : Analysis("MC_INSTANTON") {
    }
    void init() {
      declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 1000*MeV), "FS");
      book(m_histos["Mass"],        "Mass",        50, 50.0, 550.0);
      book(m_histos["Mass_Shower"], "Mass_Shower", 50, 50.0, 550.0);
      book(m_histos["Rapidity"],    "Rapidity",    10, -5.0,   5.0);
      book(m_histos["N_50"],        "N_50",        50,  5.5,  55.5);
      book(m_histos["N_100"],       "N_100",       50,  5.5,  55.5);
      book(m_histos["N_200"],       "N_200",       50,  5.5,  55.5);
      book(m_histos["N_b_50"],      "N_b_50",      10, -0.5,   9.5);
      book(m_histos["N_b_100"],     "N_b_100",     10, -0.5,   9.5);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const GenEvent   * genevent = event.genEvent();
      const double     & weight   = event.weight();
      ConstGenVertexPtr signal      = genevent->signal_process_vertex();
      ConstGenParticlePtr instanton = (*signal->particles_out_const_begin());
      ConstGenVertexPtr shower      = instanton->end_vertex();
      FourMomentum  p_instanton     = instanton->momentum();
      FourMomentum  p_total(0.,0.,0.,0.);
      size_t N_b(0);
      for (ConstGenParticlePtr pout:
	     HepMCUtils::particles(shower, Relatives::CHILDREN)) {
	p_total += pout->momentum();
	if (pout->pdg_id()==5 || pout->pdg_id()==-5) N_b++;
      }
      double        M_instanton   = p_instanton.mass();
      double        M_total       = p_total.mass();
      double        Y_instanton   = p_instanton.rapidity();
      // subtract 2 incoming partons
      size_t        N_partons     = shower->particles_out_size()-2;
      m_histos["Mass"]->fill(M_instanton,weight);
      m_histos["Mass_Shower"]->fill(M_total,weight);
      m_histos["Rapidity"]->fill(Y_instanton,weight);

      if (M_instanton>50. && M_instanton<100.)  {
	m_histos["N_50"]->fill(N_partons,weight);
	m_histos["N_b_50"]->fill(N_b,weight);
      }
      if (M_instanton>100. && M_instanton<200.) {
	m_histos["N_100"]->fill(N_partons,weight);
	m_histos["N_b_100"]->fill(N_b,weight);
      }
      if (M_instanton>200.) m_histos["N_200"]->fill(N_partons,weight);
      /*
	// some control output below.  cross-checked with sherpa event record.  
	std::cout<<"Event with weight = "<<weight<<"\n"
	<<"  signal    = "<<(*signal) <<"\n"
	<<"  instanton = "<<(*instanton)<<"\n"
	<<"  mass      = "<<M_instanton<<", "
	<<N_partons<<" outgoing particles.\n";
	for (std::vector<HepMC::GenParticle*>::const_iterator pout=
	shower->particles_out_const_begin();
	pout!=shower->particles_out_const_end();pout++) {
	std::cout<<"   --> "<<(**pout)<<"\n";
	}
      */
    }
    
    void finalize() {
      for (std::map<std::string, Histo1DPtr>::iterator hit=m_histos.begin();
	   hit!=m_histos.end();hit++) {
	if (hit->first==string("Mass") ||
	    hit->first==string("Rapidity")) {
	  // norm to cross section
	  scale(hit->second, crossSection()/picobarn/sumOfWeights());
	}
	else {
	  // normalize to unity
	  normalize(hit->second); 
	}
      }
    }
  };
  DECLARE_RIVET_PLUGIN(MC_INSTANTON);
};
