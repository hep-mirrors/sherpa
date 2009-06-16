#ifndef Analysis_Triggers_Midpoint_Cone_H
#define Analysis_Triggers_Midpoint_Cone_H

#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Histogram.H"
#include <list>
#include <vector>

namespace ANALYSIS {
  class Primitive_Analysis;

  class Midpoint_Cone : public Jet_Algorithm_Base {

    struct _Vector {
      _Vector() 
	: p(0.0,0.0,0.0,0.0), pt(0.0), eta(0.0), phi(0.0) {}
      
      _Vector(const ATOOLS::Vec4D& q) 
	: p(q), pt(q.PPerp()), eta(q.Y()), phi(q.Phi()) {}
      
      ATOOLS::Vec4D p;
      double pt, eta, phi;
    };

    struct _Proto {
      _Vector psum;
      std::list<unsigned int> pars;
    };

    void _M_do_clustering(double, double);
    void _M_check_trial_cone(unsigned int, unsigned int *, double);
    unsigned int _M_check_shared_towers(const _Proto&, const _Proto&);
    void _M_split_merge(double, std::list<_Proto>::iterator, std::list<_Proto>::iterator);
    
    void _M_assadd(_Vector&, const _Vector&) const;
    void _M_iterate_cone(double, const _Vector&);
    void _M_was_it_already_found(const _Proto *);

    //       private data members
    //   recombination scheme
    int _M_reco, _M_midpoint;
    
    //   momenta of the particles
    std::vector<_Vector> _M_pp;
    
    //   list of the protojets
    std::list<_Proto> _M_proto; 
    
    //   the list of the overlapping partons
    std::list<unsigned int> _M_overlap; 
    
    //   temporary variable to store the jet momenta
    ATOOLS::Particle_List       * p_jets;
    std::vector<double> * p_kts;
    
    static double _S_dphi(double);

    double m_f;

    void AddToKtlist(double );
    void AddToJetlist(const ATOOLS::Vec4D &, bool);

  public:
    Midpoint_Cone(ATOOLS::Particle_Qualifier_Base * const qualifier,
		  int reco = 0, double f=.5);
    ~Midpoint_Cone();

    bool ConstructJets(const ATOOLS::Particle_List * pl = 0,ATOOLS::Particle_List * jets = 0,
		       std::vector<double> * kt2 = 0,double R=.7);
    int    NumberOfJets();

  private:
    struct pT_sort {
      bool operator()(const _Proto& p1, const _Proto& p2) const {
	return p1.psum.pt > p2.psum.pt;
      }
    };
  };



  inline int    Midpoint_Cone::NumberOfJets() { return p_jets->size(); }

  //
}

#endif
