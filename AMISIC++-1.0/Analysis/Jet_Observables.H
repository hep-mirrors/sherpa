#ifndef Jet_Observables_H
#define Jet_Observables_H

#include "Jet_Multiplicity.H"
#include "Jet_PT.H"
#include "Jet_Y.H"

namespace AMISIC {

  class Jet_Observables: public ATOOLS::Primitive_Observable_Base {
  private:

    ATOOLS::Jet_Algorithm *p_jetfinder;

    Jet_Multiplicity *p_jetmulti;
    Jet_PT           *p_jetpt;
    Jet_Y            *p_jety;

  public:

    // constructor
    Jet_Observables(int _m_typept,double _m_minpt,double _m_maxpt,int _m_nbinspt,
		    int _m_typey,double _m_miny,double _m_maxy,int _m_nbisy,
		    int _m_jets,int _m_maxmult,double _m_ycut,double _m_ptcut);

    // destructor
    ~Jet_Observables();
    
    // member functions
    void Evaluate(const ATOOLS::Particle_List &particles,double weight);
    void Evaluate(const ATOOLS::Blob_List &blobs,double weight);

    void EndEvaluation();
    void Output(std::string pathname); 

  };// end of class Jet_Observables

}// end of namespace AMISIC

#endif
