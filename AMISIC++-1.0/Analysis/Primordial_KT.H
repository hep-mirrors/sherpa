#ifndef Primordial_KT_H
#define Primordial_KT_H

#include "Primitive_Observable_Base.H"

namespace AMISIC {

  class Primordial_KT: public ATOOLS::Primitive_Observable_Base {
  public:

    // constructor
    Primordial_KT(int _m_type,double _m_min,double _m_max,int _m_nbins);

    // destructor
    ~Primordial_KT();

    // member functions
    void Evaluate(const ATOOLS::Particle_List &particles,double weight);
    void Evaluate(const ATOOLS::Blob_List &blobs,double weight);

    void EndEvaluation();
    void Output(std::string pathname);

  };// end of class Primordial_KT

}// end of namespace AMISIC

#endif
