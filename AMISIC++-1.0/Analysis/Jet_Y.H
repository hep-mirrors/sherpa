#ifndef Jet_Y_H
#define Jet_Y_H

#include "Primitive_Observable_Base.H"

namespace AMISIC {

  class Jet_Y: public ATOOLS::Primitive_Observable_Base {
  private:

    std::vector<ATOOLS::Histogram*> *p_histogram;

    std::vector<int>    *p_jets;
    std::vector<double> *p_ys;

    unsigned int m_njets, m_mode;

    void SortJetPT(std::vector<ATOOLS::Vec4D> &jetmomenta);

  public:

    // constructor
    Jet_Y(int _m_type,double _m_min,double _m_max,int _m_nbins,int _m_njets);

    // destructor
    ~Jet_Y();

    // setters
    void SetMode(unsigned int _m_mode)                     { m_mode=_m_mode;           }

    // getters
    unsigned int           Mode()      { return m_mode;      }

    // member functions
    void Evaluate(const ATOOLS::Particle_List &particles,double weight);
    void Evaluate(const ATOOLS::Blob_List &blobs,double weight);

    void EndEvaluation();
    void Output(std::string pathname);

  };// end of class Jet_Y

}// end of namespace AMISIC

#endif
