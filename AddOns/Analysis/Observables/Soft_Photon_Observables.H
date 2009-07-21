#ifndef Analysis_Observables_Soft_Photon_Observables_H
#define Analysis_Observables_Soft_Photon_Observables_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Soft_Photon_Observable_Base : public Primitive_Observable_Base {
  protected:
    bool f_special;
    std::vector<ATOOLS::Flavour> m_flavs;
    virtual void Evaluate(const ATOOLS::Particle_Vector&,
                          double=1.,double=1.) = 0;
    virtual void Evaluate(double,double=1.,double=1.);
  public:
    Soft_Photon_Observable_Base(const std::vector<ATOOLS::Flavour> &,
                                int type, double xmin, double xmax, int nbins,
                                const std::string& listname,
                                const std::string& name);
    virtual void Evaluate(int nout, const ATOOLS::Vec4D *,
                          const ATOOLS::Flavour *, double=1., double=1.);
    virtual void Evaluate(const ATOOLS::Particle_List & plist,
                          double weight, double ncount);
    virtual void Evaluate(const ATOOLS::Blob_List & blobs, double weight,
                          double ncount);
  };

  //==========================================================================

  class Soft_Photon_Energy : public Soft_Photon_Observable_Base {
    void Evaluate(const ATOOLS::Particle_Vector&, double=1., double=1.);
  public:
    Soft_Photon_Energy(const std::vector<ATOOLS::Flavour>&,
                       int type, double xmin, double xmax, int nbins,
                       const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

  class Soft_Photon_Angle : public Soft_Photon_Observable_Base {
    void Evaluate(const ATOOLS::Particle_Vector&, double=1., double=1.);
  public:
    Soft_Photon_Angle(const std::vector<ATOOLS::Flavour>&,
                      int type, double xmin, double xmax, int nbins,
                      const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

}
#endif
