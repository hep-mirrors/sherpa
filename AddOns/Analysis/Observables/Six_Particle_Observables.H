#ifndef Analysis_Observables_Six_Particle_Observables_H
#define Analysis_Observables_Six_Particle_Observables_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Six_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    bool f_special;
    std::vector<ATOOLS::Flavour> m_flavs;
    void virtual Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                          const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                          const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                          double=1.,double=1.) = 0; 
    void virtual Evaluate(double,double=1.,double=1.); 
  public:
    Six_Particle_Observable_Base(const std::vector<ATOOLS::Flavour> &,
          int,double,double,int,
          const std::string& listname,
          const std::string& name);
    void virtual Evaluate(int,const ATOOLS::Vec4D *,
        const ATOOLS::Flavour *,double=1.,double=1.);
    void virtual Evaluate(const ATOOLS::Particle_List &,
        double=1.,double=1.);
  };

  //==========================================================================

  class Six_Particle_PT : public Six_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  double=1.,double=1.);
  public:
    Six_Particle_PT(const std::vector<ATOOLS::Flavour>&,
                    int,double,double,int,
                    const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

  class Six_Particle_ET : public Six_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  double=1.,double=1.);
  public:
    Six_Particle_ET(const std::vector<ATOOLS::Flavour>&,
                    int,double,double,int,
                    const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

  class Six_Particle_Mass : public Six_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
                  double=1.,double=1.);
  public:
    Six_Particle_Mass(const std::vector<ATOOLS::Flavour>&,
                      int,double,double,int,
                      const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

  class Two_Partontriplett_DeltaPhi : public Six_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  double=1.,double=1.);
  public:
    Two_Partontriplett_DeltaPhi(const std::vector<ATOOLS::Flavour>&,
                                int, double, double, int,
                                const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

  class Two_Partontriplett_DeltaEta : public Six_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  double=1.,double=1.);
  public:
    Two_Partontriplett_DeltaEta(const std::vector<ATOOLS::Flavour>&,
                                int, double, double, int,
                                const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

  class Two_Partontriplett_DR : public Six_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
                  double=1.,double=1.);
  public:
    Two_Partontriplett_DR(const std::vector<ATOOLS::Flavour>&,
                          int, double, double, int,
                          const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };


  // ======================================================================

//   class Di_Mass : public Primitive_Observable_Base {
//   public:
//     Di_Mass(unsigned int type,double xmin,double xmax,int nbins,
//      const std::string & =std::string("Jets"));
//     void Evaluate(const ATOOLS::Blob_List & blobs,double value, int);
//     Primitive_Observable_Base * Copy() const ;
//   };

  // =====================================================================

//   class Six_Particle_PlaneAngleCMS : public Six_Particle_Observable_Base {
//     void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &, 
//                   const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
//                   double=1.,int=1); 
//   public:
//     Six_Particle_PlaneAngleCMS(const std::vector<ATOOLS::Flavour>&,
//                              int,double,double,int,
//                              const std::string& listname);//=std::string("NRAngle:"));
//     Primitive_Observable_Base* Copy() const;
//   };

  // =====================================================================

//   class Six_Particle_EnergyCMS : public Six_Particle_Observable_Base {
//     void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &, 
//                   const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
//                   double=1.,int=1); 
//   public:
//     Six_Particle_EnergyCMS(const std::vector<ATOOLS::Flavour>&,
//                              int,double,double,int,
//                              const std::string& listname);//=std::string("NRAngle:"));
//     Primitive_Observable_Base* Copy() const;
//   };

}

#endif
