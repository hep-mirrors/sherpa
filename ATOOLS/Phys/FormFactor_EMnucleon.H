#ifndef ATOOLS_Phys_FormFactor_EMnucleon_H
#define ATOOLS_Phys_FormFactor_EMnucleon_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include <memory>

namespace ATOOLS {

  // Forward declaration of the generic model interface
  class FormFactor_Model_Interface;

  struct ffmodel
  {
    enum code // using code is "bad coding apparently", fix this later 
    {
      off = 0,
      kelly = 1,
      dipole = 2
    };
  };
  // Stream operators for ffmodel
  std::istream &operator>>(std::istream &str, ffmodel::code &model);
  std::ostream &operator<<(std::ostream &str, const  ffmodel::code &model);

  struct incomingboson
  {
    enum code
    {
      off = 0,
      photon = 1,
      W = 2,
      Z = 3
    };
  };

  struct incomingnucleon
  { 
    enum code
    {
      off = 0,
      proton = 1,
      neutron = 2
    };
  };
  // Stream operators for incomingnucleon
  std::istream &operator>>(std::istream &str, incomingnucleon::code &nucleon);
  std::ostream &operator<<(std::ostream &str, const incomingnucleon::code &nucleon);

  // Struct to hold boson-nucleon interaction type for printing
  struct BosonNucleonType {
    incomingboson::code boson;
    incomingnucleon::code nucleon;
    
    BosonNucleonType(incomingboson::code b = incomingboson::off, 
                     incomingnucleon::code n = incomingnucleon::off)
        : boson(b), nucleon(n) {}
  };
  // Stream operators for BosonNucleonType
  std::ostream &operator<<(std::ostream &str, const BosonNucleonType &type);
  std::istream &operator>>(std::istream &str, BosonNucleonType &type);

  // Struct to hold nucleon form factors
  struct NucleonFormFactors {
    double F1;  // Dirac (vector) form factor
    double F2;  // Pauli (anomalous magnetic) form factor
    double FA;  // Axial form factor
    double FP;  // Pseudoscalar form factor

    // J_mu = Vector_mu - Axial_mu 
    //      = F1*gamma_mu + F2*i*sigma_mu_nu*q_nu/(2M) - FA*gamma_mu*gamma5   - FP*q_mu*gamma5/(2M)

    // default values
    NucleonFormFactors(double f1 = 0.0, double f2 = 0.0, double fa = 0.0, double fp = 0.0)
        : F1(f1), F2(f2), FA(fa), FP(fp) {}
  };

  // Abstract interface that all form factor models must implement
  class FormFactor_Model_Interface {
  public:
    virtual ~FormFactor_Model_Interface() {}
    
    // Pure virtual functions, all models must implement these otherwise =0 gives error
    virtual double F1p(const double &q2) = 0;
    virtual double F2p(const double &q2) = 0;
    virtual double F1n(const double &q2) = 0;
    virtual double F2n(const double &q2) = 0;
    virtual double F1W(const double &q2) = 0;
    virtual double F2W(const double &q2) = 0;
    virtual double FAW(const double &q2) = 0;
    virtual double FPW(const double &q2) = 0;
    virtual double F1Zp(const double &q2) = 0;
    virtual double F2Zp(const double &q2) = 0;
    virtual double FAZp(const double &q2) = 0;
    virtual double FPZp(const double &q2) = 0;
    virtual double F1Zn(const double &q2) = 0;
    virtual double F2Zn(const double &q2) = 0;
    virtual double FAZn(const double &q2) = 0;
    virtual double FPZn(const double &q2) = 0;
  };

  class FormFactor_EMnucleon //base //virtual means you can override, double means no override , =0 you have to fill them 
  {
  protected:
    // Pointer to model implementation (no mention of specific model)
    std::unique_ptr<FormFactor_Model_Interface> p_model_impl;

  public:
    FormFactor_EMnucleon(incomingboson::code boson, incomingnucleon::code nucleon);
    virtual ~FormFactor_EMnucleon();

    void RegisterDefaultsOff();
    
    // Main interface: Get form factors using stored boson and nucleon types
    NucleonFormFactors GetFormFactors(const double &q2);
    inline bool On() { return m_formfactor_model; }

    // Individual form factor functions 
    // (point to implementations in derived classes)
    double F1p(const double &q2);
    double F2p(const double &q2);
    double F1n(const double &q2);
    double F2n(const double &q2);
    double F1W(const double &q2);
    double F2W(const double &q2);
    double FAW(const double &q2);
    double FPW(const double &q2);
    double F1Zp(const double &q2);
    double F2Zp(const double &q2);
    double FAZp(const double &q2);
    double FPZp(const double &q2);
    double F1Zn(const double &q2);
    double F2Zn(const double &q2);
    double FAZn(const double &q2);
    double FPZn(const double &q2);

    // Functions returning form factor structs
    NucleonFormFactors Photon_Proton(const double &q2);
    NucleonFormFactors Photon_Neutron(const double &q2);
    NucleonFormFactors W_Boson(const double &q2);
    NucleonFormFactors Z_Proton(const double &q2);
    NucleonFormFactors Z_Neutron(const double &q2);

    // global parameters
    bool m_PCAC;
    double m_sin2thetaW;
    double m_cos2thetaW;
    double m_gA;
    double m_g;
    double m_masspi;

    // Off form factor model parameters
    double m_F1p, m_F2p;
    double m_F1n, m_F2n;
    double m_F1W, m_F2W, m_FAW, m_FPW;
    double m_F1Zp, m_F2Zp, m_FAZp, m_FPZp;
    double m_F1Zn, m_F2Zn, m_FAZn, m_FPZn;

    ffmodel::code m_formfactor_model;
    incomingboson::code m_boson_type;
    incomingnucleon::code m_nucleon_type;
  };
  extern FormFactor_EMnucleon *nucleonformfactor;
}  // namespace ATOOLS

#endif
