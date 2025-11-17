#ifndef ATOOLS_Phys_FormFactor_EMnucleon_H
#define ATOOLS_Phys_FormFactor_EMnucleon_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  struct ffmodel
  {
    enum code // fix this later 
    {
      off = 0,
      kelly = 1,
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
    
    // default values
    NucleonFormFactors(double f1 = 0.0, double f2 = 0.0, double fa = 0.0, double fp = 0.0)
        : F1(f1), F2(f2), FA(fa), FP(fp) {}
  };

  class FormFactor_EMnucleon //base //virtual means you can override, double means no override , =0 you have to fill them 
  {
  public:
    FormFactor_EMnucleon(incomingboson::code boson, incomingnucleon::code nucleon);
    ~FormFactor_EMnucleon();

    void RegisterDefaultsOff();
    void RegisterDefaultsAxial();
    void RegisterDefaultsKelly();
    
    // Main interface: Get form factors using stored boson and nucleon types
    NucleonFormFactors GetFormFactors(const double &q2);
    inline bool On() { return m_formfactor_model; }

    double Q2_check(const double &Q2); 

    // Helper functions for kelly parametrisation
    double tau_eval(const double &Q2, const double &mass);
    double Kelly_func(const double &Q2, const double &a0, const double &a1,
                      const double &b1, const double &b2, const double &b3);

    // Individual form factor functions 
    // (internally changes based on m_formfactor_model)
    double F1p(const double &Q2);
    double F2p(const double &Q2);
    double F1n(const double &Q2);
    double F2n(const double &Q2);
    double F1W(const double &Q2);
    double F2W(const double &Q2);
    double FAW(const double &Q2);
    double FPW(const double &Q2);
    double F1Zp(const double &Q2);
    double F2Zp(const double &Q2);
    double FAZp(const double &Q2);
    double FPZp(const double &Q2);
    double F1Zn(const double &Q2);
    double F2Zn(const double &Q2);
    double FAZn(const double &Q2);
    double FPZn(const double &Q2);

    // Functions returning form factor structs
    NucleonFormFactors Photon_Proton(const double &Q2);
    NucleonFormFactors Photon_Neutron(const double &Q2);
    NucleonFormFactors W_Boson(const double &Q2);
    NucleonFormFactors Z_Proton(const double &Q2);
    NucleonFormFactors Z_Neutron(const double &Q2);

    // global parameters
    double m_sin2thetaW;
    double m_gA;

    // Off form factor model parameters
    double m_F1p, m_F2p;
    double m_F1n, m_F2n;
    double m_F1W, m_F2W, m_FAW, m_FPW;
    double m_F1Zp, m_F2Zp, m_FAZp, m_FPZp;
    double m_F1Zn, m_F2Zn, m_FAZn, m_FPZn;

    // Kelly parametrisation parameters
    double m_massp;   
    double m_mup;     
    double m_a0pE, m_a1pE;
    double m_b1pE, m_b2pE, m_b3pE;
    double m_a0pM, m_a1pM;
    double m_b1pM, m_b2pM, m_b3pM;
    double m_massn;   
    double m_mun;     
    double m_An, m_Bn, m_Deltasq;
    double m_a0nE, m_a1nE;
    double m_b1nE, m_b2nE, m_b3nE;
    double m_a0nM, m_a1nM;
    double m_b1nM, m_b2nM, m_b3nM;

    // Axial form factor parameters
    double m_massA;   
    double m_fA;      
    double m_masspi;  

    ffmodel::code m_formfactor_model;
    incomingboson::code m_boson_type;
    incomingnucleon::code m_nucleon_type;
  };
  extern FormFactor_EMnucleon *nucleonformfactor;
}  // namespace ATOOLS

#endif
