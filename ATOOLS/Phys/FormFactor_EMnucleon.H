#ifndef ATOOLS_Phys_FormFactor_EMnucleon_H
#define ATOOLS_Phys_FormFactor_EMnucleon_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

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

  class FormFactor_EMnucleon
  {
  public:
    FormFactor_EMnucleon(incomingboson::code boson, incomingnucleon::code nucleon);
    ~FormFactor_EMnucleon();

    void RegisterDefaultsProton();
    void RegisterDefaultsNeutron();
    void RegisterDefaultsAxial();
    void RegisterDefaultsZ();
    
    // Main interface: Get form factors using stored boson and nucleon types
    NucleonFormFactors GetFormFactors(const double &q2);

    double Q2_eval(const double &q2);
    double tau_eval(const double &q2, const double &mass);

    double GEp(const double &q2);
    double GMp(const double &q2);
    double GEn(const double &q2);
    double GMn(const double &q2);

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

    // Proton parameters
    double m_massp;   // Proton mass
    double m_mup;     // Proton magnetic moment
    // Kelly parametrisation coefficients for proton electric form factor
    double m_a0pE, m_a1pE;
    double m_b1pE, m_b2pE, m_b3pE;
    // Kelly parametrisation coefficients for proton magnetic form factor
    double m_a0pM, m_a1pM;
    double m_b1pM, m_b2pM, m_b3pM;

    // Neutron parameters
    double m_massn;   // Neutron mass
    double m_mun;     // Neutron magnetic moment
    // Galster parametrisation for neutron
    double m_An, m_Bn, m_Deltasq;
    // Kelly parametrisation coefficients for neutron electric form factor
    double m_a0nE, m_a1nE;
    double m_b1nE, m_b2nE, m_b3nE;
    // Kelly parametrisation coefficients for neutron magnetic form factor
    double m_a0nM, m_a1nM;
    double m_b1nM, m_b2nM, m_b3nM;

    // Axial form factor parameters
    double m_massA;   // Axial mass
    double m_gA;      // Axial coupling constant
    double m_fA;      // Additional axial parameter
    double m_masspi;  // Pion mass

    // Weak mixing angle
    double m_sin2thetaW;

    incomingboson::code m_boson_type;
    incomingnucleon::code m_nucleon_type;
  };
extern FormFactor_EMnucleon *nucleonformfactor;
}  // namespace ATOOLS

#endif
