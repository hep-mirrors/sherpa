#ifndef ATOOLS_Phys_FormFactor_EMnucleon_Dipole_H
#define ATOOLS_Phys_FormFactor_EMnucleon_Dipole_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"

namespace ATOOLS {

  class FormFactor_EMnucleon_Dipole : public FormFactor_Model_Interface
  {
  public:
    FormFactor_EMnucleon_Dipole(incomingboson::code boson, incomingnucleon::code nucleon, double gA, double sin2thetaW);
    virtual ~FormFactor_EMnucleon_Dipole();

    void RegisterDefaultsAxial();
    void RegisterDefaultsDipole();

    double Q2_check(const double &Q2); 

    // Helper functions for Dipole parametrisation
    double tau_eval(const double &Q2, const double &mass);
    double Dipole_func(const double &Q2, const double &Lambda2, const double &norm);

    // Individual form factor functions (override virtual functions from interface)
    virtual double F1p(const double &Q2) override;
    virtual double F2p(const double &Q2) override;
    virtual double F1n(const double &Q2) override;
    virtual double F2n(const double &Q2) override;
    virtual double F1W(const double &Q2) override;
    virtual double F2W(const double &Q2) override;
    virtual double FAW(const double &Q2) override;
    virtual double FPW(const double &Q2) override;
    virtual double F1Zp(const double &Q2) override;
    virtual double F2Zp(const double &Q2) override;
    virtual double FAZp(const double &Q2) override;
    virtual double FPZp(const double &Q2) override;
    virtual double F1Zn(const double &Q2) override;
    virtual double F2Zn(const double &Q2) override;
    virtual double FAZn(const double &Q2) override;
    virtual double FPZn(const double &Q2) override;

    // Dipole parametrisation parameters
    double m_massp;   
    double m_mup;     
    double m_Lambda2, m_norm;
    double m_massn;   
    double m_mun;     

    // Axial form factor parameters
    double m_massA;   
    double m_fA;      
    double m_masspi;
    
    // Store boson/nucleon type and global params
    incomingboson::code m_boson_type;
    incomingnucleon::code m_nucleon_type;
    double m_gA;
    double m_sin2thetaW;
    ffmodel::code m_formfactor_model;
  };

}  // namespace ATOOLS

#endif
