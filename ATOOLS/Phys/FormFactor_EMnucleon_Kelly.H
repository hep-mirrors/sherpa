#ifndef ATOOLS_Phys_FormFactor_EMnucleon_Kelly_H
#define ATOOLS_Phys_FormFactor_EMnucleon_Kelly_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"

namespace ATOOLS {

  class FormFactor_EMnucleon_Kelly : public FormFactor_Model_Interface
  {
  public:
    FormFactor_EMnucleon_Kelly(incomingboson::code boson, incomingnucleon::code nucleon, double gA, double sin2thetaW);
    virtual ~FormFactor_EMnucleon_Kelly();

    void RegisterDefaultsAxial();
    void RegisterDefaultsKelly();

    double Q2_check(const double &Q2); 

    // Helper functions for kelly parametrisation
    double tau_eval(const double &Q2, const double &mass);
    double Kelly_func(const double &Q2, const double &a0, const double &a1,
                      const double &b1, const double &b2, const double &b3);

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
    
    // Store boson/nucleon type and global params
    incomingboson::code m_boson_type;
    incomingnucleon::code m_nucleon_type;
    double m_gA;
    double m_sin2thetaW;
    ffmodel::code m_formfactor_model;
  };

}  // namespace ATOOLS

#endif
