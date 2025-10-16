#ifndef ATOOLS_Phys_FormFactor_H
#define ATOOLS_Phys_FormFactor_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"


namespace ATOOLS {

struct finalstate {
  enum code {
    off       = 0,
    pion       = 1,
    kplus       = 2
  };
};

std::istream &operator>>(std::istream &str, finalstate::code &sm);
std::ostream &operator<<(std::ostream &str, const finalstate::code &sm);

    class FormFactor
    {
    public:
    	FormFactor();
    	~FormFactor();
        
        void RegisterDefaultsPion();
        void RegisterDefaultsKaon();

        double Eval(const double &q2);

    	Complex ppi(const double &q2);
        Complex d(const double &m);
        Complex h(const double &q2);
        Complex dh(const double &q2);

        Complex f(const double &q2, const double &m, const double &g);
        Complex gamma(const double &q2, const double &m, const double &g);

        inline bool On() {return m_form_mode;}

        Complex b(const double &q2, const Flavour &fl);
        Complex BW_GS(const double &q2, const Flavour &fl);
        Complex Pion(const double &q2);
        Complex KPlus(const double &q2);
        double m_crho, m_crhop, m_crhopp, m_crhoppp;
        double m_comega,  m_comegap,  m_comegapp; 
        double m_cphi, m_cphip, m_cphipp;
        double m_prho, m_prhop, m_prhopp, m_prhoppp, m_pomega, m_pphi;
        double m_omega_mass, m_omega_g, m_omegap_mass, m_omegap_g, m_omegapp_mass, m_omegapp_g;
        double m_phi_mass, m_phi_g, m_phip_mass, m_phip_g;
        double  m_rho_mass, m_rho_g, m_rhop_mass, m_rhop_g, m_rhopp_mass, m_rhopp_g, m_rhoppp_mass, m_rhoppp_g;
        Flavour m_flv, m_rho, m_rhop, m_rhopp, m_rhoppp;
        Flavour m_omega, m_omegap, m_omegapp;
        Flavour m_phi, m_phip;

        Complex m_I;
        finalstate::code m_form_mode;	
    };
    extern FormFactor * pionform;
}

#endif
