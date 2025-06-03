#ifndef ATOOLS_Phys_Pion_FormFactor_H
#define ATOOLS_Phys_Pion_FormFactor_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"


namespace ATOOLS {

    class Pion_FormFactor
    {
    public:
    	Pion_FormFactor();
    	~Pion_FormFactor();

        double Eval(const double &q2);

    	double ppi(const double &q2);
        double d(const double &m);
        double h(const double &q2);
        double dh(const double &q2);

        double f(const double &q2, const double &m, const double &g);
        double gamma(const double &q2, const double &m, const double &g);

        Complex b(const double &q2, const Flavour &fl);
        Complex GS(const double &q2, const Flavour &fl);

    	double m_crho, m_crhop, m_crhopp, m_crhoppp, m_comega, m_cphi;
        double m_prho, m_prhop, m_prhopp, m_prhoppp, m_pomega, m_pphi;
        Flavour m_flv, m_rho, m_rhop, m_rhopp, m_rhoppp, m_omega, m_phi;

        Complex m_I;
        int m_form_mode;	
    };
    extern Pion_FormFactor * pionform;
}

#endif
