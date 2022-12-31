#ifndef Virtual_Amp_H
#define Virtual_Amp_H

#include <complex>
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ggvvamp.h"
#include "ggvvampprec.h"
#include "Precise_Amp.H"
#include "Grids.H"

namespace DiBosonsNLO{

class Virtual_Amp: public PHASIC::Virtual_ME2_Base{

    const std::complex<double> I = std::complex<double>(0.,1.); // it is needed to define purely imaginary obj

    // Members needed for off-shell amplitudes squaring
    ATOOLS::Vec4D_Vector m_P;
    std::complex<double> M[2][4][4];
    std::complex<double> kin_factors[2][9];
    static const int m_nHel = 2;
    static const int m_nCoeff = 9;
    static const int m_nLoop = 2;
    static const int m_nreim = 2;
    int m_hel_independent_piece;
    //std::vector<int> m_flavs;
    int m_current_1_type[3];
    int m_current_2_type[3];
    int m_process_type[4];
    // This object is used to evaluate the 4 lepton final state form factors
    // using the precision check implemented by default in ggvvamp
    ggvvampprec amp;
    double m_E[m_nHel][m_nCoeff][m_nLoop][m_nreim];

    Precision::Precise_Amp F;
    int m_process_type_on_shell;
    double m_on_shell_factor;

    Grids::Grids interp_amp;
    bool m_use_grid;

    // General parameters
    //TODO: they might be moved into a Interface file/class
    ATOOLS::Flavour Z;
    ATOOLS::Flavour W;
    const double m_MZ = Z.Mass();
    const double m_WZ = Z.Width();
    const double m_MW = W.Mass();
    const double m_WW = W.Width();
    std::complex<double> m_cos2w = (m_MW*m_MW - I*m_WW*m_MW)/(m_MZ*m_MZ - I*m_WZ*m_MZ);
    std::complex<double> m_sin2w = 1.-m_cos2w;


    // Default values, are lated overwritten with Sherpa's values
    double alphas=0.118;
    double alphaqed=1./128.;


    double m_symfac;
    double m_beta0;
    bool m_gluon_induced, m_off_shell;



    std::map<int, PHASIC::Process_Info> m_procmap;
    ATOOLS::Flavour_Vector m_flavs;


    public:

    Virtual_Amp(const PHASIC::Process_Info& pi,
                const ATOOLS::Flavour_Vector& flavs,
                const bool off_shell, const bool use_grid);
    ~Virtual_Amp();


    // Methods for the off-shell case
    std::complex<double> calc_hel_amp(const int hel, const int nLoop);
    std::complex<double> hm(char hel, int i, int j);
    void get_form_factors();
    void get_interp_form_factors();
    inline ATOOLS::Vec4D hel_to_mom(int i) {return m_P[i];}
    void kin_factors_calc(const char hs[], int i0, int i1, int i2, int i3, int i4, int i5);
    void calc_hel_amp();
    double hel_squared(int i, int j);
    double born_hel_squared(int i, int j);
    std::complex<double> Denominator(const ATOOLS::Vec4D& Q, double M, double W);
    void Get_Current_Types();
    void Select_Processes();
    std::complex<double> Current_Factor(int helicity);
    double Couplings();
    double Sym_Factor();
    std::complex<double> LJy(const ATOOLS::Vec4D& Q, const int which_boson);
    std::complex<double> RJy(const ATOOLS::Vec4D& Q, const int which_boson);
    std::complex<double> LJZ(const ATOOLS::Vec4D& Q, const int which_boson);
    std::complex<double> RJZ(const ATOOLS::Vec4D& Q, const int which_boson);
    std::complex<double> LJW(const ATOOLS::Vec4D& Q, const int which_boson);
    double Calc_OffShell(const ATOOLS::Vec4D_Vector& P, const double& Born);


    double Interp_OffShell(const ATOOLS::Vec4D_Vector& P, const double& Born);


    // Methods for the on-shell case
    double Calc_OnShell(const ATOOLS::Vec4D_Vector& P, const double& Born);
    void Form_Factors_Computation(const ATOOLS::Vec4D_Vector& P);
    int Process_Type();
    double EWCouplings(int process_type);


    // General members needed for the Virtual_ME2_Base object
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& p);
    void Calc(const ATOOLS::Vec4D_Vector& P);
    void Calc(const ATOOLS::Vec4D_Vector& P, const double& Born);
    inline bool SetColours(const ATOOLS::Vec4D_Vector& momenta){return true;}
    bool IsMappableTo(const PHASIC::Process_Info& pi);

};
}


#endif
