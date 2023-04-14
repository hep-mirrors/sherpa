#ifndef YFS_RealVirtual_Virtual_H
#define YFS_RealVirtual_Virtual_H
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Color_Correlated_ME2.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "YFS/Main/YFS_Base.H"


// namespace PHASIC{
//   class Color_Correlated_ME2;
// }

namespace YFS{

    class RealVirtual : public YFS_Base
    {
    public:
        RealVirtual(const PHASIC::Process_Info& pi);
        ~RealVirtual();
        double Calc(const ATOOLS::Vec4D_Vector momenta, double b);    
        PHASIC::Virtual_ME2_Base     *p_loop_me;
        PHASIC::Color_Correlated_ME2* p_corr_me;
        double Calc_V(const ATOOLS::Vec4D_Vector& p,
                    const double B,
                    const double mur);
        MODEL::Coupling_Map m_cpls;
        PHASIC::Scale_Setter_Base *p_scale;
        PHASIC::Tree_ME2_Base* p_born_me;
        double m_sym;

    };
}



#endif
