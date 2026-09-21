#ifndef YFS_RealVirtual_Virtual_H
#define YFS_RealVirtual_Virtual_H
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Color_Correlated_ME2.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "YFS/Main/YFS_Base.H"
#include <memory>


// namespace PHASIC{
//   class Color_Correlated_ME2;
// }

namespace YFS{

    class RealVirtual : public YFS_Base
    {
    public:
        RealVirtual(const PHASIC::Process_Info& pi);
        ~RealVirtual();
        double Calc(const ATOOLS::Vec4D_Vector &momenta, double b);    
        std::unique_ptr<PHASIC::Virtual_ME2_Base> p_loop_me;
        double Calc_V(const ATOOLS::Vec4D_Vector& p,
                    const double B,
                    const double mur);
        MODEL::Coupling_Map m_cpls;
        //! Not owned; set by SetProc, may be null.
        PHASIC::Process_Base *p_rvproc = nullptr;
        ATOOLS::Flavour_Vector m_flavs;
        double m_sym, m_factor;
        //! RV_Generator as configured. Unlike Real/RealReal there is no
        //! internal-vs-external choice here: p_loop_me (the value this class
        //! actually returns) is always built from RV_Generator directly, and
        //! p_rvproc is used only for Trigger()/cuts, not for the value -- so
        //! it must NOT be used as a generator label.
        std::string m_gen;
        inline std::string ActiveGenName() const { return m_gen.empty()?"none":m_gen; }
        inline void SetProc(PHASIC::Process_Base *realproc) {p_rvproc=realproc;}
        inline bool FailCut() {return m_failcut;}
    };
}



#endif
