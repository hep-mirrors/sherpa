#ifndef YFS_Virtual_Virtual_H
#define YFS_Virtual_Virtual_H
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

    class VirtualVirtual : public YFS_Base
    {
    public:
        VirtualVirtual(const PHASIC::Process_Info& pi);
        ~VirtualVirtual();
        double Calc(const ATOOLS::Vec4D_Vector momenta, double b);    
        std::unique_ptr<PHASIC::Virtual_ME2_Base> p_loop_me;
        double Calc_V(const ATOOLS::Vec4D_Vector& p,
                    const double B,
                    const double mur);
        MODEL::Coupling_Map m_cpls;
        double m_factor;
        //! Not owned; set by SetProc, may be null.
        PHASIC::Process_Base *p_vproc = nullptr;
        std::string m_gen;
        inline std::string ActiveGenName() const { return m_gen.empty()?"none":m_gen; }
        inline void SetProc(PHASIC::Process_Base *proc) {p_vproc=proc;}
        inline bool FailCut() {return m_failcut;}
    };
}



#endif
