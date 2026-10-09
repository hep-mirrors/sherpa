 #ifndef YFS_Virtual_H
#define YFS_Virtual_H
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Color_Correlated_ME2.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/NLO/Massification.H"
#include <memory>

// namespace PHASIC{
//   class Color_Correlated_ME2;
// }

namespace YFS{

    class Virtual : public YFS_Base
    {
    public:
        Virtual(const PHASIC::Process_Info& pi);
        ~Virtual();
        double Calc(const ATOOLS::Vec4D_Vector momenta, double b);    
        std::unique_ptr<PHASIC::Virtual_ME2_Base> p_loop_me;
        //! YFS: VIRTUAL_MASSIFICATION (null when off)
        std::unique_ptr<Massification> p_massif;
        //! the provider's 1/eps coefficient (ME_E1 units) of the last call,
        //! the massive one under massification
        double E1() const;
        //! the 1/eps^2 coefficient, 0 under massification (checked per call)
        double E2() const;
        //! YFS: EW_HIGHER_ORDER griffin - GRIFFIN's NNLO - NLO at the Born
        //! point (null otherwise), added by NLO_Base::CalculateVirtual.
        std::unique_ptr<PHASIC::Virtual_ME2_Base> p_ewho_me;
        inline bool HasEWHigherOrder() const { return p_ewho_me != nullptr; }
        double EWHigherOrder(const ATOOLS::Vec4D_Vector &momenta);
        //! peval (optional): the frame the loop is evaluated in (NLO cuts
        //! and the running-alpha scale still use p). See YFS:
        //! VIRTUAL_CANONICAL_FRAME in NLO_Base::CalculateVirtual.
        double Calc_V(const ATOOLS::Vec4D_Vector& p,
                    const double B,
                    const double mur,
                    const ATOOLS::Vec4D_Vector *peval = nullptr);
        double CalcInFrame(const ATOOLS::Vec4D_Vector &momenta,
                           const ATOOLS::Vec4D_Vector &peval, double b);
        MODEL::Coupling_Map m_cpls;
        double m_factor;
        //! Not owned; set by SetProc, may be null.
        PHASIC::Process_Base *p_vproc = nullptr;
        inline void SetProc(PHASIC::Process_Base *proc) {p_vproc=proc;}
        inline bool FailCut() {return m_failcut;}
    };
}



#endif
