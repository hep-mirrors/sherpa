#ifndef YFS_RealVirtual_Virtual_H
#define YFS_RealVirtual_Virtual_H
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

    class RealVirtual : public YFS_Base
    {
    public:
        RealVirtual(const PHASIC::Process_Info& pi);
        ~RealVirtual();
        double Calc(const ATOOLS::Vec4D_Vector &momenta, double b);    
        std::unique_ptr<PHASIC::Virtual_ME2_Base> p_loop_me;
        //! YFS: VIRTUAL_MASSIFICATION (null when off); LoopOverTree only
        std::unique_ptr<Massification> p_massif;
        //! the provider's 1/eps coefficient (ME_E1 units) of the last call
        double E1() const;
        double Calc_V(const ATOOLS::Vec4D_Vector& p,
                    const double B,
                    const double mur);
        /*!
          YFS: RV_MODE 1. The provider's IR-finite one-loop over its own tree,
          lt = (alpha/2pi) V_fin/T, at the (n+1)-body point with renormalisation
          scale mur2. False when the point fails the NLO cuts, the tree
          vanishes, or the provider's finite part is not tree-normalised.
          p: the point in the lab (NLO cuts), peval: the same point in the
          frame the loop is evaluated in.
        */
        bool LoopOverTree(const ATOOLS::Vec4D_Vector &p,
                          const ATOOLS::Vec4D_Vector &peval, const double mur2,
                          double &lt);
        MODEL::Coupling_Map m_cpls;
        //! Not owned; set by SetProc, may be null.
        PHASIC::Process_Base *p_rvproc = nullptr;
        ATOOLS::Flavour_Vector m_flavs;
        double m_factor;
        //! RV_Generator as configured (p_rvproc is used only for cuts, not the value)
        std::string m_gen;
        inline std::string ActiveGenName() const { return m_gen.empty()?"none":m_gen; }
        inline void SetProc(PHASIC::Process_Base *realproc) {p_rvproc=realproc;}
        inline bool FailCut() {return m_failcut;}
        //! LoopOverTree calls, and those whose loop was exactly 0
        long m_ltcalls = 0, m_ltzero = 0;
    };
}



#endif
