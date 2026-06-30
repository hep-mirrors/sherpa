#ifndef EXTRAXS_XS_QUARKONIA_H
#define EXTRAXS_XS_QUARKONIA_H

#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"


class QuarkoniumState {
    public:
        size_t index;
        kf_code kfc;
        double mass;
        double mass2;
        double ldme;
        double R02; 
        double R12;
        
        QuarkoniumState(const ATOOLS::Flavour_Vector& fl);
};

namespace Quarkonium{
    const ATOOLS::Flavour* FindState(const ATOOLS::Flavour_Vector& fl);

    namespace ColourFlowSetter{
        void qg_qSinglet(size_t m_q, size_t m_a, size_t m_g, size_t m_S, std::vector<std::vector<int>>& m_colours);
        void qqbar_gSinglet(size_t m_S, size_t m_a, std::vector<std::vector<int>>& m_colours);
        void gg_gSinglet(size_t m_S, std::vector<std::vector<int>>& m_colours);
        void qg_qOctet(size_t m_a, size_t m_g, size_t m_S, std::vector<std::vector<int>>& m_colours);
        void qqbar_gOctet(size_t m_S, size_t m_a, std::vector<std::vector<int>>& m_colours);
        void gg_gOctet(size_t m_S, std::vector<std::vector<int>>& m_colours);

    }
    namespace IsCandidate{
        bool qg_qQuarkonium(const ATOOLS::Flavour_Vector& fl, std::initializer_list<kf_code> allowed_fls);
        bool qqbar_gQuarkonium(const ATOOLS::Flavour_Vector& fl, std::initializer_list<kf_code> allowed_fls);
        bool gg_gQuarkonium(const ATOOLS::Flavour_Vector& fl, std::initializer_list<kf_code> allowed_fls);
    }
    namespace Calc{
        double qg_q1S0(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g1S0(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g1S0(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3S1(const double & s, const double & t, const double & u, const double & mass2);
        double qg_q3P0(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g3P0(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3P0(const double & s, const double & t, const double & u, const double & mass2);
        double qg_q3P1(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g3P1(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3P1(const double & s, const double & t, const double & u, const double & mass2);
        double qg_q3P2(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g3P2(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3P2(const double & s, const double & t, const double & u, const double & mass2);

        double qg_q1S0_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g1S0_oct(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g1S0_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qg_q3S1_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g3S1_oct(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3S1_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qg_q3P0_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g3P0_oct(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3P0_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qg_q3P1_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g3P1_oct(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3P1_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qg_q3P2_oct(const double & s, const double & t, const double & u, const double & mass2);
        double qqbar_g3P2_oct(const double & s, const double & t, const double & u, const double & mass2);
        double gg_g3P2_oct(const double & s, const double & t, const double & u, const double & mass2);
    }
}

namespace EXTRAXS {
  class XS_qg_q1S0 : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q1S0(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3S1 : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3S1(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3P0 : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3P0(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3P1 : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3P1(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3P2 : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3P2(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g1S0 : public ME2_Base {
    private:
        size_t m_a;
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g1S0(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g3S1 : public ME2_Base {
    private:
        size_t m_a;
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3S1(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

   class XS_qqbar_g3P0 : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3P0(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P1 : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3P1(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P2 : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3P2(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_gg_g1S0 : public ME2_Base {
    private:
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g1S0(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_gg_g3S1 : public ME2_Base {
    private:
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3S1(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };

  class XS_gg_g3P0 : public ME2_Base {
    private:
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3P0(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };

  class XS_gg_g3P1 : public ME2_Base {
    private:
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3P1(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };

  class XS_gg_g3P2 : public ME2_Base {
    private:
        double m_alphaS,  m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3P2(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };

  // OCTETS

  class XS_qg_q1S0_oct : public ME2_Base {
    private:
        size_t m_g, m_q, m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q1S0_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3S1_oct : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3S1_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3P0_oct : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3P0_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3P1_oct : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3P1_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qg_q3P2_oct : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qg_q3P2_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g1S0_oct : public ME2_Base {
    private:
        size_t  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g1S0_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g3S1_oct : public ME2_Base {
    private:
        size_t  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3S1_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

   class XS_qqbar_g3P0_oct : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3P0_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P1_oct : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3P1_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_qqbar_g3P2_oct : public ME2_Base {
    private:
        size_t m_g, m_q,  m_a;
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_qqbar_g3P2_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_gg_g1S0_oct : public ME2_Base {
    private:
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g1S0_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);
  };

  class XS_gg_g3S1_oct : public ME2_Base {
    private:
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3S1_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };

  class XS_gg_g3P0_oct : public ME2_Base {
    private:
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3P0_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };

  class XS_gg_g3P1_oct : public ME2_Base {
    private:
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3P1_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };

  class XS_gg_g3P2_oct : public ME2_Base {
    private:
        double m_alphaS, m_pref;
        QuarkoniumState Q;
    public:
        XS_gg_g3P2_oct(const PHASIC::External_ME_Args& args);
        double operator()(const ATOOLS::Vec4D_Vector& mom);
        bool SetColours(const ATOOLS::Vec4D_Vector& mom);

  };
}

#endif