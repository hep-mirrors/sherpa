#ifndef LL_Sudakov_H
#define LL_Sudakov_H

#include "LL_Branching.H"
#include "LL_Single_Sudakov.H"
#include "Message.H"
#include <vector>
#include <map>

namespace SHERPA {

  class NLL_Sudakov_Base;

}

namespace SHERPA {

  class LL_Sudakov : public ATOOLS::Function_Base {
  public:

    typedef std::map<ATOOLS::Flavour,NLL_Sudakov_Base*> Sudakov_Map;

  private:

    std::vector<NLL_Sudakov_Base*> m_all_suds;
    Sudakov_Map m_sud_map;

    double m_lambda;
    double m_mu2,m_asmu;

    MODEL::Running_AlphaS    *p_alphas;

    double m_q2min, m_q2max;

    double Delta(const ATOOLS::Flavour &, double ,double );
    void   FixLambda2();    

    std::string m_outpath;

  public:

    LL_Sudakov(MODEL::Running_AlphaS *_p_alpha);
    ~LL_Sudakov();

    void AssignKeys(ATOOLS::Integration_Info *const info);

    double Lambda2() { return m_lambda; }

    inline void SetQ2Min(const double _m_q2min) { m_q2min=_m_q2min; }
    inline void SetQ2Max(const double _m_q2max) { m_q2max=_m_q2max; }
    
    void SetOutPath(const std::string _m_outpath);

    void Initialize();

    NLL_Sudakov_Base &Delta(const ATOOLS::Flavour &);

  };

}

#endif
