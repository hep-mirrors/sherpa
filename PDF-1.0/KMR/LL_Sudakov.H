#ifndef LL_Sudakov_H
#define LL_Sudakov_H

#include "LL_Branching.H"
#include "LL_Single_Sudakov.H"
#include "Message.H"
#include <vector>
#include <map>

namespace SHERPA { class NLL_Sudakov_Base; }

namespace PDF {

  class LL_Sudakov : public ATOOLS::Function_Base {
  public:

    typedef std::map<ATOOLS::Flavour,SHERPA::NLL_Sudakov_Base*> Sudakov_Map;

  private:

    std::vector<SHERPA::NLL_Sudakov_Base*> m_all_suds;
    Sudakov_Map m_sud_map;

    double m_lambda;
    double m_mu2,m_asmu;

    MODEL::Running_AlphaS    *p_alphas;

    double m_q2min, m_q2max;

    double Delta(const ATOOLS::Flavour &, double ,double );
    void   FixLambda2();    

    std::string m_outpath;

  public:

    // constructor
    LL_Sudakov(MODEL::Running_AlphaS *_p_alpha);

    // destructor
    ~LL_Sudakov();

    // member functions
    void AssignKeys(ATOOLS::Integration_Info *const info);

    void Initialize();

    SHERPA::NLL_Sudakov_Base &Delta(const ATOOLS::Flavour &);

    // inline functions
    inline double Lambda2() { return m_lambda; }

    inline double Q2Min() const { return m_q2min; }
    inline double Q2Max() const { return m_q2max; }
    
    inline const std::string &OutPath() const 
    { return m_outpath; }

    inline void SetQ2Min(const double q2min) { m_q2min=q2min; }
    inline void SetQ2Max(const double q2max) { m_q2max=q2max; }
    
    inline void SetOutPath(const std::string outpath) 
    { m_outpath=outpath; }

  };// end of class LL_Sudakov

}// end of namspace PDF

#endif
