#ifndef LL_Branching_H
#define LL_Branching_H

#include "NLL_Branching_Probability_Base.H"
#include "Running_AlphaS.H"
#include "Splitting_Group.H"

namespace PDF {

  class LL_Branching: public SHERPA::NLL_Branching_Probability_Base {
  public:

    typedef std::set<APACIC::Splitting_Function*> SF_Set;

    const static unsigned int s_nf;
    static SF_Set             s_splitting;
    
  protected:

    ATOOLS::Flavour        m_flavour;
    MODEL::Running_AlphaS *p_alphas;

    std::vector<APACIC::Splitting_Group*> m_splitting;
    std::vector<double>                   m_const;

    void GenerateName();
    void Insert(APACIC::Splitting_Function *splitting);

    SF_Set::iterator Find(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b);

  public:

    // constructor
    LL_Branching();
    LL_Branching(const ATOOLS::Flavour flavour,MODEL::Running_AlphaS *alphas);

    // destructor
    ~LL_Branching();

    // getters
    const SF_Set &AllSplittings() const;

    // member functions
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);

    static void DeleteSplittings();

  };// end of class LL_Branching_Probabilities

  extern LL_Branching LLB;

}// end of namespace PDF

#endif
