#ifndef LL_Branching_H
#define LL_Branching_H

#include "NLL_Branching_Probability_Base.H"
#include "Running_AlphaS.H"

namespace APACIC { 

  class Splitting_Function;
  class Splitting_Group;

}

namespace PDF {

  class LL_Branching: public SHERPA::NLL_Branching_Probability_Base {
  public:

    class Initializer_Function {
    public:
      
      // destructor
      ~Initializer_Function();
      
      // member functions
      void Initialize() const;
      
    };// end of class Initializer_Function
    
    typedef std::set<APACIC::Splitting_Function*> SF_Set;

    const static double s_nfmax, s_nc, s_cf, s_ca;

    static SF_Set s_splittings;
    
  private:

    static Initializer_Function s_initializer;

    static void Insert(APACIC::Splitting_Function *splitting);

    static SF_Set::iterator Find(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b);

    friend class Initializer_Function;

  protected:

    ATOOLS::Flavour        m_flavour;
    MODEL::Running_AlphaS *p_alphas;

    std::vector<APACIC::Splitting_Group*> m_splittings;
    std::vector<double>                   m_const;

    void GenerateName();

  public:

    // constructor
    LL_Branching(const ATOOLS::Flavour flavour,MODEL::Running_AlphaS *alphas);

    // destructor
    ~LL_Branching();

    // member functions
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);

    // inline functions
    inline static const SF_Set &AllSplittings() { return s_splittings; }
    
    inline static double NfMax() { return s_nfmax; }
    inline static double NC()    { return s_nc;    }
    inline static double CF()    { return s_cf;    }
    inline static double CA()    { return s_ca;    }

  };// end of class LL_Branching_Probabilities

}// end of namespace PDF

#endif
