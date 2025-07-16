#ifndef PHASIC_Process_Color_Correlated_ME2_H
#define PHASIC_Process_Color_Correlated_ME2_H

#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Math/Vector.H"
#include "MODEL/Main/Model_Base.H"

namespace MODEL { struct Coupling_Map; }

namespace PHASIC {

  /* Color-correlated matrix elements [arXiv:hep-ph/9605323]
     of the form <1,2,3,... | TiTj |...,3,2,1 > */

  class Color_Correlated_ME2 {

  protected:
    
    typedef std::vector< std::vector<double> > CCMatrix;

    MODEL::Coupling_Data* p_aqcd, * p_aqed;

  public:

    Color_Correlated_ME2 (const External_ME_Args& args);

    virtual ~Color_Correlated_ME2() {}

    /* Calculate and store full color-correlation matrix */
    virtual void Calc(const ATOOLS::Vec4D_Vector &p) = 0;
    
    /* Get element of color-correlation matrix:
       GetValue(i,j) = <1,2,3,... | TiTj | ...,3,2,1 > */
    virtual double GetValue(const size_t& i, const size_t& j) const = 0;

    /* Get squared uncorrelated amplitude (comes for free in typical
       implementations) */
    virtual double GetBorn2() const = 0;

    virtual void SetCouplings(const MODEL::Coupling_Map& cpls);
    
    virtual double AlphaQCD() const;
    virtual double AlphaQED() const;

    static Color_Correlated_ME2 *GetME2(const PHASIC::External_ME_Args& args);
    static Color_Correlated_ME2 *GetME2(const std::string& tag,
					const PHASIC::External_ME_Args& args);

  };

}


#define DECLARE_COLORCORRELATEDME2_GETTER(NAME,TAG)		                \
  DECLARE_GETTER(NAME,TAG,PHASIC::Color_Correlated_ME2,PHASIC::External_ME_Args);   \
  void ATOOLS::Getter<PHASIC::Color_Correlated_ME2,PHASIC::External_ME_Args,NAME>:: \
  PrintInfo(std::ostream &str,const size_t width) const		                \
  {                                                                             \
    str<<#TAG;                                                                  \
  }

#endif
