#ifndef XS_Base_H
#define XS_Base_H


#include "Phase_Space_Handler.H"
#include "Integrable_Base.H"

#include "ISR_Handler.H"
#include "Beam_Spectra_Handler.H"
#include "Selector.H"
#include "Combined_Selector.H"
#include "Standard_Selector.H"

#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"

#include <string>

namespace EXTRAXS { class XS_Base; }

namespace EXTRAXS {
  class XS_Base : public PHASIC::Integrable_Base {
  protected:
    int                    m_nin, m_nout;
    ATOOLS::Flavour      * p_fl;
    ATOOLS::Vec4D        * p_moms;
    std::string            m_name;

    double                 m_s, m_t, m_u;
    double                 m_scale,m_scalefactor,m_thres;

    int                    m_scalescheme,m_kfactorscheme;
    long int               m_n;
    double                 m_last,m_lastlumi,m_lastdxs;
    double                 m_totalxs,m_totalsum,m_totalsumsqr,m_totalerr, m_max;

    int m_swaped;

    int                 ** p_colours;
    XS_Base              * p_selected;

    PHASIC::Phase_Space_Handler    * p_ps; 
    BEAM::Beam_Spectra_Handler     * p_beam;
    PDF::ISR_Handler               * p_isr;
    ATOOLS::Selector_Base       * p_sel;

    void                             Init(ATOOLS::Flavour *);
    void                             GenerateName();
   public:
    XS_Base(int,int,ATOOLS::Flavour *,
	    PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,ATOOLS::Selector_Data *,
	    int,int,double);
    XS_Base(int,int,ATOOLS::Flavour *);
    XS_Base() {   m_name     = std::string("Empty XS"); }

    
    virtual ~XS_Base();

    virtual int    Size() { return 0; }
    virtual double Differential(ATOOLS::Vec4D *);
    virtual double Differential2()                              = 0;
    virtual double Scale(ATOOLS::Vec4D *);
    virtual bool   SetColours(ATOOLS::Vec4D *);
    virtual double KFactor(double);

    int                        Nin();
    int                        Nout();
    int                        Nvec();
    std::string                Name();
    ATOOLS::Flavour       * Flavs();
    double                     ISRThreshold();
    double                     Max();
    double                     Total();
    double                     Last();
    double                     LastXS();
    double                     LastLumi();
    double                     Scale();
    ATOOLS::Vec4D          * Momenta(); 
    int                     ** Colours();
    XS_Base                  * Selected();
    ATOOLS::Cut_Data      * Cuts();
    ATOOLS::Selector_Base * Selector();
    void                       UpdateCuts(double,double)  {};
    void                       SetMax(double);
    void                       SetTotal(double);
    void                       SetName(std::string);
    virtual void               SetISR(PDF::ISR_Handler *)          = 0;
    virtual void               SetTotalXS()                        = 0;
    virtual double             Differential(double,double,double)  = 0;
    virtual double             operator()(double,double,double)    = 0;
    virtual bool               SetColours(double, double, double)  = 0;
    virtual XS_Base          * operator[](int idx)                 = 0;
    virtual bool               CalculateTotalXSec()                =0;
    virtual void               AddPoint(double)                    = 0;
    virtual void               DeSelect()                          = 0;
    virtual void               SelectOne()                         = 0;
    virtual void               SetMomenta(ATOOLS::Vec4D *) {} //      = 0;
    virtual bool               OneEvent()                          = 0;
    //virtual bool               SameEvent()                         = 0;
    //virtual double             WeightedEvent()                     = 0;

    virtual void SwapInOrder();
    virtual void RestoreInOrder();
    int InSwaped();
  };

  inline int                        XS_Base::Nin()                   { return m_nin; }
  inline int                        XS_Base::Nout()                  { return m_nout; }  
  inline int                        XS_Base::Nvec()                  { return m_nin+m_nout; }  
  inline std::string                XS_Base::Name()                  { return m_name; }
  inline ATOOLS::Flavour       * XS_Base::Flavs()                 { return p_fl; }
  inline double                     XS_Base::ISRThreshold()          { return m_thres; }
  inline double                     XS_Base::Max()                   { return m_max; }
  inline double                     XS_Base::Total()                 { return m_totalxs; }
  inline double                     XS_Base::Last()                  { return m_last; }
  inline double                     XS_Base::LastXS()                { return m_lastdxs; }
  inline double                     XS_Base::LastLumi()              { return m_lastlumi; }
  inline double                     XS_Base::Scale()                 { return m_scale; }   
  inline ATOOLS::Vec4D          * XS_Base::Momenta()               { return p_moms; } 
  inline int                     ** XS_Base::Colours()               { return p_colours; }
  inline XS_Base                  * XS_Base::Selected()              { return p_selected; }
  inline ATOOLS::Cut_Data      * XS_Base::Cuts()                  { return NULL; }
  inline ATOOLS::Selector_Base * XS_Base::Selector()              { return p_sel; }

  inline void                       XS_Base::SetMax(double _max)             { m_max     = _max;   }
  inline void                       XS_Base::SetTotal(double _total)         { m_totalxs = _total; }
  inline void                       XS_Base::SetName(std::string _name)      { m_name    = _name; }
  inline int XS_Base::InSwaped() { return m_swaped; }
}

#endif

/*! 
  \class XS_Base 
  \brief The mother class of all processes : Single_XS, XS_Group, All_XSes.
  
  Every class, which is to be integrated by the Phase_Space_Integrator is derived 
  from this class or from Process_Base in AMEGIC. 
*/

/*! 
  \var   double   totalxs, totalerr, totalsum, totalsumsqr, max, long int n
  \brief Cross sections and their errors. 

  There's a double strategy : 
  First of all the relative contribution of individual channels is optimized within 
  the multichannel such that the overall variance is minimized. To do this, individual 
  results, variance for each /sa Single_Channel are to be stored within each 
  Single_Channel. 
  The total result and its error is accumulated in the /sa Phase_Space_Integrator and 
  in the full /sa Multi_Channel as well. It should be noted that this /sa Multi_Channel
  can also be a /sa Single_Channel only, both are derived from the /sa Channel_Base. 
  It should be stressed here that, if ISR is on, there are two Multi_Channels,
  one for th final state momenta and one for the initial state piece of the
  integration.
  
  On the other hand, results have to be distributed over the processes that are
  integrated. This is especially true for /sa XS_Groups. There the total result
  coincides with the one of the integration stored in /sa Phase_Space_Integrator
  and its specific /sa Channel_Base. However, this total result has to be distributed
  over the /sa Single_XS(es), therefore they have to accumulate their own
  result, variance, number of points etc..
*/

/*!
  \var double last, lastdxs, lastlumi;
  The results for the partial cross section at the last phase space point.
  This is important for instance for the integration, if one /sa Single_XS
  can be mapped onto another one. Then the last amplitude squared lastdxs of this 
  partner will be used and folded with the PDF (lastlumi) that depends on the 
  particular initial state. We have   last = lastlumi * lastdxs.
*/


