#ifndef Integrable_Base_H
#define Integrable_Base_H

#include "Selector.H"
#include "Cut_Data.H"
#include "Flavour.H"
#include "Selector.H"
#include "Vector.H"
#include "Exception.H"

#include <string>

namespace ATOOLS { class Blob_Data_Base; }
namespace BEAM   { class Beam_Spectra_Handler; }
namespace PDF    { class ISR_Handler; }

namespace PHASIC {

  struct stp {

    enum id {
      ren     =  1,
      fac     =  2,
      kp21    =  3,
      kp22    =  4,
      as      = 10,
      unknown = 99
    };

  };// end of struct stp

  class Phase_Space_Handler;

  class Integrable_Base: public ATOOLS::Terminator_Object {
  protected:

    std::string m_name, m_resultpath, m_resultfile;
    size_t      m_nin, m_nout, m_nvector;

    ATOOLS::Flavour *p_flavours;
    ATOOLS::Vec4D   *p_momenta;

    int    m_scalescheme, m_kfactorscheme;
    double m_scalefactor, m_threshold, m_overflow;
    bool   m_updatescales;

    std::map<stp::id,double> m_scale;

    std::vector<double> m_xinfo;

    long unsigned int m_n;

    double m_last, m_lastlumi, m_lastdxs, m_max;
    double m_totalxs, m_totalsum, m_totalsumsqr, m_totalerr;
    double m_ssum, m_ssumsqr,m_smax,m_ssigma2;

    std::vector<double> m_vsmax;
    std::vector<long unsigned int> m_vsn;

    long unsigned int m_sn;

    bool m_swaped;

    Integrable_Base *p_selected;

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;

    PHASIC::Phase_Space_Handler *p_pshandler, *p_activepshandler; 
    ATOOLS::Selector_Base       *p_selector;
    ATOOLS::Cut_Data            *p_cuts;

    bool m_ownselector;

  protected:

    double TotalResult();
    double TotalVar();

  public:

    // constructor
    Integrable_Base(const size_t nin=0,const size_t nout=0,
		    const ATOOLS::Flavour *flavours=NULL,
		    const int scalescheme=0,const int kfactorscheme=0,
		    const double scalefactor=1.,
		    BEAM::Beam_Spectra_Handler *const beamhandler=NULL,
		    PDF::ISR_Handler *const isrhandler=NULL,
		    ATOOLS::Selector_Data *const selectordata=NULL);
  
    // destructor
    virtual ~Integrable_Base();
  
    // member functions
    void SetMomenta(const ATOOLS::Vec4D *momenta); 
    bool OneEvent(const double mass,const int mode);

    virtual bool Trigger(const ATOOLS::Vec4D *const momenta);

    virtual void SwapInOrder() = 0;
    virtual void RestoreInOrder() = 0;

    virtual double Differential(const ATOOLS::Vec4D *momenta) = 0;
    virtual double Differential2() = 0;
    virtual void   AddPoint(const double xs) = 0;

    virtual double CalculateScale(const ATOOLS::Vec4D *momenta) = 0;
    virtual double KFactor(const double scale) = 0;

    virtual void DeSelect() = 0;
    virtual void SelectOne() = 0;

    virtual void SetScale(const double scale);
    virtual void SetMax(const double max, int depth=0);
    virtual void SetMax();  
    virtual void ResetMax(int);  
    virtual bool OneEvent();
    virtual bool SameEvent();
    virtual void OptimizeResult();

    virtual ATOOLS::Blob_Data_Base *WeightedEvent(const int mode=0);
    virtual ATOOLS::Blob_Data_Base *SameWeightedEvent();

    virtual void SetPSHandler(Phase_Space_Handler *const pshandler);
    
    Integrable_Base *const Selected();

    // inline functions
    inline void SetName(const std::string &name) { m_name=name; }

    inline void SetScaleScheme(const int scalescheme)     
    { m_scalescheme=scalescheme;     }
    inline void SetKFactorScheme(const int kfactorscheme) 
    { m_kfactorscheme=kfactorscheme; }
    inline void SetScaleFactor(const double scalefactor)  
    { m_scalefactor=scalefactor;     }

    inline void SetTotalXS(const double totalxs)   { m_totalxs=totalxs;    }
    inline void SetTotalError(const double error)  { m_totalerr=error;     }
    inline void SetSum(const double sum)           { m_totalsum=sum;       }
    inline void SetSumSqr(const double sumsqr)     { m_totalsumsqr=sumsqr; }
    inline void SetPoints(const long int points)   { m_n=points;           }
    inline void SetSSum(const double sum)           { m_ssum=sum;       }
    inline void SetSSumSqr(const double sumsqr)     { m_ssumsqr=sumsqr; }
    inline void SetSigmaSum(const double sigma2)    { m_ssigma2=sigma2; }
    inline void SetSPoints(const long int points)   { m_sn=points;           }

    inline void SetOverflow(const double overflow) { m_overflow=overflow;  }

    inline void SetScale(const double scale,const stp::id type) 
    { m_scale[type]=scale; }

    inline void SetSelected(Integrable_Base *const selected) 
    { p_selected=selected; }
    
    inline size_t NIn() const     { return m_nin;     }
    inline size_t NOut() const    { return m_nout;    }
    inline size_t NVector() const { return m_nvector; }

    inline const std::string     &Name() const     { return m_name;     }
    inline const ATOOLS::Flavour *Flavours() const { return p_flavours; }
    inline const ATOOLS::Vec4D   *Momenta() const  { return p_momenta;  }

    inline double Last() const     { return m_last;     }
    inline double LastXS() const   { return m_lastdxs;  }
    inline double LastLumi() const { return m_lastlumi; }

    inline double   Sum() const    { return m_totalsum;    }
    inline double   SumSqr() const { return m_totalsumsqr; }
    inline long int Points() const { return m_n;           }

    inline double TotalXS() const    { return m_totalxs;  }
    inline double TotalError() const { return m_totalerr; }
    inline double Max() const        { return m_max;      }
    inline double Overflow() const   { return m_overflow; }

    inline double ISRThreshold() const      { return m_threshold;   }

    inline int    KFactorScheme() const { return m_kfactorscheme; }
    inline int    ScaleScheme() const   { return m_scalescheme;   }
    inline double ScaleFactor() const   { return m_scalefactor;   }

    inline const double Scale(const stp::id type=stp::as) 
    { return m_scale[type]; }

    inline int InSwaped() const { return m_swaped; }

    inline std::vector<double> &XInfo() { return m_xinfo; }

    inline ATOOLS::Cut_Data      *const Cuts() const     { return p_cuts;     }
    inline ATOOLS::Selector_Base *const Selector() const { return p_selector; }

    inline BEAM::Beam_Spectra_Handler *const Beam() const 
    { return p_beamhandler; }
    inline PDF::ISR_Handler           *const ISR() const  
    { return p_isrhandler;  }

    inline Phase_Space_Handler *const PSHandler(const bool own=true) 
    { if (own) return p_pshandler; return p_activepshandler; }
    
    inline void CreateMomenta(const size_t n)
    { delete p_momenta; p_momenta = new ATOOLS::Vec4D[n]; m_nvector=n; }

  };// end of class Integrable_Base

}// end of namespace PHASIC

#endif


