#ifndef Integrable_Base_H
#define Integrable_Base_H

#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Selector.H"
#include "Cut_Data.H"
#include "Flavour.H"
#include "Selector.H"
#include "Vector.H"
#include "Exception.H"

#include <string>

namespace ATOOLS { class Blob_Data_Base; }

namespace PHASIC {

  struct stp {

    enum id {
      as      =  1,
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
    double m_scale, m_scalefactor, m_threshold, m_overflow;

    long int m_n;
    double   m_last, m_lastlumi, m_lastdxs, m_max;
    double   m_totalxs, m_totalsum, m_totalsumsqr, m_totalerr;

    bool m_swaped;

    Integrable_Base *p_selected;

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;

    PHASIC::Phase_Space_Handler *p_pshandler; 
    ATOOLS::Selector_Base       *p_selector;
    ATOOLS::Cut_Data            *p_cuts;

  public:

    // constructor
    Integrable_Base(const size_t nin=0,const size_t nout=0,const ATOOLS::Flavour *flavours=NULL,
		    const int scalescheme=0,const int kfactorscheme=0,const double scalefactor=0.,
		    BEAM::Beam_Spectra_Handler *const beamhandler=NULL,
		    PDF::ISR_Handler *const isrhandler=NULL,
		    ATOOLS::Selector_Data *const selectordata=NULL);
  
    // destructor
    virtual ~Integrable_Base();
  
    // member functions
    void SetMomenta(const ATOOLS::Vec4D *momenta); 

    virtual void SwapInOrder() = 0;
    virtual void RestoreInOrder() = 0;

    virtual double Differential(const ATOOLS::Vec4D *momenta) = 0;
    virtual double Differential2() = 0;
    virtual void   AddPoint(const double xs) = 0;

    virtual double Scale(const ATOOLS::Vec4D *momenta) = 0;
    virtual double KFactor(const double scale) = 0;

    virtual void DeSelect() = 0;
    virtual void SelectOne() = 0;

    virtual void SetMax();  
    virtual bool OneEvent();
    virtual bool SameEvent();

    virtual ATOOLS::Blob_Data_Base *WeightedEvent();
    virtual ATOOLS::Blob_Data_Base *SameWeightedEvent();

    Integrable_Base *const Selected();

    // inline functions
    inline void SetName(const std::string &name) { m_name=name; }

    inline void SetScaleScheme(const int scalescheme)     { m_scalescheme=scalescheme;     }
    inline void SetKFactorScheme(const int kfactorscheme) { m_kfactorscheme=kfactorscheme; }

    inline void SetTotalXS(const double totalxs)   { m_totalxs=totalxs;   }
    inline void SetTotalError(const double error)  { m_totalerr=error;    }
    inline void SetMax(const double max)           { m_max=max;           }
    inline void SetOverflow(const double overflow) { m_overflow=overflow; }

    inline void SetScale(const double scale) { m_scale=scale; }

    inline void SetSelected(Integrable_Base *const selected) { p_selected=selected; }

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

    inline double Scale(const stp::id type) const { return m_scale;     }
    inline double ISRThreshold() const            { return m_threshold; }

    inline int    KFactorScheme() const { return m_kfactorscheme; }
    inline int    ScaleScheme() const   { return m_scalescheme;   }
    inline double ScaleFactor() const   { return m_scalefactor;   }

    inline int InSwaped() const { return m_swaped; }

    inline ATOOLS::Cut_Data      *const Cuts() const     { return p_cuts;     }
    inline ATOOLS::Selector_Base *const Selector() const { return p_selector; }

    inline BEAM::Beam_Spectra_Handler *const Beam() const { return p_beamhandler; }
    inline PDF::ISR_Handler           *const ISR() const  { return p_isrhandler;  }

  };// end of class Integrable_Base

}// end of namespace PHASIC

#endif


