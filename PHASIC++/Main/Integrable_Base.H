#ifndef Integrable_Base_H
#define Integrable_Base_H

#include "Selector.H"
#include "Jet_Finder_Base.H"
#include "Cut_Data.H"
#include "Flavour.H"
#include "Vector.H"
#include "Exception.H"
#include "Histogram.H"
#include "Spin_Correlation_Tensor.H"
#include "Info_Key.H"
#include "Algebra_Interpreter.H"

#include <string>

namespace ATOOLS { class Blob_Data_Base; }
namespace HELICITIES { class Amplitude_Tensor; }
namespace BEAM   { class Beam_Spectra_Handler; }
namespace PDF    { class ISR_Handler; }

namespace PHASIC {

  class Multi_Channel;

  struct stp {

    enum id {
      unknown =  0,
      ren     =  1,
      fac     =  2,
      kp21    =  3,
      kp22    =  4,
      sfs     =  8,
      sis     = 16
    };

  };// end of struct stp

  struct scl {

    enum scheme {
      unknown   = 0,
      ckkw      = 1,
      fixed     = 2,
      shat      = 4,
      ameanpt   = 8,
      gmeanpt   = 16,
      minpt     = 32,
      sumpt     = 64,
      strong    = 128,
      mass      = 256,
      div_by_2  = 512,
      mult_by_2 = 1024, 
      updf      = 2048,
      bfkl      = 4096,
      reggeise  = 8192
    };

  };// end of struct scl

  struct cls {

    enum scheme {
      unknown   = 0,
      sum       = 1,
      sample    = 2
    };

  };// end of struct cls

  std::ostream &operator<<(std::ostream &str,const cls::scheme &s);

  struct hls {

    enum scheme {
      unknown   = 0,
      sum       = 1,
      sample    = 2
    };

  };// end of struct hls

  std::ostream &operator<<(std::ostream &str,const hls::scheme &s);

  inline const scl::scheme operator|(const scl::scheme c1,const scl::scheme c2)
  { return (scl::scheme)((int)c1|(int)c2); }
  inline const scl::scheme operator&(const scl::scheme c1,const scl::scheme c2)
  { return (scl::scheme)((int)c1&(int)c2); }
  inline const scl::scheme operator+(const scl::scheme c1,const scl::scheme c2)
  { return (scl::scheme)((int)c1|(int)c2); }
  inline const scl::scheme operator-(const scl::scheme c1,const scl::scheme c2)
  { return (scl::scheme)((int)c1&(int)~c2);}

  inline scl::scheme & operator+=(scl::scheme &c1,const scl::scheme c2)
  { return c1=c1+c2; }
  inline scl::scheme & operator-=(scl::scheme &c1,const scl::scheme c2)
  { return c1=c1-c2; }

  class Integrable_Base;

  class Tag_Setter: public ATOOLS::Tag_Replacer {
  private:

    Integrable_Base *p_ib;

    ATOOLS::Algebra_Interpreter *p_calc;

  public:
    
    // constructor
    inline Tag_Setter(Integrable_Base *const ib): p_ib(ib), p_calc(NULL) {}
    
    // member functions
    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    // inline functions
    void SetCalculator(ATOOLS::Algebra_Interpreter *const calc) { p_calc=calc; }

  };// end of class Tag_Setter

  class Regulator_Base;
  class Phase_Space_Handler;

  class Integrable_Base : public ATOOLS::Terminator_Object {
  protected:

    std::string m_name, m_resultpath, m_resultfile, m_histofile;
    size_t      m_nin, m_nout, m_nvector, m_corenout;

    ATOOLS::Flavour *p_flavours;
    ATOOLS::Vec4D   *p_momenta;

    scl::scheme m_scalescheme;
    cls::scheme m_colorscheme;
    hls::scheme m_helicityscheme;
    int         m_kfactorscheme, m_maxjetnumber, m_coremaxjetnumber;
    int         m_nstrong, m_neweak, m_orderQCD, m_orderEW, m_usepi;
    double      m_threshold, m_overflow, m_enhancefac;
    double      m_rfactor, m_ps_kfactor, m_ps_cpl_factor;

    std::map<stp::id,double>     m_scale;
    std::vector<ATOOLS::Flavour> m_resonances;
    std::vector<double>          m_xinfo;

    long unsigned int m_n;

    long int m_expevents, m_dicedevents, m_accevents;

    double m_anasum, m_validanasum, m_gmin;

    double m_last, m_lastlumi, m_lastdxs, m_max;
    double m_totalxs, m_totalsum, m_totalsumsqr, m_totalerr;
    double m_ssum, m_ssumsqr,m_smax,m_ssigma2,m_wmin;

    double m_procweight, m_me_as_factor, m_ycut, m_cycut;

    std::string m_cuttag;

    std::vector<double> m_vsmax;
    std::vector<long unsigned int> m_vsn;

    long unsigned int m_sn,m_son;

    bool m_swaped;

    Integrable_Base *p_selected, *p_parent;
    Regulator_Base  *p_regulator;

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;

    PHASIC::Phase_Space_Handler *p_pshandler, *p_activepshandler; 
    ATOOLS::Selector_Base       *p_selector;
    ATOOLS::Cut_Data            *p_cuts;
    ATOOLS::Histogram           *p_whisto;

    ATOOLS::Jet_Finder_Base *p_jf;

    bool m_ownselector;

    std::string m_efunc, m_muf2tag, m_mur2tag;

    ATOOLS::Info_Key m_kfkey;

    ATOOLS::Algebra_Interpreter *p_muf2calc, *p_mur2calc;

    Tag_Setter m_muf2tagset, m_mur2tagset;

  public:

    double TotalResult();
    double TotalVar();
    double RemainTimeFactor(double maxerr); 

  public:

    // constructor
    Integrable_Base(const size_t nin=0,const size_t nout=0,
		    const scl::scheme scalescheme=scl::unknown,
		    const int kfactorscheme=0,
		    BEAM::Beam_Spectra_Handler *const beamhandler=NULL,
		    PDF::ISR_Handler *const isrhandler=NULL,
		    ATOOLS::Selector_Data *const selectordata=NULL,
		    const cls::scheme &clsc=cls::sum,
		    const hls::scheme &hlsc=hls::sum);
  
    // destructor
    virtual ~Integrable_Base();
  
    // member functions
    void SetMomenta(const ATOOLS::Vec4D *momenta); 

    ATOOLS::Blob_Data_Base *OneEvent(const double mass,const int mode);

    virtual bool Trigger(const ATOOLS::Vec4D *const momenta);

    virtual void SwapInOrder() = 0;
    virtual void RestoreInOrder() = 0;

    virtual double Differential(const ATOOLS::Vec4D *momenta) = 0;
    virtual double Differential2() = 0;
    virtual void   AddPoint(const double xs);
    virtual ATOOLS::Spin_Correlation_Tensor* GetSpinCorrelations();
    virtual void   FillAmplitudes(HELICITIES::Amplitude_Tensor*,double=1.);

    virtual double CalculateScale(const ATOOLS::Vec4D *momenta);
    virtual double KFactor();

    virtual void DeSelect() = 0;
    virtual bool ReSelect(int) = 0;
    virtual bool SelectOne() = 0;

    virtual void GetGMin(double &gmin, double &meff);
    virtual void SetMax(const double max, int depth=0);
    virtual void SetMax();  
    virtual void ResetMax(int);  
    virtual void OptimizeResult();

    virtual ATOOLS::Blob_Data_Base *OneEvent();
    virtual ATOOLS::Blob_Data_Base *SameEvent();

    virtual ATOOLS::Blob_Data_Base *WeightedEvent(const int mode=0);
    virtual ATOOLS::Blob_Data_Base *SameWeightedEvent();

    virtual void SetPSHandler(Phase_Space_Handler *const pshandler);
    virtual void SetISRThreshold(const double threshold);
    
    virtual void InitWeightHistogram();
    virtual void ReadInHistogram(std::string);
    virtual void WriteOutHistogram(std::string);

    virtual double TriggerEfficiency();
    
    double GetMaxEps(double);
    void   SetMomenta();

    Integrable_Base *const Selected();
    Integrable_Base *const Parent();

    void CreateMomenta(const size_t n);

    static std::map<std::string,std::string> ScaleTags();
    static std::map<std::string,std::string> ColorSchemeTags();
    static std::map<std::string,std::string> HelicitySchemeTags();

    void SetScaleScheme(const scl::scheme s);

    virtual void SetFactorizationScale(const std::string &muf2);
    virtual void SetRenormalizationScale(const std::string &mur2);

    void SetFactorizationScale();
    void SetRenormalizationScale();

    virtual bool FillSIntegrator(Multi_Channel *&mc);
    virtual void UpdateIntegrator(Multi_Channel *&mc);

    // inline functions
    inline void SetName(const std::string &name) { m_name=name; }

    inline void SetColorScheme(const cls::scheme &s)    { m_colorscheme=s;    }
    inline void SetHelicityScheme(const hls::scheme &s) { m_helicityscheme=s; }

    inline void SetKFactorScheme(const int k) { m_kfactorscheme=k; }
    inline void SetUsePI(const int usepi)     { m_usepi=usepi;     }

    inline void SetTotalXS(const double totalxs)   { m_totalxs=totalxs;    }
    inline void SetTotalError(const double error)  { m_totalerr=error;     }
    inline void SetSum(const double sum)           { m_totalsum=sum;       }
    inline void SetSumSqr(const double sumsqr)     { m_totalsumsqr=sumsqr; }
    inline void SetPoints(const long int points)   { m_n=points;           }
    inline void SetSSum(const double sum)           { m_ssum=sum;       }
    inline void SetSSumSqr(const double sumsqr)     { m_ssumsqr=sumsqr; }
    inline void SetSigmaSum(const double sigma2)    { m_ssigma2=sigma2; }
    inline void SetSPoints(const long int points)   { m_sn=points;      }
    inline void SetWMin(const double wmin)          { m_wmin=wmin;      }
    inline void SetOptCounter(const long int son)   { m_son=son;        }

    inline void SetOverflow(const double overflow) { m_overflow=overflow;  }

    inline void SetScale(const double scale,const stp::id type) 
    { m_scale[type]=scale; }

    inline void SetSelected(Integrable_Base *const selected) 
    { p_selected=selected; }
    inline void SetParent(Integrable_Base *const parent) 
    { p_parent=parent; }
    
    inline void SetEnhanceFunction(const std::string &efunc)
    { m_efunc=efunc; }

    inline size_t NIn() const     { return m_nin;     }
    inline size_t NOut() const    { return m_nout;    }
    inline size_t NVector() const { return m_nvector; }

    inline const std::string     &Name() const     { return m_name;     }
    inline const ATOOLS::Flavour *Flavours() const { return p_flavours; }
    inline const ATOOLS::Vec4D   *Momenta() const  { return p_momenta;  }

    inline double Last() const     { return m_last;     }
    inline double LastXS() const   { return m_lastdxs;  }
    inline double LastLumi() const { return m_lastlumi; }

    inline double   Sum() const     { return m_totalsum;    }
    inline double   SumSqr() const  { return m_totalsumsqr; }
    inline long int Points() const  { return m_n;           }
    inline long int SPoints() const { return m_sn;          }

    inline long int ExpectedEvents() const { return m_expevents;   }
    inline long int DicedEvents() const    { return m_dicedevents; }
    inline double   ProcWeight() const     { return m_procweight;  }
    inline void     SetPW(const double w)  { m_procweight=w;       }

    inline double AnalysisSum() const      { return m_anasum;      }
    inline double ValidAnalysisSum() const { return m_validanasum; }
    
    inline double TotalXS() const    { return m_totalxs;  }
    inline double TotalError() const { return m_totalerr; }
    inline double Max() const        { return m_max;      }
    inline double Overflow() const   { return m_overflow; }

    inline double ISRThreshold() const      { return m_threshold;   }
    inline double EnhanceFactor() const     { return m_enhancefac;  }

    inline cls::scheme ColorScheme() const    { return m_colorscheme;    }
    inline hls::scheme HelicityScheme() const { return m_helicityscheme; }

    inline int         KFactorScheme() const { return m_kfactorscheme; }
    inline scl::scheme ScaleScheme()  const  { return m_scalescheme;   }

    inline std::string FactorizationScale() const { return m_muf2tag; }

    inline int NStrong() const { return m_nstrong; }
    inline int NEWeak() const  { return m_neweak;  }

    inline int OrderEWeak() const  { return m_orderEW;  }
    inline int OrderStrong() const { return m_orderQCD; }

    inline const double Scale(const stp::id type=stp::ren) 
    { return m_scale[type]; }

    inline int InSwaped() const { return m_swaped; }

    inline std::vector<double> &XInfo() { return m_xinfo; }

    inline ATOOLS::Cut_Data      *const Cuts() const     
    { return p_cuts;     }
    inline ATOOLS::Selector_Base *const Selector() const 
    { return p_selector; }
    virtual void FillOnshellConditions() {}

    inline BEAM::Beam_Spectra_Handler *const Beam() const 
    { return p_beamhandler; }
    inline PDF::ISR_Handler           *const ISR() const  
    { return p_isrhandler;  }

    inline Phase_Space_Handler *const PSHandler(const bool own=true) 
    { if (own) return p_pshandler; return p_activepshandler; }

    inline ATOOLS::Jet_Finder_Base *JetFinder() const { return p_jf; }
    
    inline std::string EnhanceFunction() const { return m_efunc; }

    inline double GMin() const { return m_gmin; }
    inline double Ycut() const { return m_ycut; }

    inline std::string YcutTag() const { return m_cuttag; }

  };// end of class Integrable_Base

}// end of namespace PHASIC

#endif


