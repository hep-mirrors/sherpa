#ifndef PHASIC_Main_Process_Integrator_H
#define PHASIC_Main_Process_Integrator_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Info_Key.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Terminator_Objects.H"

#include <string>

namespace ATOOLS {
  class Blob_Data_Base;
  class Cluster_Amplitude;
}
namespace BEAM   { class Beam_Spectra_Handler; }
namespace PDF    { class ISR_Handler; }
namespace YFS    { class YFS_Handler; }

namespace PHASIC {

  class Multi_Channel;
  class Color_Integrator;
  class Helicity_Integrator;

  class Phase_Space_Handler;

  class Process_Integrator: public ATOOLS::Terminator_Object {
  protected:

    Process_Base *p_proc;

    /// handler created by this instance, automatically deleted with this
    std::unique_ptr<Phase_Space_Handler> p_ownpshandler;
    /// pointer to actually used phase space handler,
    /// might be owned by different instance and not deleted by this
    Phase_Space_Handler * p_pshandler;

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;
    YFS::YFS_Handler           *p_yfshandler;

    std::string m_resultpath;
    size_t      m_nin, m_nout, m_smode, m_swmode;

    ATOOLS::Vec4D_Vector p_momenta;

    cls::scheme m_colorscheme;
    hls::scheme m_helicityscheme;

    double m_threshold, m_enhancefac, m_maxeps, m_rsfac;

    long unsigned int m_n, m_itmin, m_itmax;

    double m_max, m_totalxs, m_totalsum, m_totalsumsqr, m_totalerr;
    double m_ssum, m_ssumsqr, m_smax, m_ssigma2, m_wmin;
    double m_mssum, m_mssumsqr, m_msn;

    std::vector<double> m_vsmax, m_vsum;
    std::vector<long unsigned int> m_vsn;

    long unsigned int m_sn, m_son;

    bool m_writeout;

    ATOOLS::Histogram *p_whisto;

    std::shared_ptr<Color_Integrator>    p_colint;
    std::shared_ptr<Helicity_Integrator> p_helint;

    double Sigma2() const;

  public:

    double TotalSigma2() const;
    double TotalResult() const;
    double TotalVar() const;

  public:

    // constructor
    Process_Integrator(Process_Base *const proc);

    // destructor
    ~Process_Integrator();
  
    // member functions
    bool Initialize(BEAM::Beam_Spectra_Handler *const beamhandler=NULL,
		    PDF::ISR_Handler *const isrhandler=NULL,
		    YFS::YFS_Handler *const yfshandler=NULL);
  
    void SetMomenta(const ATOOLS::Vec4D_Vector &p);
    void SetMomenta(const ATOOLS::Cluster_Amplitude &ampl);

    void AddPoint(const double xs);
    void Reset(const int mode=1);

    void SetMax(const double max);
    void SetTotal(const int mode=1);

    void ResetMax(int);  
    void OptimizeSubResult(const double &s2);

    //void SetPSHandler(const std::shared_ptr<Phase_Space_Handler> &pshandler);
    void SetPSHandler(Phase_Space_Handler * pshandler);
    void SetPSHandler(const double &maxerr,const std::string eobs,const std::string efunc);
    void SetISRThreshold(const double threshold);
    
    void InitWeightHistogram();
    void ReadInHistogram(std::string);
    void WriteOutHistogram(std::string);
    bool ReadInXSecs(const std::string &path);
    void WriteOutXSecs(const std::string &path);

    void OptimizeResult();
    void EndOptimize();

    double GetMaxEps(double);
    void   SetUpEnhance(const int omode=0);

    void ReadResults();
    void StoreResults(const int mode=0);
    void StoreBackupResults();
    void PrepareTerminate();

    void SetEnhanceFactor(const double &efac);
    void SetMaxEpsilon(const double &maxeps);
    void SetRSEnhanceFactor(const double &rsfac);

    double SelectionWeight(const int mode) const;

    void MPISync(const int mode=0);
    void MPICollect(std::vector<double> &sv,std::vector<double> &mv,size_t &i);
    void MPIReturn(std::vector<double> &sv,std::vector<double> &mv,size_t &i);

    // inline functions
    inline Process_Base *Process() const { return p_proc; }

    inline void SetColorScheme(const cls::scheme &s)    { m_colorscheme=s;    }
    inline void SetHelicityScheme(const hls::scheme &s) { m_helicityscheme=s; }

    inline void SetItMin(const long unsigned int &itmin) { m_itmin=itmin;  }
    inline void SetItMax(const long unsigned int &itmax) { m_itmax=itmax;  }

    inline size_t NIn() const     { return m_nin;     }
    inline size_t NOut() const    { return m_nout;    }

    inline const ATOOLS::Vec4D_Vector &Momenta() const { return p_momenta;  }
    inline ATOOLS::Vec4D_Vector       &Momenta()       { return p_momenta;  }

    inline double   Sum() const     { return m_totalsum;    }
    inline double   SumSqr() const  { return m_totalsumsqr; }
    inline long int Points() const  { return m_n+m_sn;      }
    inline long int SPoints() const { return m_sn;          }

    inline double TotalXS() const    { return m_totalxs;  }
    inline double TotalError() const { return m_totalerr; }
    inline double Max() const        { return m_max;      }

    inline double ISRThreshold() const      { return m_threshold;   }
    inline double EnhanceFactor() const     { return m_enhancefac;  }
    inline double RSEnhanceFactor() const   { return m_rsfac;       }

    inline cls::scheme ColorScheme() const    { return m_colorscheme;    }
    inline hls::scheme HelicityScheme() const { return m_helicityscheme; }

    inline void SetBeam(BEAM::Beam_Spectra_Handler *const beam) 
    { p_beamhandler=beam; }
    inline void SetISR(PDF::ISR_Handler *const isr)  
    { p_isrhandler=isr;  }
    inline void SetYFS(YFS::YFS_Handler *const yfs)
    { p_yfshandler=yfs;  }

    inline BEAM::Beam_Spectra_Handler *Beam() const 
    { return p_beamhandler; }
    inline PDF::ISR_Handler           *ISR() const  
    { return p_isrhandler;  }

    inline YFS::YFS_Handler           *YFS() const
    { return p_yfshandler;  }
    //inline std::shared_ptr<Phase_Space_Handler> PSHandler() const
    //{ return p_pshandler; }

    inline Phase_Space_Handler * PSHandler() const { return p_pshandler; }

    inline void SetColorIntegrator(const std::shared_ptr<Color_Integrator> colint)
    { p_colint=colint; }
    inline void SetHelicityIntegrator(const std::shared_ptr<Helicity_Integrator> helint)
    { p_helint=helint; }

    inline std::shared_ptr<Color_Integrator> ColorIntegrator() const
    { return p_colint; }
    inline std::shared_ptr<Helicity_Integrator> HelicityIntegrator() const
    { return p_helint; }

    inline long unsigned int ItMin() const { return m_itmin; }
    inline long unsigned int ItMax() const { return m_itmax; }

    inline void SetResultPath(const std::string &path) { m_resultpath=path; }

    inline std::string ResultPath() const { return m_resultpath; }

    inline ATOOLS::Histogram *WeightHisto() const { return p_whisto; }

  };// end of class Process_Integrator

}// end of namespace PHASIC

#endif


