#ifndef XS_Base_H
#define XS_Base_H

namespace EXTRAXS { class XS_Base; }

#include "Integrable_Base.H"

namespace PHASIC { class Regulator_Base; }

namespace EXTRAXS {

  class XS_Base: public PHASIC::Integrable_Base {
  protected:

    int **p_colours;

    bool m_channels, m_psmc;

    double m_s, m_t, m_u;

    std::string m_gpath;

    static bool s_sortflavours;

    static std::string GenerateName(const size_t nin,const size_t nout,
				    const ATOOLS::Flavour *flavours);

    static void SortFlavours(ATOOLS::Flavour *flavs,const size_t &n);

    void Init(const ATOOLS::Flavour *flavours);

    friend class XS_Selector;
    friend class XS_Group;

    void PrepareTerminate();

   public:

    // constructors
    XS_Base(); 
    XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
	    const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
	    BEAM::Beam_Spectra_Handler *const beamhandler,
	    PDF::ISR_Handler *const isrhandler,
	    ATOOLS::Selector_Data *const selectordata);

    // destructor
    ~XS_Base();

    // member functions
    void Initialize(const PHASIC::scl::scheme scalescheme,
		    const int kfactorscheme,
		    BEAM::Beam_Spectra_Handler *const beamhandler,
		    PDF::ISR_Handler *const isrhandler,
		    ATOOLS::Selector_Data *const selectordata);

    void AssignRegulator(const std::string &regulator,
			 const std::vector<double> &parameters);

    virtual size_t Size() const;

    virtual XS_Base *const operator[](const size_t i) const;

    void SetSTU(const ATOOLS::Vec4D *momenta);

    double Differential(const ATOOLS::Vec4D *momenta);
    void   SetMax();

    void SwapInOrder();
    void RestoreInOrder();

    bool SelectOne();
    bool SelectOneFromList();
    void DeSelect();
    bool ReSelect(int);

    void FillSIntegrator(PHASIC::Multi_Channel *&mc);
    void UpdateIntegrator(PHASIC::Multi_Channel *&mc);

    virtual void Print();

    ATOOLS::Blob_Data_Base *SameWeightedEvent();

    virtual void Reset();
    virtual void ResetSelector(ATOOLS::Selector_Data *const selectordata);
    virtual void CreateISRChannels();
    
    virtual bool SetColours(const ATOOLS::Vec4D *momenta);
    virtual void SetMax(const double max,const int flag);

    virtual void SetISR(PDF::ISR_Handler *const isrhandler) = 0;
    virtual void SetTotal(const int mode=1) = 0;

    virtual double Differential(const double s,const double t,const double u) = 0;
    virtual bool   SetColours(const double s,const double t,const double u) = 0;
    virtual double operator()(const double s,const double t,const double u) = 0;

    virtual bool CalculateTotalXSec(const std::string &resultpath,
				    const bool create) = 0;
    virtual void WriteOutXSecs(std::ofstream &outfile) = 0;

    virtual bool SetPDFFlavours() = 0;
    virtual bool Tests() = 0;
    
    // inline functions
    inline int **const Colours() const { return p_colours; }

    inline const std::vector<ATOOLS::Flavour> &Resonances() const 
    { return m_resonances; }

    inline PHASIC::Regulator_Base *const Regulator() const 
    { return p_regulator; }

    inline void SetPSMC(const bool psmc) { m_psmc=psmc; }

    inline void SetGPath(const std::string &gpath) { m_gpath=gpath; }

    inline static void SetSortFlavours(const bool &sort)
    { s_sortflavours=sort; }

  };// end of class XS_Base

}// end of namespace EXTRAXS

#endif



