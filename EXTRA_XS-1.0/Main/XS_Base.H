#ifndef XS_Base_H
#define XS_Base_H

namespace EXTRAXS { class XS_Base; }

#include "Integrable_Base.H"
#include "Regulator_Base.H"

namespace EXTRAXS {

  class XS_Base: public PHASIC::Integrable_Base {
  protected:

    Regulator_Base *p_regulator;

    int **p_colours;

    double m_s, m_t, m_u;

    std::vector<ATOOLS::Flavour> m_resonances;

    static std::string GenerateName(const size_t nin,const size_t nout,
				    const ATOOLS::Flavour *flavours);

    void Init(const ATOOLS::Flavour *flavours);

   public:

    // constructors
    XS_Base(); 
    XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    XS_Base(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
	    const int scalescheme,const int kfactorscheme,const double scalefactor,
	    BEAM::Beam_Spectra_Handler *const beamhandler,PDF::ISR_Handler *const isrhandler,
	    ATOOLS::Selector_Data *const selectordata);

    // destructor
    ~XS_Base();

    // member functions
    void AssignRegulator(const std::string &regulator,
			 const std::vector<double> &parameters);

    virtual size_t Size() const;

    virtual XS_Base *const operator[](const size_t i) const;

    void SetSTU(const ATOOLS::Vec4D *momenta);

    double Differential(const ATOOLS::Vec4D *momenta);
    void   SetMax();

    double CalculateScale(const ATOOLS::Vec4D *momenta);
    double KFactor(const double scale);

    void SwapInOrder();
    void RestoreInOrder();

    void SelectOne();
    void DeSelect();

    virtual void Print();

    ATOOLS::Blob_Data_Base *SameWeightedEvent();

    virtual void Reset();
    virtual void ResetSelector(ATOOLS::Selector_Data *const selectordata);
    
    virtual bool SetColours(const ATOOLS::Vec4D *momenta);
    virtual void SetMax(const double max,const int flag);

    virtual void SetISR(PDF::ISR_Handler *const isrhandler) = 0;
    virtual void SetTotal() = 0;

    virtual double Differential(const double s,const double t,const double u) = 0;
    virtual bool   SetColours(const double s,const double t,const double u) = 0;
    virtual double operator()(const double s,const double t,const double u) = 0;

    virtual bool CalculateTotalXSec(const std::string &resultpath,
				    const bool create) = 0;
    virtual void WriteOutXSecs(std::ofstream &outfile) = 0;

    // inline functions
    inline int **const Colours() const { return p_colours; }

    inline const std::vector<ATOOLS::Flavour> &Resonances() const { return m_resonances; }

  };// end of class XS_Base

}// end of namespace EXTRAXS

#endif



