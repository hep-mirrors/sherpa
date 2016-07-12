#ifndef MCATNLO_Showers_Sudakov_H
#define MCATNLO_Showers_Sudakov_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PDF/Main/PDF_Base.H"
#include "MCATNLO/Showers/Splitting_Function_Group.H"
#include "MCATNLO/Tools/Parton.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Phys/Variations.H"

namespace MCATNLO {
  typedef std::map<ATOOLS::Flavour,Splitting_Function_Base*> SF_E_Map;
  typedef std::map<ATOOLS::Flavour,SF_E_Map> SF_EE_Map;
  typedef std::map<ATOOLS::Flavour,SF_EE_Map> SF_EEE_Map;

  class Shower;

  class Sudakov : public Splitting_Function_Group {
    double m_k0sqi, m_k0sqf, m_phi;
    double m_deltaz, m_kperp2, m_z, m_y, m_x, m_ymin, m_ymax, m_scale;
    double m_zmin,m_zmax, m_weight;
    int m_ewmode, m_scalescheme;
    ATOOLS::Mass_Selector *  p_rms;
    ATOOLS::Flavour m_cfl, m_flspec;
    Parton * p_split, * p_spect;
    SF_EEE_Map m_fffmap, m_ffimap, m_iffmap, m_ifimap, m_fifmap, m_fiimap;
    SF_EEE_Map m_sffmap, m_sfimap, m_sifmap, m_siimap;
    std::vector<Splitting_Function_Base *> m_addsplittings;
    Shower *p_shower;
    std::vector<SFC_Getter*> m_cgets;
    ATOOLS::Variation_Weights *p_variationweights;
    bool m_reweightalphas, m_reweightpdfs;
    std::vector<size_t> m_disallowflav;
    void ProduceT();
    bool Veto(double,double);
    bool Splitting(double,double);
    bool DefineFFBoundaries(double,double);
    bool DefineFIBoundaries(double,double,int);
    bool DefineIFBoundaries(double,double,int);
    bool DefineIIBoundaries(double,double,int);
    void Add(Splitting_Function_Base * split);
    void AddToMaps(Splitting_Function_Base * split,const int mode=1);
    double Reweight(ATOOLS::Variation_Parameters *,
                    ATOOLS::Variation_Weights *,
                    const bool &success);
  public:
    Sudakov(PDF::ISR_Handler *isr,const int qed);
    ~Sudakov();

    void InitSplittingFunctions(MODEL::Model_Base *md,const int kfmode);
    void SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    void SetReweightAlphaS(bool v) { m_reweightalphas = v; }
    void SetReweightPDFs(bool v) { m_reweightpdfs = v; }
    bool   Generate(Parton *);
    double OverIntegrated(const double,const double,const double,
			  const double=1.,const int=-1);

    inline void GetSplittingParameters(double &,double &,double &,double &);

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_rms=ms; }

    inline ATOOLS::Mass_Selector *MS() const { return p_rms; }

    inline const SF_EEE_Map &FFFMap() const { return m_fffmap; }
    inline const SF_EEE_Map &FFIMap() const { return m_ffimap; }
    inline const SF_EEE_Map &IFFMap() const { return m_iffmap; }
    inline const SF_EEE_Map &IFIMap() const { return m_ifimap; }
    inline const SF_EEE_Map &FIFMap() const { return m_fifmap; }
    inline const SF_EEE_Map &FIIMap() const { return m_fiimap; }

    const SF_E_Map *HasKernel(const ATOOLS::Flavour &fli,
			      const ATOOLS::Flavour &flj,
			      const cstp::code type) const;

    int HasKernel(const ATOOLS::Flavour &fli,
                  const ATOOLS::Flavour &flj,
                  const ATOOLS::Flavour &flk,
                  const cstp::code type) const;

    inline double FSPT2Min() const { return m_k0sqf; }
    inline double ISPT2Min() const { return m_k0sqi; }

    inline void SetShower(Shower *const shower) { p_shower=shower; }
    inline double Weight() const { return m_weight; }

    inline void SetVariationWeights(ATOOLS::Variation_Weights * varwgts)
    { p_variationweights = varwgts; }

    inline void SetDisallowFlavour(const std::vector<size_t>& df)
    { m_disallowflav=df; }
  };
  
  void Sudakov::GetSplittingParameters
  (double & kt2,double & z,double & y,double & phi) {
    kt2 = m_kperp2; z = m_z; y = m_y; phi = m_phi;
  }
}

#endif
