#ifndef CSSHOWER_Showers_Sudakov_H
#define CSSHOWER_Showers_Sudakov_H

#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/PDF_Base.H"
#include "CSSHOWER++/Showers/Splitting_Function_Group.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Phys/Variations.H"


namespace CSSHOWER {
  typedef std::map<ATOOLS::Flavour,Splitting_Function_Base*> SF_E_Map;
  typedef std::map<ATOOLS::Flavour,SF_E_Map> SF_EE_Map;
  typedef std::map<ATOOLS::Flavour,SF_EE_Map> SF_EEE_Map;

  class Shower;

  class Sudakov : public Splitting_Function_Group {
    double m_k0sqi, m_k0sqf, m_k0sq_gsplit_fac, m_sphi;
    double m_deltaz, m_st, m_sz, m_sy, m_ymin, m_ymax, m_scale;
    double m_zmin,m_zmax, m_weight;
    int m_qcdmode,m_ewmode;
    ATOOLS::Flavour m_cfl, m_flspec;
    Parton * p_split, * p_spect;
    SF_EEE_Map m_fffmap, m_ffimap, m_iffmap, m_ifimap, m_fifmap, m_fiimap;
    SF_EEE_Map m_sffmap, m_sfimap, m_sifmap, m_siimap;
    std::vector<Splitting_Function_Base *> m_addsplittings;
    Shower *p_shower;
    std::vector<SFC_Getter*> m_cgets;
    bool m_keeprewinfo, m_forced_splittings;
    double m_reweightscalecutoff, m_gluon_xscaling_in_forced_splittings;
    std::pair<double, double> m_pdfmin;
    static bool s_init;
    double ProduceT(double t);
    bool Veto(double,double,double,double,double);
    bool Splitting(double,double,double,double,double);
    bool DefineFFBoundaries(double,double);
    bool DefineFIBoundaries(double,double,int);
    bool DefineIFBoundaries(double,double,int);
    bool DefineIIBoundaries(double,double,int);
    void Add(Splitting_Function_Base * split);
    void AddToMaps(Splitting_Function_Base * split,const int mode=1);
  public:
    Sudakov(PDF::ISR_Handler *isr,const int qcd,const int qed);
    ~Sudakov();

    void InitSplittingFunctions(MODEL::Model_Base *md,const int kfmode);
    void SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac,
		     const double &k0sq_gsplit_fac);
    void SetKeepReweightingInfo(bool on)  { m_keeprewinfo = on; }
    void SetReweightScaleCutoff(double v) { m_reweightscalecutoff = v; }
    void SetPDFMin(std::pair<double, double>& m) { m_pdfmin = m; }
    void SetForcedHQSplittings(bool on, const double &expo) {
      m_forced_splittings                   = on;
      m_gluon_xscaling_in_forced_splittings = expo;
    }
    bool   Generate(Parton *,double);
    int    Generate(Parton *,Parton *,double,double,double&,double&,double&,double&);
    double OverIntegrated(const double,const double,const double,
			  const double=1.,const int=-1);

    inline void GetSplittingParameters(double &,double &,double &,double &);

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
    double CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
                  const ATOOLS::Flavour &flk,const cstp::code type,
		  const int cpl,const double &mu2) const;

    inline double FSPT2Min() const { return m_k0sqf; }
    inline double ISPT2Min() const { return m_k0sqi; }

    inline void SetShower(Shower *const shower) { p_shower=shower; }
    inline double Weight() const { return m_weight; }
  };

  void Sudakov::GetSplittingParameters
  (double & kt2,double & z,double & y,double & phi) {
    kt2 = m_st; z = m_sz; y = m_sy; phi = m_sphi;
  }
}

#endif
