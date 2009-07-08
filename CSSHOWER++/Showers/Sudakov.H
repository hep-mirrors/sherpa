#ifndef CSSHOWER_Showers_Sudakov_H
#define CSSHOWER_Showers_Sudakov_H

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PDF/Main/PDF_Base.H"
#include "CSSHOWER++/Showers/Splitting_Function_Group.H"
#include "CSSHOWER++/Tools/Parton.H"
#include "MODEL/Main/Running_AlphaS.H"

namespace CSSHOWER {
  typedef std::map<ATOOLS::Flavour,Splitting_Function_Base*> SF_E_Map;
  typedef std::map<ATOOLS::Flavour,SF_E_Map> SF_EE_Map;
  typedef std::map<ATOOLS::Flavour,SF_EE_Map> SF_EEE_Map;

  class Shower;

  class Sudakov : public Splitting_Function_Group {
    double m_k0sq, m_nspect, m_phi;
    double m_deltaz, m_kperp2, m_z, m_y, m_ymin, m_ymax, m_scale;
    double m_as_is_fac,m_as_fs_fac,m_rbmax,m_zmin,m_zmax;
    PDF::PDF_Base         *  p_pdf[2];
    ATOOLS::Mass_Selector *  p_rms;
    ATOOLS::Flavour m_cfl, m_flspec;
    Parton * p_split, * p_spec;
    SF_EEE_Map m_ffmap, m_fimap, m_ifmap, m_iimap;
    SF_EEE_Map m_sffmap, m_sfimap, m_sifmap, m_siimap;
    std::vector<Splitting_Function_Base *> m_addsplittings;
    Shower *p_shower;
    void ProduceT();
    bool Veto(double,double);
    bool KinCheck(double,double);
    bool Splitting(double,double);
    bool DefineFFBoundaries(double,double);
    bool DefineFIBoundaries(double,double,int);
    bool DefineIFBoundaries(double,double,int);
    bool DefineIIBoundaries(double,double,int);
    void Add(Splitting_Function_Base * split);
    void AddToMaps(Splitting_Function_Base * split,const int mode=1);
  public:
    Sudakov(PDF::ISR_Handler *isr,const int qed);
    ~Sudakov();

    void InitSplittingFunctions(MODEL::Model_Base *md);
    void SetCoupling(MODEL::Model_Base *md,const double &k0sq,
		     const double &isfac,const double &fsfac);
    bool   Dice(Parton *);
    double OverIntegrated(const double,const double,const double,
			  const double=1.,const int=-1);

    inline void GetSplittingParameters(double &,double &,double &,double &);
    inline void SetPT2Min(const double &pt2) {m_k0sq=pt2;}     
    inline double FSCplFac() const { return m_as_fs_fac; }
    inline double ISCplFac() const { return m_as_is_fac; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_rms=ms; }
    inline void SetRBMax(const double &rbmax) { m_rbmax=rbmax; }

    inline ATOOLS::Mass_Selector *MS() const { return p_rms; }

    inline const SF_EEE_Map &FFMap() const { return m_ffmap; }
    inline const SF_EEE_Map &FIMap() const { return m_fimap; }
    inline const SF_EEE_Map &IFMap() const { return m_ifmap; }
    inline const SF_EEE_Map &IIMap() const { return m_iimap; }

    bool HasKernel(const ATOOLS::Flavour &fli,
		   const ATOOLS::Flavour &flj,
		   const cstp::code type) const;

    inline double PT2Min() const { return m_k0sq; }

    inline void SetShower(Shower *const shower) { p_shower=shower; }
  };
  
  void Sudakov::GetSplittingParameters
  (double & kt2,double & z,double & y,double & phi) {
    kt2 = m_kperp2; z = m_z; y = m_y; phi = m_phi;
  }
}

#endif
