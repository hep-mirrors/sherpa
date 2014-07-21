#ifndef COMIX_Amplitude_Amplitude_H
#define COMIX_Amplitude_Amplitude_H

#include "PHASIC++/Process/Subprocess_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Coupling_Data.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace PHASIC {
  class Color_Integrator;
  class Helicity_Integrator;
}

namespace METOOLS {
  class Dipole_Info;
  class Spin_Amplitudes;
}

using namespace METOOLS;

namespace COMIX {

  struct Coupling_Info {
    METOOLS::Vertex *p_v;
    size_t m_oqcd, m_oew;
    MODEL::Coupling_Data *p_aqcd, *p_aqed;
    inline Coupling_Info(METOOLS::Vertex *const v,
			 const size_t &oqcd,const size_t &oew,
			 MODEL::Coupling_Data *const aqcd,
			 MODEL::Coupling_Data *const aqed):
      p_v(v), m_oqcd(oqcd), m_oew(oew), p_aqcd(aqcd), p_aqed(aqed) {}
  };// end of struct Coupling_Info

  typedef std::vector<Coupling_Info> CouplingInfo_Vector;

  typedef std::vector<Current_Vector> Current_Matrix;

  typedef std::vector<long unsigned int> LongInt_Vector;

  typedef std::set<std::pair<size_t,size_t> >     Combination_Set;
  typedef std::map<size_t,ATOOLS::Flavour_Vector> CFlavVector_Map;

  class Amplitude {
  public:

    typedef std::complex<double> DComplex;

    typedef std::vector<DComplex> DComplex_Vector;

  private:

    MODEL::Model_Base *p_model;

    ATOOLS::Flavour_Vector m_fl, m_ndc;
    ATOOLS::Vec4D_Vector   m_p;

    ATOOLS::DecayInfo_Vector m_decid;

    CouplingInfo_Vector m_cpls;

    Int_Vector   m_ch, m_dirs, m_sid, m_cchirs;
    Int_Matrix   m_cl;

    size_t m_nin, m_nout;
    size_t m_n, m_wfmode, m_pgmode, m_ngpl, m_sccmur;
    size_t m_oew, m_oqcd, m_maxoew, m_maxoqcd;
    size_t m_minntc, m_maxntc;
    char   m_pmode;

    Spin_Structure<DComplex> m_ress, m_rress, m_cress;

    std::vector<std::vector<double> >     m_dsij;
    std::map<size_t,std::pair<int,int> >  m_dsm;

    std::vector<double> m_dsf;

    double m_res, m_cmur[2];
    double m_smth, m_smpow;

    Current_Matrix m_cur;
    Current_Vector m_scur;

    Combination_Set m_combs;
    CFlavVector_Map m_flavs;

    METOOLS::Dipole_Info *p_dinfo;

    PHASIC::Color_Integrator    *p_colint;
    PHASIC::Helicity_Integrator *p_helint;

    double m_sf, m_fsf;
    bool   m_trig;

    ATOOLS::NLO_subevtlist    m_subs;
    PHASIC::Virtual_ME2_Base *p_loop;  

    std::vector<std::vector<long int> > m_affm;

    ATOOLS::NLO_subevt *p_sub;

    static size_t     MakeId(const Int_Vector &ids,const int t);
    static Int_Vector MakeId(const size_t &id,const size_t &n);

    void CleanUp();
    void FillCombinations();
    void ConstructNLOEvents();
    void ConstructDSijMap();

    int  CheckDecay(const ATOOLS::Flavour &fl,const Int_Vector &ids) const;
    bool MatchIndices(const Int_Vector &ids,const Vertex_Key &vkey) const;
    Vertex *AddCurrent(const Current_Key &ckey,Vertex_Key &vkey,
		       const size_t &n,const int dec,
		       size_t &oewmax,size_t &oqcdmax,
		       std::map<std::string,Current*> &curs);
    void AddCurrent(const Int_Vector &ids,const size_t &n,
		    const ATOOLS::Flavour &fl,const int dir);

    Current *CopyCurrent(Current *const c);
    bool AddRSDipole(Current *const c,Current *const sc,Current_Vector &scur);
    bool AddRSDipoles();
    bool AddVIDipole(Current *const c,Current *const sc,Current_Vector &scur);
    bool AddVIDipoles();

    bool Construct(ATOOLS::Flavour_Vector &fls,
		   Int_Vector ids,const size_t &n);
    bool Construct(const ATOOLS::Flavour_Vector &flavs);
    bool Construct(const Int_Vector &incs,
		   const ATOOLS::Flavour_Vector &flavs,
		   MODEL::Model_Base *const model,
		   MODEL::Coupling_Map *const cpls);

    bool ConstructChirs();

    void SetGauge(const size_t &n);

    void CalcJL();

    void SetCouplings() const;

    double EpsSchemeFactor(const ATOOLS::Vec4D_Vector &mom) const;

    void WriteOutGraph(std::ostream &str,Graph_Node *graph,
		       size_t &ng,std::set<std::string> &cvs) const;

  public:
    
    // constructor
    Amplitude();

    // destructor
    ~Amplitude();

    // member functions
    bool Initialize(const size_t &nin,const size_t &nout,
		    const std::vector<ATOOLS::Flavour> &flavs,
		    const double &isf,const double &fsf,
		    MODEL::Model_Base *const model,
		    MODEL::Coupling_Map *const cpls,const int smode,
		    const size_t &oew,const size_t &oqcd,
		    const size_t &maxoew,const size_t &maxoqcd,
		    const size_t &minntc,const size_t &maxntc,
		    const std::string &name);

    void ResetZero();

    bool Evaluate(const Int_Vector &chirs);
    bool EvaluateAll();

    double Differential(NLO_subevt *const sub=NULL);
    double Differential(const Int_Vector &ci,const Int_Vector &cj,
			const int set=0);

    bool Map(const Amplitude &ampl,ATOOLS::Flavour_Map &flmap);
    bool GaugeTest(const ATOOLS::Vec4D_Vector &moms,const int mode=0);

    bool SetMomenta(const ATOOLS::Vec4D_Vector &moms);
    void SetColors(const Int_Vector &rc,
		   const Int_Vector &ac,const int set=0);

    bool JetTrigger(PHASIC::Combined_Selector *const sel,const int mode);

    bool Combinable(const size_t &idi,const size_t &idj) const;

    const ATOOLS::Flavour_Vector &
    CombinedFlavour(const size_t &idij) const;

    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes> &amps,
			std::vector<std::vector<Complex> > &cols);

    void WriteOutGraphs(const std::string &file) const;

    void WriteOutAmpFile(const std::string &name);
    bool ReadInAmpFile(const std::string &name);

    void PrintStatistics(std::ostream &str,const int mode=0) const;

    double Coupling(const int mode) const;
    void   FillMEWeights(ATOOLS::ME_wgtinfo &wgtinfo) const;

    double KT2Trigger(ATOOLS::NLO_subevt *const sub,const int mode);

    // inline functions
    inline void SetColorIntegrator(PHASIC::Color_Integrator *const cint)
    { p_colint=cint; }
    inline void SetHelicityIntegrator(PHASIC::Helicity_Integrator *const hint)
    { p_helint=hint; }

    inline PHASIC::Color_Integrator *ColorIntegrator() const
    { return p_colint; }
    inline PHASIC::Helicity_Integrator *HelicityIntegrator() const
    { return p_helint; }

    inline double FSSymmetryFactor() const { return m_fsf;      }
    inline double ISSymmetryFactor() const { return m_sf/m_fsf; }

    inline size_t MaxOrderEW() const  { return m_maxoew;  }
    inline size_t MaxOrderQCD() const { return m_maxoqcd; }

    inline void SetMinNTChannel(const size_t &ntc) { m_minntc=ntc; }

    inline const ATOOLS::Vec4D_Vector   &Momenta() const  { return m_p;  }
    inline const ATOOLS::Flavour_Vector &Flavours() const { return m_fl; }

    inline const Current_Matrix &Currents() const { return m_cur; }

    inline void SetDecayInfos(const ATOOLS::DecayInfo_Vector &id) 
    { m_decid=id; }
    inline const ATOOLS::DecayInfo_Vector &DecayInfos() const
    { return m_decid; }

    inline void SetNoDecays(const ATOOLS::Flavour_Vector &id) { m_ndc=id; }

    inline char PMode() const { return m_pmode; }

    inline ATOOLS::NLO_subevtlist &SubEvts() { return m_subs; }

    inline void SetLoopME(PHASIC::Virtual_ME2_Base *const me) { p_loop=me; }

    inline METOOLS::Dipole_Info *DInfo() const { return p_dinfo; }

    inline const std::vector<std::vector<double> > &DSij() const
    { return m_dsij; }

  };// end of class Amplitude

}// end of namespace COMIX

#endif
