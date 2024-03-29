#ifndef COMIX_Amplitude_Amplitude_H
#define COMIX_Amplitude_Amplitude_H

#include "PHASIC++/Process/Subprocess_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Coupling_Data.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__Threading 
#include <pthread.h> 
#endif 

namespace PHASIC {
  class Color_Integrator;
  class Helicity_Integrator;
}

namespace METOOLS {
  class Dipole_Info;
  class Spin_Amplitudes;
}

namespace PDF {
  class NLOMC_Base;
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

  typedef std::map<size_t,size_t> SizeT_Map;

#ifdef USING__Threading
  class Amplitude;

  struct CDBG_ME_TID {
    pthread_t m_id;
    Amplitude *p_ampl;
    size_t m_s, m_n, m_b, m_e, m_i;
    pthread_mutex_t m_s_mtx, m_t_mtx;
    pthread_cond_t m_s_cnd, m_t_cnd;
    CDBG_ME_TID(): p_ampl(NULL), m_s(2), m_b(0), m_e(0) {}
  };// end of struct CDBG_ME_TID

  typedef std::vector<CDBG_ME_TID*> CDBG_ME_TID_Vector; 
#endif
 
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

    size_t m_nin, m_nout, m_dec;
    size_t m_n, m_wfmode, m_pgmode, m_ngpl;
    bool m_sccmur, m_murcoeffvirt;
    size_t m_minntc, m_maxntc, m_stype;
    char   m_pmode;

    std::vector<int> m_maxcpl, m_mincpl, m_maxacpl, m_minacpl;

    std::vector<Spin_Structure<DComplex> > m_ress, m_cress;
    std::vector<std::pair<size_t,size_t> > m_on, m_son;

    std::vector<std::vector<double> >     m_dsij;
    std::map<size_t,std::pair<int,int> >  m_dsm;

    std::vector<std::pair<size_t,double> > m_dsf;

    double m_res, m_born, m_cmur[2];
    double m_smth, m_smpow;

    Current_Matrix m_cur;
    Current_Vector m_scur;

    METOOLS::Dipole_Info *p_dinfo;

    PHASIC::Color_Integrator    *p_colint;
    PHASIC::Helicity_Integrator *p_helint;

    double m_sf, m_fsf;
    bool   m_trig;

    ATOOLS::NLO_subevtlist    m_subs;
    PHASIC::Virtual_ME2_Base *p_loop;  

    std::vector<std::map<size_t,std::vector<long int> > > m_affm;

    ATOOLS::NLO_subevt *p_sub;

    size_t BornID(const size_t &id,const NLO_subevt *sub) const;

#ifdef USING__Threading 
    CDBG_ME_TID_Vector *p_cts;
#endif 

    void CleanUp();
    void Prune();
    void ConstructNLOEvents();
    void ConstructDSijMap();

    int  CheckDecay(const ATOOLS::Flavour &fl,const Int_Vector &ids) const;
    bool MatchDecay(const Vertex_Key &vkey) const;
    Vertex *AddCurrent(const Current_Key &ckey,Vertex_Key &vkey,
		       const size_t &n,const int dec,
		       std::vector<int> &maxcpl,std::vector<int> &mincpl,
		       std::map<std::string,Current*> &curs);
    void AddCurrent(const Int_Vector &ids,const size_t &n,
		    const ATOOLS::Flavour &fl,const int dir);

    Current *CopyCurrent(Current *const c);
    bool AddRSDipole(Current *const c,Current *const sc,
		     Current_Vector &scur,int stype);
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

    bool ConstructCouplings(MODEL::Coupling_Map *const cpls);
    bool CheckOrders();
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
		    My_In_File &ampfile,
		    MODEL::Model_Base *const model,
		    MODEL::Coupling_Map *const cpls,
		    const int stype,const int smode,
		    const ATOOLS::cs_itype::type itype,
		    const std::vector<int> &maxcpl,
		    const std::vector<int> &mincpl,
		    const std::vector<int> &maxacpl,
		    const std::vector<int> &minacpl,
		    const size_t &minntc,const size_t &maxntc,
		    const std::string &name);

    void ResetJ();
    void ResetZero();

    bool Evaluate(const Int_Vector &chirs);
    bool EvaluateAll(const bool &mode=false);

    double Differential(NLO_subevt *const sub=NULL);
    double Differential(const Int_Vector &ci,const Int_Vector &cj,
			const int set=0);

    bool Map(const Amplitude &ampl,ATOOLS::Flavour_Map &flmap);
    bool GaugeTest(const ATOOLS::Vec4D_Vector &moms,const int mode=0);

    bool SetMomenta(const ATOOLS::Vec4D_Vector &moms);
    void SetColors(const Int_Vector &rc,
		   const Int_Vector &ac,const int set=0);

    bool RSTrigger(PHASIC::Combined_Selector *const sel,const int mode);

    void FillCombinations(Combination_Set &combs,CFlavVector_Map &flavs,
			  SizeT_Map *brs=NULL,const NLO_subevt *sub=NULL) const;
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes> &amps,
			std::vector<std::vector<Complex> > &cols);

    void WriteOutGraphs(const std::string &file) const;

    void WriteOutAmpFile(const std::string &name);
    bool ReadInAmpFile(const std::string &name,My_In_File &ampfile);

    void PrintStatistics(std::ostream &str,const int mode=0) const;

    double Coupling(const int mode) const;
    void   FillMEWeights(ATOOLS::ME_Weight_Info &wgtinfo) const;

    double KT2Trigger(ATOOLS::NLO_subevt *const sub,const int mode);
    void   SetNLOMC(PDF::NLOMC_Base *const mc);

    void SetCTS(void *const cts);

#ifdef USING__Threading 
    static void *TCalcJL(void *arg);
#endif 

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

    inline const std::vector<int> &MaxCpl() const { return m_maxcpl; }
    inline const std::vector<int> &MinCpl() const { return m_mincpl; }

    inline void SetMinNTChannel(const size_t &ntc) { m_minntc=ntc; }

    inline const ATOOLS::Vec4D_Vector   &Momenta() const  { return m_p;  }
    inline const ATOOLS::Flavour_Vector &Flavours() const { return m_fl; }

    inline const Current_Matrix &Currents() const { return m_cur; }

    inline void SetDecayInfos(const ATOOLS::DecayInfo_Vector &id) 
    { m_decid=id; m_dec=m_decid.size(); }
    inline const ATOOLS::DecayInfo_Vector &DecayInfos() const
    { return m_decid; }

    inline void SetNoDecays(const ATOOLS::Flavour_Vector &id) { m_ndc=id; }

    inline char PMode() const { return m_pmode; }

    inline ATOOLS::NLO_subevtlist &SubEvts() { return m_subs; }

    inline void SetLoopME(PHASIC::Virtual_ME2_Base *const me) { p_loop=me; }

    inline METOOLS::Dipole_Info *DInfo() const { return p_dinfo; }

    inline const std::vector<std::vector<double> > &DSij() const
    { return m_dsij; }

    inline double Born() const { return m_born; }

  };// end of class Amplitude

}// end of namespace COMIX

#endif
