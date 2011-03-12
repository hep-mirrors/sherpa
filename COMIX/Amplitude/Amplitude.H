#ifndef COMIX_Amplitude_Amplitude_H
#define COMIX_Amplitude_Amplitude_H

#include "ATOOLS/Phys/Spinor.H"
#include "METOOLS/Main/C_Vector.H"
#include "METOOLS/Main/C_Spinor.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "COMIX/Amplitude/Vertex_Base.H"
#include "MODEL/Main/Coupling_Data.H"
#include "METOOLS/Main/C_Scalar.H"
#include "METOOLS/Main/C_Spinor.H"
#include "METOOLS/Main/C_Vector.H"
#include "METOOLS/Main/C_Tensor.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__Threading 
#include <pthread.h> 
#endif 

namespace COMIX {

  typedef std::vector<Current_Base*>  Current_Vector;
  typedef std::vector<Current_Vector> Current_Matrix;

  typedef std::vector<ATOOLS::Vec4D> Vec4D_Vector;

  typedef std::vector<long unsigned int> LongInt_Vector;

  typedef std::set<std::pair<size_t,size_t> >     Combination_Set;
  typedef std::map<size_t,ATOOLS::Flavour_Vector> CFlavVector_Map;

#ifdef USING__Threading
  class Amplitude;

  struct CDBG_ME_TID {
    pthread_t m_id;
    Amplitude *p_ampl;
    size_t m_s, m_n, m_b, m_e, m_i;
    pthread_mutex_t m_s_mtx, m_t_mtx;
    pthread_cond_t m_s_cnd, m_t_cnd;
    CDBG_ME_TID(Amplitude *const ampl): 
      p_ampl(ampl), m_s(2), m_b(0), m_e(0) {}
  };// end of struct CDBG_ME_TID

  typedef std::vector<CDBG_ME_TID*> CDBG_ME_TID_Vector; 
#endif
 
  class Amplitude {
  public:

    typedef std::complex<double>      DComplex;
    typedef std::complex<long double> QComplex;

    typedef std::vector<DComplex> DComplex_Vector;
    typedef std::vector<QComplex> QComplex_Vector;

  private:

    Model *p_model;

    MODEL::Coupling_Data *p_aqcd, *p_aqed;

    ATOOLS::Flavour_Vector m_fl;

    ATOOLS::DecayInfo_Vector m_decid;

    Vec4D_Vector m_p;
    Int_Vector   m_ch, m_dirs; 
    Int_Matrix   m_cl, m_chirs;

    size_t m_n, m_nf, m_pgmode, m_ngpl;
    size_t m_oew, m_oqcd, m_maxoew, m_maxoqcd, m_minntc;
    char   m_pmode;

    QComplex_Vector m_ress;

    double m_res;

    Current_Matrix m_cur;

    LongUIntMap_Matrix m_cchirs;

    Combination_Set m_combs;
    CFlavVector_Map m_flavs;

#ifdef USING__Threading 
    CDBG_ME_TID_Vector m_cts;

    static void *TCalcJL(void *arg);
#endif 

    static size_t     MakeId(const Int_Vector &ids,const int t);
    static Int_Vector MakeId(const size_t &id,const size_t &n);

    void CleanUp();
    void FillCombinations();

    int  CheckDecay(const ATOOLS::Flavour &fl,const Int_Vector &ids) const;
    bool MatchIndices(const Int_Vector &ids,const size_t &n,
		      const size_t &i,const size_t &j,const size_t &k);
    void AddCurrent(const Int_Vector &ids,const size_t &n,
		    const ATOOLS::Flavour &fl,const int dir);

    bool Construct(ATOOLS::Flavour_Vector &fls,
		   Int_Vector ids,const size_t &n);
    bool Construct(const ATOOLS::Flavour_Vector &flavs);

    bool CheckChirs(const Int_Vector &chirs);
    bool ConstructChirs(Int_Vector chirs,const size_t &i);

    void SetGauge(const size_t &n);

    void CalcJL();

    double CouplingFactor(Vertex_Base *const v) const;

    void WriteOutGraph(std::ostream &str,Graph_Node *graph,
		       size_t &ng,std::set<std::string> &cvs) const;

  public:
    
    // constructor
    Amplitude();

    // destructor
    ~Amplitude();

    // member functions
    bool Evaluate(const Int_Vector &chirs);
    bool EvaluateAll();

    void ResetZero();

    bool Construct(const Int_Vector &incs,
		   const ATOOLS::Flavour_Vector &flavs,
		   Model *const model,MODEL::Coupling_Map *const cpls);
    bool Map(const Amplitude &ampl,ATOOLS::Flavour_Map &flmap);
    bool GaugeTest(const Vec4D_Vector &moms);

    void SetMomenta(const Vec4D_Vector &moms);
    void SetColors(const Int_Vector &rc,
		   const Int_Vector &ac,const bool set=false);

    bool Combinable(const size_t &idi,const size_t &idj) const;

    const ATOOLS::Flavour_Vector &
    CombinedFlavour(const size_t &idij) const;

    void WriteOutGraphs(const std::string &file) const;

    void PrintStatistics(std::ostream &str,const int mode=0) const;

    // inline functions
    inline void SetOrderEW(const size_t &oew)   { m_oew=oew;   }
    inline void SetOrderQCD(const size_t &oqcd) { m_oqcd=oqcd; }

    inline void SetMaxOrderEW(const size_t &oew)   { m_maxoew=oew;   }
    inline void SetMaxOrderQCD(const size_t &oqcd) { m_maxoqcd=oqcd; }

    inline size_t MaxOrderEW() const  { return m_maxoew;  }
    inline size_t MaxOrderQCD() const { return m_maxoqcd; }

    inline void SetMinNTChannel(const size_t &ntc) { m_minntc=ntc; }

    inline const Int_Matrix &Chiralities() const { return m_chirs; }

    inline double Result() const { return m_res; }

    inline const Vec4D_Vector           &Momenta() const  { return m_p;  }
    inline const ATOOLS::Flavour_Vector &Flavours() const { return m_fl; }

    inline const Current_Matrix &Currents() const { return m_cur; }

    inline void SetDecayInfos(const ATOOLS::DecayInfo_Vector &id) 
    { m_decid=id; }
    inline const ATOOLS::DecayInfo_Vector &DecayInfos() const
    { return m_decid; }

    inline char PMode() const { return m_pmode; }

  };// end of class Amplitude

}// end of namespace COMIX

#endif
