#ifndef COMIX_Phasespace_PS_Generator_H
#define COMIX_Phasespace_PS_Generator_H

#include "COMIX/Amplitude/Amplitude.H"

#ifdef USING__Threading 
#include <pthread.h> 
#endif 

namespace PHASIC { class Cut_Data; }

namespace COMIX {

  class Process_Base;

  typedef std::vector<double> Double_Vector;

#ifdef USING__Threading
  class PS_Generator;

  struct CDBG_PG_TID {
    pthread_t m_id;
    PS_Generator *p_psg;
    size_t m_s, m_n, m_b, m_e, m_i;
    pthread_mutex_t m_s_mtx, m_t_mtx;
    pthread_cond_t m_s_cnd, m_t_cnd;
    CDBG_PG_TID(PS_Generator *const psg): 
      p_psg(psg), m_s(2), m_b(0), m_e(0) {}
  };// end of struct CDBG_PG_TID

  typedef std::vector<CDBG_PG_TID*> CDBG_PG_TID_Vector; 
#endif

  class PS_Generator {
  public:

    typedef std::map<Current_Base*,Current_Base*>      CB_Map;
    typedef std::multimap<Current_Base*,Current_Base*> CB_MMap;
    typedef std::pair<Current_Base*,Current_Base*>     CB_Pair;

    typedef std::map<size_t,Current_Vector> TCC_Map;

    void AddSTCC();

  private:

    Process_Base *p_xs;

    Int_Matrix m_cl;

    Current_Matrix m_cur;
    Current_Vector m_ctt;

    size_t m_n, m_itmin, m_itmax;

    int m_zmode, m_ecmode, m_pmsinit;

    double m_thmass, m_chmass;

    Double_Vector m_smasses, m_swidths;

    CB_MMap m_cmap;
    CB_Map  m_cbmap;

    TCC_Map m_tccs;

#ifdef USING__Threading 
    CDBG_PG_TID_Vector m_cts;

    static void *TCalcJL(void *arg);
#endif 

    void CleanUp();
    void CalcJL();

    bool AddCurrent(Current_Base *const ref,const ATOOLS::Flavour &fl,
		    const size_t &n,const int mode=0,
		    const double &m=0.0,const double &w=0.0,
		    Current_Base *const scc=NULL);
    void AddExtraCurrent(Current_Base *const cur,const size_t &n,
			 const double &m,const double &w,
			 Current_Base *const scc=NULL);

  public:
    
    // constructor
    PS_Generator(Process_Base *const xs=NULL);

    // destructor
    ~PS_Generator();

    // member functions
    bool Construct(Amplitude *const ampl);
    void SetPrefMasses(PHASIC::Cut_Data *const cuts);

    bool Evaluate();

    void SetColors(const Int_Vector &rc,
		   const Int_Vector &ac);

    // inline functions
    inline const Current_Matrix &Graphs() const { return m_cur; }

    inline void SetZMode(const int &mode) { m_zmode=mode; }

    inline double ThresholdMass() const { return m_thmass; }

    inline const Double_Vector &ISRMasses() const { return m_smasses; }
    inline const Double_Vector &ISRWidths() const { return m_swidths; }

  };// end of class PS_Generator

}// end of namespace COMIX

#endif
