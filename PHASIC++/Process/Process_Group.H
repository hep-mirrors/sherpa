#ifndef PHASIC_Process_Process_Group_H
#define PHASIC_Process_Process_Group_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__Threading
#include <pthread.h> 
#endif 

namespace PHASIC {

#ifdef USING__Threading
  class Process_Group;

  struct PG_TID {
    pthread_t m_id;
    Process_Group *p_proc;
    double m_d;
    const ATOOLS::Vec4D_Vector *p_p;
    size_t m_s, m_m, m_b, m_e, m_i;
    pthread_mutex_t m_s_mtx, m_t_mtx;
    pthread_cond_t m_s_cnd, m_t_cnd;
    PG_TID(Process_Group *const proc): 
      p_proc(proc), m_s(2), m_b(0), m_e(0) {}
  };// end of struct PS_TID

  typedef std::vector<PG_TID*> PG_TID_Vector; 
#endif

  class Process_Group: public Process_Base {
  protected:

    std::vector<Process_Base*>          m_procs;
    std::map<std::string,Process_Base*> m_procmap;

#ifdef USING__Threading 
    std::vector<Process_Base *> m_umprocs, m_mprocs;
    PG_TID_Vector m_cts;

    static void *TDifferential(void *arg);
#endif 

    Process_Base *GetProcess(const std::string &name);

    bool CheckFlavours(const Subprocess_Info &ii,const Subprocess_Info &fi) const;
    void SetFlavour(Subprocess_Info &cii,Subprocess_Info &cfi,
		    const ATOOLS::Flavour &fl,const size_t i) const;
    bool ConstructProcesses(Process_Info pi,const size_t &ci);

  public:

    ~Process_Group();

    size_t Size() const;
    Process_Base *operator[](const size_t &i);

    virtual void Add(Process_Base *const proc);
    virtual bool Remove(Process_Base *const proc);
    virtual bool Delete(Process_Base *const proc);
    virtual void Clear();
    
    virtual Process_Base *GetProcess(const PHASIC::Process_Info &pi) const = 0;
    virtual bool Initialize(Process_Base *const proc) = 0;

    void DeSelect();
    bool SelectOne();

    Weight_Info *OneEvent();
    Weight_Info *WeightedEvent(const int mode=0);

    double Differential(const ATOOLS::Vec4D_Vector &p);
    double Differential2();

    void SetScale(const std::string &scale,
		  const std::string &mur2tag,
		  const std::string &muf2tag);
    void SetKFactor(const std::string &kfactor,
		    const size_t &oqcdlo,
		    const size_t &oewlo);

    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create);
    void SetLookUp(const bool lookup);

    bool IsGroup() const;
    
    bool ConstructProcesses();
    void SetUpThreading();

    void SetGenerator(ME_Generator_Base *const gen);
    void SetShower(PDF::Shower_Base *const ps);
    void SetSelector(const Selector_Key &key);

  };// end of class Process_Group

}// end of namespace ATOOLS

#endif
