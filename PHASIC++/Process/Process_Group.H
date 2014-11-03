#ifndef PHASIC_Process_Process_Group_H
#define PHASIC_Process_Process_Group_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace PHASIC {

  class Process_Group: public Process_Base {
  protected:

    std::vector<Process_Base*>          m_procs;
    std::map<std::string,Process_Base*> m_procmap;

    Process_Base *GetProcess(const std::string &name);

    bool CheckFlavours(const Subprocess_Info &ii,
		       const Subprocess_Info &fi,int mode=0) const;
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
    bool SelectOne(const int mode);

    Weight_Info *OneEvent(const int wmode,const int mode=0);

    double Differential(const ATOOLS::Vec4D_Vector &p);

    void SetScale(const Scale_Setter_Arguments &args);
    void SetKFactor(const KFactor_Setter_Arguments &args);

    void SetFixedScale(const std::vector<double> &s);
    void SetSelectorOn(const bool on);
    void SetUseBIWeight(bool on);

    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create);
    void SetLookUp(const bool lookup);

    bool IsGroup() const;
    
    bool ConstructProcesses();

    void SetGenerator(ME_Generator_Base *const gen);
    void SetShower(PDF::Shower_Base *const ps);
    void SetSelector(const Selector_Key &key);
    bool Trigger(const ATOOLS::Vec4D_Vector &p);

    void FillOnshellConditions();

  };// end of class Process_Group

}// end of namespace ATOOLS

#endif
