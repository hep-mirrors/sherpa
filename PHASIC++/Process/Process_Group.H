#ifndef PHASIC_Process_Process_Group_H
#define PHASIC_Process_Process_Group_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace PHASIC {

  class Process_Group: public Process_Base {
  protected:

    std::vector<Process_Base*>          m_procs;
    std::map<std::string,Process_Base*> m_procmap;

    std::vector<int> m_blocks;

    Process_Base *GetProcess(const std::string &name);

    bool CheckFlavours(const Subprocess_Info &ii,
		       const Subprocess_Info &fi,int mode=0) const;
    void SetFlavour(Subprocess_Info &cii,Subprocess_Info &cfi,
		    const ATOOLS::Flavour &fl,const size_t i) const;
    bool ConstructProcesses(Process_Info &pi,const size_t &ci);
    bool ConstructProcess(Process_Info &pi);

  public:

    ~Process_Group();

    size_t Size() const override;
    Process_Base *operator[](const size_t &i) override;

    virtual void Add(Process_Base *const proc,const int mode=0);
    virtual bool Remove(Process_Base *const proc);
    virtual bool Delete(Process_Base *const proc);
    virtual void Clear();
    
    virtual Process_Base *GetProcess(const PHASIC::Process_Info &pi) const = 0;
    virtual bool Initialize(Process_Base *const proc) = 0;

    void ConstructColorMatrix() override;

    void DeSelect();
    bool SelectOne(const int mode);

    ATOOLS::Weight_Info *OneEvent(
        const int wmode,
        ATOOLS::Variations_Mode varmode=ATOOLS::Variations_Mode::all,
        const int mode=0) override;

    ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector& p,
                                     ATOOLS::Variations_Mode) override;

    bool InitScale() override;

    void SetScale(const Scale_Setter_Arguments &args) override;
    void SetKFactor(const KFactor_Setter_Arguments &args) override;
    void SetEWSudakovKFactor(const KFactor_Setter_Arguments &args) override;

    void SetFixedScale(const std::vector<double> &s) override;
    void SetSelectorOn(const bool on) override;
    void SetUseBIWeight(bool on) override;

    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create) override;
    void SetLookUp(const bool lookup) override;

    bool IsGroup() const override;
    
    bool ConstructProcesses();

    void SetGenerator(ME_Generator_Base *const gen) override;
    void SetShower(PDF::Shower_Base *const ps) override;
    void SetNLOMC(PDF::NLOMC_Base *const mc) override;
    void SetSelector(const Selector_Key &key) override;
    bool Trigger(const ATOOLS::Vec4D_Vector &p) override;

    void FillOnshellConditions() override;

    int PerformTests() override;

  };// end of class Process_Group

}// end of namespace PHASIC

#endif
