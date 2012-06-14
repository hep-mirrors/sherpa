#ifndef AMEGIC_DipoleSubtraction_Single_Real_Correction_H
#define AMEGIC_DipoleSubtraction_Single_Real_Correction_H

#include "AMEGIC++/Main/Process_Base.H"

#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
#include "AMEGIC++/DipoleSubtraction/Single_OSTerm.H"

#include "AMEGIC++/Amplitude/Amplitude_Handler.H"
#include <vector>

namespace AMEGIC {
  class Single_DipoleTerm;
  class Helicity;
  struct DM_Info {
    ATOOLS::Flavour_Vector *p_fl;
    std::string m_pname;
    DM_Info(ATOOLS::Flavour_Vector *const fl=NULL,
	    const std::string &pname=""):
      p_fl(fl), m_pname(pname) {}
  };
  class Single_Real_Correction : public Process_Base {
  private:
    bool                    m_newlib;
    double                  m_iresult;
    int                     m_libnumb;
    bool                    m_valid;
    bool                    m_ossubon;

    Single_Real_Correction     * p_partner;

    Process_Base          * p_tree_process;
    ATOOLS::Vec4D_Vector    m_real_momenta;
    std::vector<Single_DipoleTerm*> m_subtermlist;
    std::vector<Single_OSTerm*> m_subostermlist;
    ATOOLS::NLO_subevtlist          m_subevtlist;
    ATOOLS::NLO_subevt              m_realevt;

    std::vector<size_t> m_sids;
    std::map<void*,DM_Info> m_dfmap;
    void ReMapFlavs(ATOOLS::NLO_subevt *const sub,const int mode=0);
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Real_Correction();
    ~Single_Real_Correction();
 
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
    bool Combinable(const size_t &idi,const size_t &idj);
    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

    void SetFixedScale(const std::vector<double> &s);
    void SetSelectorOn(const bool on);
    void SetGenerator(PHASIC::ME_Generator_Base *const gen);
    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return m_newlib;}
    int                 InitAmplitude(MODEL::Model_Base *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    Amplitude_Handler * GetAmplitudeHandler()          { return p_tree_process->GetAmplitudeHandler();}
    Helicity *          GetHelicity()                  { return p_tree_process->GetHelicity(); }    
    double              Result()                       { return m_iresult; } 
    bool                IsValid()                      { return m_valid; }

    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    void             SetLookUp(const bool lookup);
    std::string      LibName()                          { return p_tree_process->LibName();     }
    std::string      PSLibName()                        { return p_tree_process->PSLibName();   }
    Process_Base   * Partner() const                    { return p_partner;     }
    void             Minimize();

    void SetSelector(const PHASIC::Selector_Key &key);
    void SetShower(PDF::Shower_Base *const ps);

    inline ATOOLS::NLO_subevtlist*             GetSubevtList()
    { return &m_subevtlist; }
    inline Single_DipoleTerm*                  GetSubTerm(size_t i)
    { return m_subtermlist[i]; }
    inline size_t                              GetSubTermNumber()
    { return m_subtermlist.size(); }

    int Type() { return 100; }
    
    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double         Partonic(const ATOOLS::Vec4D_Vector &,const int mode);
    double         operator()(const ATOOLS::Vec4D_Vector &,const int mode);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);

    void SetScale(const PHASIC::Scale_Setter_Arguments &args);
    void SetKFactor(const PHASIC::KFactor_Setter_Arguments &args);

    void           FillAmplitudes(METOOLS::Amplitude_Tensor*, double =1.);
    int            NumberOfDiagrams();
    Point        * Diagram(int i);

    bool PerformTests() { return true; }

   /*------------------------------------------------------------------------------

      Helpers
      
      ------------------------------------------------------------------------------*/
  public:
    void           PrintProcessSummary(int=0); 
    void           PrintSubevtSummary();
    size_t         NumberOfDipoles() {return m_subtermlist.size();} 
  };
}



#endif

