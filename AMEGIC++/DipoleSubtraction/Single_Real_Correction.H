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
    std::vector<size_t> *p_id;
    std::string m_pname;
    DM_Info(ATOOLS::Flavour_Vector *const fl=NULL,
	    std::vector<size_t> *const id=NULL,
	    const std::string &pname=""):
      p_fl(fl), p_id(id), m_pname(pname) {}
  };

  class Single_Real_Correction : public Process_Base {
  private:
    bool                    m_newlib, m_no_tree, m_listdips;
    double                  m_iresult;
    double                  m_smear_threshold, m_smear_power;
    int                     m_libnumb;
    bool                    m_ossubon;
    ATOOLS::sbt::subtype    m_user_stype;
    size_t                  m_pspisrecscheme,m_pspfsrecscheme;

    Single_Real_Correction     * p_partner;

    Process_Base          * p_tree_process;
    ATOOLS::Vec4D_Vector    m_real_momenta;
    std::vector<Single_DipoleTerm*> m_subtermlist;
    std::vector<Single_OSTerm*>     m_subostermlist;
    ATOOLS::NLO_subevtlist          m_subevtlist;
    ATOOLS::NLO_subevt              m_realevt;

    std::vector<size_t> m_sids;
    std::map<void*,DM_Info> m_dfmap;

    void FillProcessMap(PHASIC::NLOTypeStringProcessMap_Map *apmap);
    void ReMapFlavs(ATOOLS::NLO_subevt *const sub,const int mode=0);
    void SmearCounterEvents(ATOOLS::NLO_subevtlist& subevtlist);
    bool AllowAsSpecInISPFF(const size_t &k);
    bool AllowAsSpecInFSPFF(const size_t &k);

    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Real_Correction();
    ~Single_Real_Correction();
 
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
    bool Combinable(const size_t &idi,const size_t &idj);
    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);
    void SetCaller(PHASIC::Process_Base *const proc);

    void SetFixedScale(const std::vector<double> &s);
    void SetSelectorOn(const bool on);
    void SetGenerator(PHASIC::ME_Generator_Base *const gen);
    size_t SetMCMode(const size_t mcmode);
    size_t SetClusterMode(const size_t cmode);
    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  public:
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return m_newlib;}
    int                 InitAmplitude(Amegic_Model *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    bool                SetUpIntegrator();
    Amplitude_Handler * GetAmplitudeHandler()          { return p_tree_process->GetAmplitudeHandler();}
    Helicity *          GetHelicity()                  { return p_tree_process->GetHelicity(); }    
    double              Result()                       { return m_iresult; } 

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
    void SetNLOMC(PDF::NLOMC_Base *const mc);

    inline ATOOLS::NLO_subevtlist*             GetSubevtList()
    { return &m_subevtlist; }
    inline Single_DipoleTerm*                  GetSubTerm(size_t i)
    { return m_subtermlist[i]; }
    inline size_t                              GetSubTermNumber()
    { return m_subtermlist.size(); }

    int Type() { return 100; }

    inline void SetNoTree(const int notree) { m_no_tree=notree; }    
    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    double Partonic(const ATOOLS::Vec4D_Vector&,
                    ATOOLS::Variations_Mode varmode,
                    int mode);
    double         operator()(const ATOOLS::Vec4D_Vector &,const int mode);
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                        std::vector<std::vector<Complex> >& cols);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);

    void SetScale(const PHASIC::Scale_Setter_Arguments &args);
    void SetKFactor(const PHASIC::KFactor_Setter_Arguments &args);

    int            NumberOfDiagrams();
    Point        * Diagram(int i);

    ATOOLS::Flavour ReMap(const ATOOLS::Flavour &fl,const size_t &id) const;

   /*------------------------------------------------------------------------------

      Helpers
      
      ------------------------------------------------------------------------------*/
  public:
    void           PrintProcessSummary(int=0); 
    void           PrintSubevtSummary();
    size_t         NumberOfDipoles() {return m_subtermlist.size();} 
    Process_Base *GetReal();

  };
}



#endif

