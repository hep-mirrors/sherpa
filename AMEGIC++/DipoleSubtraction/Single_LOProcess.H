#ifndef AMEGIC_DipoleSubtraction_Single_LOProcess_H
#define AMEGIC_DipoleSubtraction_Single_LOProcess_H

#include "AMEGIC++/Main/Process_Base.H"

#include "AMEGIC++/Main/Polarisation.H"
#include "AMEGIC++/Amplitude/Amplitude_Handler.H"
// #include "PHASIC++/Process/Single_Process.H"
#include <vector>


namespace AMEGIC {
  class Helicity;

  typedef std::set<std::pair<size_t,size_t> >     Combination_Set;
  typedef std::map<size_t,ATOOLS::Flavour_Vector> CFlavVector_Map;

  class Single_LOProcess : public Process_Base {
  protected:
    int                     m_gen_str;
    std::string             m_ptypename,m_libname,m_pslibname;
    ATOOLS::sbt::subtype    m_stype;
    Polarisation            m_pol;
    bool                    m_newlib;
    double                  m_iresult;
    int                     m_emit,m_spect,m_type;

    Helicity              * p_hel;
    Basic_Sfuncs          * p_BS;
    Amplitude_Handler     * p_ampl;
    String_Handler        * p_shand;

    Single_LOProcess      * p_partner;
    std::vector<size_t>     m_partonlistqcd, m_partonlistqed;
    std::vector<int>        m_rsmap,m_srmap;

    size_t                  m_pspissplscheme, m_pspfssplscheme;

    std::vector<double>     m_maxcpliqcd, m_mincpliqcd,
                            m_maxcpliew, m_mincpliew;

    std::map<int, size_t>   m_flavrestrictions;

    std::vector<ATOOLS::Vec4D> m_epol;

    ATOOLS::NLO_subevt *p_sub;

    Combination_Set m_ccombs;
    CFlavVector_Map m_cflavs;

    void FillCombinations(Point *const p,size_t &id);
    void FillCombinations();    
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                        std::vector<std::vector<Complex> >& cols,
                        double sfactor);

    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

   Single_LOProcess(const PHASIC::Process_Info &pi,
                    BEAM::Beam_Spectra_Handler *const beam,
                    PDF::ISR_Handler *const isr,
                    YFS::YFS_Handler *const yfs,
                    const ATOOLS::sbt::subtype& st);
    ~Single_LOProcess();

    bool Combinable(const size_t &idi,const size_t &idj);
    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);
    
    const MODEL::Coupling_Map *CouplingMap() const { return &m_cpls; }
    void SetProcMap(PHASIC::NLOTypeStringProcessMap_Map *apmap) { p_apmap=apmap; }
    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  private :
    int                 Tests(std::vector<double>* =NULL);
    int                 TestLib(std::vector<double>* =NULL);
    int                 CheckLibraries(std::vector<double>* =NULL);
    int                 CheckStrings(Single_LOProcess*,std::vector<double>* =NULL);
    void                WriteLibrary();
    bool                FoundLib(std::string& pID);
  protected:
    void                CreateMappingFile(Single_LOProcess*);
    void                UpdateMappingFile(std::string,std::map<std::string,Complex> &); 
    void                WriteAlternativeName(std::string aname);
    bool                CheckAlternatives(std::vector<Process_Base *> & links,
					  std::string procname);
    bool                CompareTestMoms(const ATOOLS::Vec4D*);
  public:
    std::string         CreateLibName();
    void                AddChannels(std::list<std::string>*);
    bool                NewLibs() {return m_newlib;}
    bool                FoundMappingFile(std::string &,std::string &);
    int                 InitAmplitude(Amegic_Model *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &);
    virtual int         InitAmplitude(Amegic_Model *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &,int);
    virtual int         InitAmplitude(Amegic_Model *,Topology *,
				      std::vector<Process_Base *> &,
				      std::vector<Process_Base *> &,
				      std::vector<ATOOLS::Vec4D>*,std::vector<double>*);
    bool                SetUpIntegrator();
    String_Handler    * GetStringHandler();
    Amplitude_Handler * GetAmplitudeHandler();
    Helicity *          GetHelicity();
    double              Result()                       { return m_iresult; } 

    int                             NumberOfDiagrams();
    AMEGIC::Point                 * Diagram(int i);
    bool                IsFreeOfFourVertex(Point * _p) { return 1-p_ampl->ExistFourVertex(_p); }

    Basic_Sfuncs      * GetBS() {return p_BS;} 
    double              Norm()  {return m_Norm;}

    inline size_t       GetEmit()  { return m_emit; }
    inline size_t       GetSpect() { return m_spect; }
    inline void         SetSubEvt(ATOOLS::NLO_subevt *sub) { p_sub=sub; }
    inline void         SetPhotonSplittingModes(const size_t &im,
                                                const size_t &fm)
                        { m_pspissplscheme=im; m_pspfssplscheme=fm; }
    inline void         SetMaxOrdersIQCD(const std::vector<double> &o)
                        { m_maxcpliqcd=o; }
    inline void         SetMinOrdersIQCD(const std::vector<double> &o)
                        { m_mincpliqcd=o; }
    inline void         SetMaxOrdersIEW(const std::vector<double> &o)
                        { m_maxcpliew=o; }
    inline void         SetMinOrdersIEW(const std::vector<double> &o)
                        { m_mincpliew=o; }
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    std::string      LibName()                          { return m_libname;     }
    std::string      PSLibName()                        { return m_pslibname;   }
    Process_Base   * Partner()   const                  { return p_partner;     }
    void             Minimize();
    inline int  Type() { return m_type; }
    inline const ATOOLS::sbt::subtype& GetSubType() const   { return m_stype; }
    inline const std::vector<size_t>& PartonListQCD() const { return m_partonlistqcd; }
    inline const std::vector<size_t>& PartonListQED() const { return m_partonlistqed; }
    inline const bool HasInitialStateQCDEmitter() const
    { return m_partonlistqcd.size()?m_partonlistqcd[0]<m_nin:false; }
    inline const bool HasInitialStateQEDEmitter() const
    { return m_partonlistqed.size()?m_partonlistqed[0]<m_nin:false; }
    inline const bool HasInitialStateEmitter() const
    { return HasInitialStateQCDEmitter() || HasInitialStateQEDEmitter(); }
    inline Polarisation *GetPolarisation()             { return &m_pol; }
    bool CheckIQCDMappability() const;
    bool CheckIQEDMappability() const;
    bool IsValid();

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  private:
    double         Calc_M2ik(const int& ci, const int& ck,
                             const std::vector<double>&,
                             const std::vector<double>&);
  public:
    double Partonic(const ATOOLS::Vec4D_Vector&,
                    ATOOLS::Variations_Mode,
                    int mode);
    virtual double operator()(const ATOOLS::Vec4D_Vector &,const ATOOLS::Vec4D *,
			      std::vector<double>*,std::vector<ATOOLS::Vec4D>*,const int);
    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                        std::vector<std::vector<Complex> >& cols);
    virtual void   Calc_AllXS(const ATOOLS::Vec4D_Vector &,
                              const ATOOLS::Vec4D *,
                              std::vector<std::vector<double> > &,
                              std::vector<std::vector<double> > &,
                              const int);
    inline const std::vector<int> &RSMap() const { return m_rsmap; }
    inline const std::vector<int> &SRMap() const { return m_srmap; }
  };
}



#endif

