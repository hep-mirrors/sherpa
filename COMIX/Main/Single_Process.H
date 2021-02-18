#ifndef COMIX_Main_Single_Process_H
#define COMIX_Main_Single_Process_H

#include "ATOOLS/Phys/NLO_Types.H"
#include "COMIX/Main/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/KP_Terms.H"
#include "COMIX/Amplitude/Amplitude.H"

namespace COMIX {

  class Amplitude;

  typedef std::set<PS_Generator*> PSGenerator_Set;

  class Single_Process: 
    public PHASIC::Single_Process, public Process_Base {

  protected:

    Amplitude      *p_bg;
    Hard_Matrix    *p_hc;
    Single_Process *p_map;

    PHASIC::Virtual_ME2_Base *p_loop;
    PHASIC::KP_Terms         *p_kpterms;

    ATOOLS::Vec4D_Vector m_p;
    ATOOLS::Flavour_Map  m_fmap;

    Combination_Set m_ccombs;
    CFlavVector_Map m_cflavs;
    PSGenerator_Set m_psset;

    double m_w, m_dxs, m_x[2];
    bool   m_checkpoles, m_allowmap, m_hasmapfile;

    ATOOLS::cs_itype::type m_itype;
    ATOOLS::NLO_subevtlist m_subs;

    void InitPSGenerator(const size_t &ismode);

    void MapSubEvts(const int mode);

    double SetZero();

    void FillMEWeights(ATOOLS::ME_Weight_Info &wgtinfo) const;

    void UpdateKPTerms(const int mode);
    double KPTerms(const int mode,PDF::PDF_Base *pdfa,
		   PDF::PDF_Base *pdfb,double scalefac2=1.0);

  public:

    // constructors
    Single_Process();

    // destructor
    ~Single_Process();

    // member functions
    bool Initialize(std::map<std::string,std::string> *const pmap,
		    std::vector<Single_Process*> *const procs,
		    const std::vector<int> &blocks,size_t &nproc);
    bool MapProcess();

    int PerformTests();

    void ConstructPSVertices(PS_Generator *ps);

    Amplitude *GetAmplitude() const;

    double Partonic(const ATOOLS::Vec4D_Vector &p,int mode);
    ATOOLS::Weights_Map Differential(const ATOOLS::Cluster_Amplitude& ampl,
                                     ATOOLS::Variations_Mode type=ATOOLS::Variations_Mode::nominal_only,
                                     int mode=0);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);

    bool GeneratePoint();
    bool Tests();

    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);

    bool Combinable(const size_t &idi,const size_t &idj);

    const ATOOLS::Flavour_Vector &CombinedFlavour(const size_t &idij);

    void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes> &amps,
			std::vector<std::vector<Complex> > &cols);

    ATOOLS::Flavour ReMap(const ATOOLS::Flavour &fl,
			  const size_t &id) const;

    ATOOLS::NLO_subevtlist *GetSubevtList();
    ATOOLS::NLO_subevtlist *GetRSSubevtList();

    void SetScale(const PHASIC::Scale_Setter_Arguments &args);
    void SetShower(PDF::Shower_Base *const ps);
    void SetNLOMC(PDF::NLOMC_Base *const mc);

    size_t SetMCMode(const size_t mcmode);
    void   SetLookUp(const bool lookup);

  };// end of class Single_Process

}// end of namespace COMIX

#endif



