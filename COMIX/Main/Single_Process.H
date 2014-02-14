#ifndef COMIX_Main_Single_Process_H
#define COMIX_Main_Single_Process_H

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
    Single_Process *p_map;

    PHASIC::Virtual_ME2_Base *p_loop;
    PHASIC::KP_Terms         *p_kpterms;

    ATOOLS::Vec4D_Vector m_p;
    ATOOLS::Flavour_Map  m_fmap;

    CFlavVector_Map m_cfmap;
    PSGenerator_Set m_psset;

    double m_w, m_dxs, m_x[2];
    bool   m_checkpoles;

    ATOOLS::NLO_subevtlist m_subs;

    void InitPSGenerator(const size_t &ismode);

    void MapSubEvts(const int mode);

    double SetZero();

    double GetKPTerms(const ATOOLS::Flavour_Vector &fl,const int mode);

    void FillMEWeights(ATOOLS::ME_wgtinfo &wgtinfo) const;

  public:

    // constructors
    Single_Process();

    // destructor
    ~Single_Process();

    // member functions
    bool Initialize(std::map<std::string,std::string> *const pmap,
		    std::vector<Single_Process*> *const procs);
    bool MapProcess();

    void ConstructPSVertices(PS_Generator *ps);

    Amplitude *GetAmplitude() const;
    
    double Partonic(const ATOOLS::Vec4D_Vector &p,const int mode);
    double Differential(const ATOOLS::Cluster_Amplitude &ampl,
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

    void SetScale(const PHASIC::Scale_Setter_Arguments &args);
    void SetShower(PDF::Shower_Base *const ps);

    size_t SetMCMode(const size_t mcmode);
    void   SetLookUp(const bool lookup);

  };// end of class Single_Process

}// end of namespace COMIX

#endif



