#ifndef COMIX_Main_Single_Process_H
#define COMIX_Main_Single_Process_H

#include "COMIX/Main/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"

namespace COMIX {

  class Matrix_Element;

  class Single_Process: 
    public PHASIC::Single_Process, public Process_Base {
  protected:

    Matrix_Element *p_bg;
    Single_Process *p_map;

    ATOOLS::Vec4D_Vector m_p;
    ATOOLS::Flavour_Map  m_fmap;

    double m_etime, m_en;

    void InitPSGenerator(const size_t &ismode);

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
    void SetKFactorOn(const bool on);

    Matrix_Element *GetME() const;
    Amplitude      *GetAmplitude() const;
    
    double Differential(const ATOOLS::Vec4D_Vector &p);
    double Differential(const ATOOLS::Cluster_Amplitude &ampl);
    double Differential2();

    bool GeneratePoint();
    bool Tests();

    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
    void UpdateIntegrator(PHASIC::Phase_Space_Handler *const psh);

    ATOOLS::Flavour ReMap(const ATOOLS::Flavour &fl) const;

  };// end of class Single_Process

}// end of namespace COMIX

#endif



