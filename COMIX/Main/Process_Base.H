#ifndef COMIX_Main_Process_Base_H
#define COMIX_Main_Process_Base_H

#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Smart_Pointer.H"

namespace BEAM { class Beam_Spectra_Handler; }

namespace PDF { class ISR_Handler; }

namespace PHASIC { 
  class Process_Base; 
  class Phase_Space_Handler;
}

namespace COMIX {

  class PS_Generator;
  class Amplitude;
  class Single_Process;

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;

  class Process_Base {
  protected:

    PHASIC::Process_Base *p_proc;
    MODEL::Model_Base    *p_model;
    
    SP(PS_Generator) p_psgen;

    Int_Matrix m_colours;

    int    m_cls, m_hls;

    std::string m_gpath;

    std::map<std::string,std::string> *p_pmap;
    std::vector<Single_Process*>      *p_umprocs;

  public:

    // constructors
    Process_Base(PHASIC::Process_Base *const p_proc,
		 MODEL::Model_Base *const model=NULL);

    // destructor
    virtual ~Process_Base();

    // member functions
    virtual bool Initialize(std::map<std::string,std::string> *const pmap,
			    std::vector<Single_Process*> *const procs);
    virtual bool MapProcess() = 0;

    virtual void InitPSGenerator(const size_t &ismode) = 0;

    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);

    virtual bool Tests() = 0;
    virtual void ConstructPSVertices(PS_Generator *ps) = 0;

    // inline functions
    inline PHASIC::Process_Base *Process() const { return p_proc; }

    inline Int_Matrix Colours() const { return m_colours; }

    inline MODEL::Model_Base *GetModel() const { return p_model; }

    inline void SetGPath(const std::string &gpath) { m_gpath=gpath; }

    inline SP(PS_Generator) PSGenerator() const { return p_psgen; }

    inline void SetModel(MODEL::Model_Base *const model) { p_model=model; }

  };// end of class Process_Base

  std::string ComixLogo();

}// end of namespace COMIX

#endif



