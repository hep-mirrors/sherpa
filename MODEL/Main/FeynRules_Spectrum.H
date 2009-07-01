#ifndef MODEL_Main_FeynRules_Spectrum_h
#define MODEL_Main_FeynRules_Spectrum_h

#include "MODEL/Main/Spectrum_Generator_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace MODEL {
  class FeynRules_Spectrum : public Spectrum_Generator_Base, 
			       public ATOOLS::Terminator_Object {
    
    void SetUpTools();
    void FillAlgebra();
    void PrepareReader(std::string);
    void SetExternalParameters();
    void SetInternalParameters();

    std::string                  m_dir, m_identfile,m_paramdeffile,m_paramfile;
    ATOOLS::Data_Reader         *p_reader;
    ATOOLS::Algebra_Interpreter *p_algebra;

    void PrepareTerminate();

  public:
    FeynRules_Spectrum(ATOOLS::Data_Reader *,Model_Base *,std::string);
    ~FeynRules_Spectrum();
  
    void Run();
  };
}
#endif
