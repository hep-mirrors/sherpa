#ifndef MODEL_Main_LesHouches_Interface_H
#define MODEL_Main_LesHouches_Interface_H

#include "MODEL/Main/Spectrum_Generator_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Org/Exception.H"

namespace MODEL {
  class LesHouches_Interface : public Spectrum_Generator_Base, 
			       public ATOOLS::Terminator_Object {
    void Info();
    void SetSMInput(PDF::ISR_Handler *const isr); 
    void SetMasses(); 
    void SetWidths(); 
    void SetHiggsParameters();
    void SetInoParameters();    
    void SetSquarkParameters();    
    void SetSleptonParameters();    

    std::string          m_dir, m_inputfile;
    ATOOLS::Data_Reader *p_reader;

    void PrepareTerminate();

  public:
    LesHouches_Interface(ATOOLS::Data_Reader *,Model_Base *,std::string);
    ~LesHouches_Interface();
    void Run(PDF::ISR_Handler *const isr);
    void FillMasses();
  };
}
#endif
