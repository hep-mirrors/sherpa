#ifndef LesHouches_Interface_h
#define LesHouches_Interface_h

#include "Spectrum_Generator_Base.H"
#include "Data_Read.H"
#include "Data_Reader.H"
#include "Matrix.H"
#include "Exception.H"

namespace MODEL {
  class LesHouches_Interface : public Spectrum_Generator_Base, 
			       public ATOOLS::Terminator_Object {
    void Info();
    void SetSMInput(); 
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
    LesHouches_Interface(ATOOLS::Data_Read *,Model_Base *,std::string);
    ~LesHouches_Interface();
    void Run(std::string);
    void FillMasses();
    void FillDecays() {}
    bool FillDecay(ATOOLS::Decay_Table *) {return true;}
  };
}
#endif
