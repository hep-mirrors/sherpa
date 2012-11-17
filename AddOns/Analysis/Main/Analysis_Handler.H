#ifndef Analysis_Main_Analysis_Handler_H
#define Analysis_Main_Analysis_Handler_H

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/File_IO_Base.H"
#include "ATOOLS/Org/Exception.H"

namespace ANALYSIS {

  typedef ANALYSIS::Argument_Matrix Argument_Matrix;
  
  class Analysis_Handler: public ATOOLS::Tester_Object,
			  public ATOOLS::File_IO_Base,
			  public SHERPA::Analysis_Interface {
  public:

    typedef ANALYSIS::Primitive_Observable_Base::Getter_Function 
    Observable_Getter_Function;
    typedef ANALYSIS::Analysis_Object::Getter_Function 
    Object_Getter_Function;
    typedef std::vector<ANALYSIS::Primitive_Analysis*> Analyses_Vector;
    
  private:

    static size_t s_maxanalyses;

    Analyses_Vector m_analyses;

    int    m_weighted; 
    bool   m_write;

    void Clean();
    bool ApproveTerminate();

    void DoAnalysis(const ATOOLS::Blob_List *bloblist,const double weight);

  public:

    // constructor
    Analysis_Handler();

    // destructor
    ~Analysis_Handler();

    // member functions
    bool Init();

    bool Run(ATOOLS::Blob_List *const bl);
    bool WriteOut();
    bool Finish();
    void CleanUp();
    
    void ShowSyntax(const int mode);

    static Argument_Matrix FindArguments(const Argument_Matrix &strings,
					 size_t &starty,size_t &startx);

    // inline functions
    inline void SetWeighted(const int weighted) { m_weighted=weighted; }
    
    inline int Weighted() const { return m_weighted; }

    inline SHERPA::Event_Handler *EventHandler() { return p_eventhandler; }

  };// end of class Analysis_Handler   

}// end of namespace ANALYSIS

#endif
