#ifndef Analysis_Handler_H
#define Analysis_Handler_H

#include "Primitive_Analysis.H"
#include "Primitive_Observable_Base.H"
#include "File_IO_Base.H"
#include "Exception.H"

namespace ANALYSIS {
  
  typedef ANALYSIS::String_Matrix String_Matrix;

  class Analysis_Handler: public ATOOLS::Terminator_Object, 
			  public ATOOLS::File_IO_Base {
  public:

    typedef ANALYSIS::Primitive_Observable_Base::Getter_Function 
    Getter_Function;
    typedef std::vector<ANALYSIS::Primitive_Analysis*> Analyses_Vector;
    
  private:

    static size_t s_maxanalyses;

    Analyses_Vector m_analyses;

    String_Matrix FindArguments(const String_Matrix &strings,
				const size_t start) const;

    void Clean();
    void PrepareTerminate();

  public:

    // constructor
    Analysis_Handler();

    // destructor
    ~Analysis_Handler();

    // member functions
    bool ReadIn();

    void DoAnalysis(const ATOOLS::Blob_List *bloblist,const double weight);
    void Finish(const std::string &path="");
    void Clear();
    
    static void ShowSyntax(const size_t i);

  };// end of class Analysis_Handler   

}// end of namespace ANALYSIS

#endif
