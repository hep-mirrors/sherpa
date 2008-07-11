#ifndef Analysis_Handler_H
#define Analysis_Handler_H

#include "Primitive_Analysis.H"
#include "Primitive_Observable_Base.H"
#include "File_IO_Base.H"
#include "Exception.H"

namespace ANALYSIS {

  typedef ANALYSIS::Argument_Matrix Argument_Matrix;
  
  class Analysis_Handler: public ATOOLS::Terminator_Object,
			  public ATOOLS::File_IO_Base {
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
    bool   m_initialized;

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
    void WriteOut(const std::string &path="",const bool cont=true);
    void Finish(const std::string &path="");
    void Clear();
    
    static void ShowSyntax(const size_t i);

    void Test(const int mode=0);

    static Argument_Matrix FindArguments(const Argument_Matrix &strings,
					 size_t &starty,size_t &startx);

    // inline functions
    inline void SetWeighted(const int weighted) { m_weighted=weighted; }
    
    inline int Weighted() const { return m_weighted; }

  };// end of class Analysis_Handler   

}// end of namespace ANALYSIS

#endif
