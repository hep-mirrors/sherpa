#ifndef Model_Base_h
#define Model_Base_h

#include "Data_Read.H"
#include "Flavour.H"
#include "Matrix.H"
#include "Function_Base.H"
#include "MathTools.H"
#include "MyComplex.H"
#include <map>

namespace MODEL { class Spectrum_Generator_Base; }

namespace MODEL {
  typedef std::map<std::string,int>                       ScalarNumbersMap; 
  typedef std::map<std::string,double>                    ScalarConstantsMap; 
  typedef std::map<std::string,AMATOOLS::Function_Base *> ScalarFunctionsMap; 
  typedef std::map<std::string,AMATOOLS::CMatrix>         ComplexMatricesMap; 


  class Model_Base {
  protected:
    std::string               m_dir,m_file;
    std::string               m_name;

    ScalarNumbersMap        * p_numbers;
    ScalarConstantsMap      * p_constants;
    ScalarFunctionsMap      * p_functions;
    ComplexMatricesMap      * p_matrices; 

    AORGTOOLS::Data_Read    * p_dataread;

    int                       m_spectrum;
    Spectrum_Generator_Base * p_spectrumgenerator;

    virtual void ReadInFile() = 0;
  public:
    Model_Base(std::string,std::string);
    ~Model_Base();
    virtual bool RunSpectrumGenerator() = 0;
    virtual void FillMasses()           = 0;
    virtual void FillWidths()           = 0;

    virtual std::string          Name()               { return m_name; } 
    virtual ScalarNumbersMap   * GetScalarNumbers()   { return p_numbers; }
    virtual ScalarFunctionsMap * GetScalarFunctions() { return p_functions; }
    virtual ScalarConstantsMap * GetScalarConstants() { return p_constants; }
    virtual ComplexMatricesMap * GetComplexMatrices() { return p_matrices; }

    virtual int                       ScalarNumber(const std::string);
    virtual double                    ScalarConstant(const std::string);
    virtual AMATOOLS::Function_Base * GetScalarFunction(const std::string);
    virtual double                    ScalarFunction(const std::string,double);
    virtual double                    ScalarFunction(const std::string);
    virtual AMATOOLS::CMatrix         ComplexMatrix(const std::string);
    virtual Complex                   ComplexMatrixElement(const std::string,const int,const int);
 };
}

#endif


