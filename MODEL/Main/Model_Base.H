#ifndef MODEL_Main_Model_Base_H
#define MODEL_Main_Model_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Main/Coupling_Data.H"
#include <map>

namespace ATOOLS { 
  class Decay_Table; 
  class Data_Reader;
} 

namespace MODEL {

  class Spectrum_Generator_Base;
  class Vertex;
  class Single_Vertex;

  typedef  std::vector<Single_Vertex *>  Vertex_List;
  typedef  std::map<ATOOLS::Flavour, Vertex_List> Vertex_Table;

  typedef std::map<std::string,int>                     ScalarNumbersMap; 
  typedef std::map<std::string,double>                  ScalarConstantsMap; 
  typedef std::map<std::string,Complex>                 ComplexConstantsMap; 
  typedef std::map<std::string,ATOOLS::Function_Base *> ScalarFunctionsMap; 
  typedef std::map<std::string,ATOOLS::CMatrix>         ComplexMatricesMap; 

  class Interaction_Model_Base;

  struct Model_Arguments {
    std::string m_path, m_file;
    bool m_elementary;
    Model_Arguments(const std::string &path,const std::string &file,const bool &elementary):
      m_path(path), m_file(file), m_elementary(elementary) {}
  };// end of struct Model_Arguments

  class Model_Base {
  public:

    typedef ATOOLS::Getter_Function<Model_Base,Model_Arguments> 
    Model_Getter_Function;

  protected:
    Interaction_Model_Base  * p_model;

    std::string               m_dir,m_file;
    std::string               m_name;
    bool                      m_elementary;

    ATOOLS::Data_Reader     * p_dataread;
    ScalarNumbersMap        * p_numbers;
    ScalarConstantsMap      * p_constants;
    ComplexConstantsMap     * p_complexconstants;
    ScalarFunctionsMap      * p_functions;
    ComplexMatricesMap      * p_matrices; 

    Spectrum_Generator_Base * p_spectrumgenerator;
    Vertex                  * p_vertex;
    Vertex_Table            * p_vertextable;

    void CustomContainerInit();

    virtual void ParticleInit() = 0;
    virtual void FillSpectrum() = 0;
    void ReadParticleData();
  public:
    Model_Base(std::string,std::string,bool);
    virtual ~Model_Base();

    virtual bool ModelInit() = 0;

    static void ShowSyntax(const size_t mode);
    void InitializeInteractionModel();
    void GetCouplings(Coupling_Map &cpls) const;

    std::string           Name() const          { return m_name; } 
    ScalarNumbersMap    * GetScalarNumbers()    { return p_numbers; }
    ScalarFunctionsMap  * GetScalarFunctions()  { return p_functions; }
    ScalarConstantsMap  * GetScalarConstants()  { return p_constants; }
    ComplexConstantsMap * GetComplexConstants() { return p_complexconstants; }
    ComplexMatricesMap  * GetComplexMatrices()  { return p_matrices; }

    inline Interaction_Model_Base *GetInteractionModel() { return p_model; }

    inline ScalarNumbersMap * ExtractScalarNumbers() {
      ScalarNumbersMap * temp(p_numbers);
      p_numbers=NULL;
      return temp;
    }
    inline ScalarFunctionsMap * ExtractScalarFunctions() {
      ScalarFunctionsMap * temp(p_functions);
      p_functions=NULL;
      return temp;
    }
    inline ScalarConstantsMap * ExtractScalarConstants() {
      ScalarConstantsMap * temp(p_constants);
      p_constants=NULL;
      return temp;
    }
    inline ComplexConstantsMap * ExtractComplexConstants() {
      ComplexConstantsMap * temp(p_complexconstants);
      p_complexconstants=NULL;
      return temp;
    }
    inline ComplexMatricesMap * ExtractComplexMatrices() {
      ComplexMatricesMap * temp(p_matrices);
      p_matrices=NULL;
      return temp;
    }

    virtual int                     ScalarNumber(const std::string);
    virtual double                  ScalarConstant(const std::string);
    virtual Complex                 ComplexConstant(const std::string);
    virtual ATOOLS::Function_Base * GetScalarFunction(const std::string);
    virtual double                  ScalarFunction(const std::string,double);
    virtual double                  ScalarFunction(const std::string);
    virtual ATOOLS::CMatrix         ComplexMatrix(const std::string);
    virtual Complex                 ComplexMatrixElement(const std::string,const int,const int);

    inline  Vertex       * GetVertex()            { return p_vertex; }
    inline  Vertex_Table * GetVertexTable() const { return p_vertextable; }

    virtual bool CheckFlavours(int nin, int nout, ATOOLS::Flavour* flavs);
  };

  extern Model_Base *s_model;

  /*!
    \namespace MODEL
    \brief The namespace to host and evaluate parameters of physical models
  */
  /*!
    \class Model_Base
    \brief is the base class for physical models.

    This class is the base class for all physical models within the SHERPA framework.
    It contains (and evaluates where needed) all physical parameters like couplings and masses
    of the particles. In some cases, also branching widths are evaluated to some higher order here
    and might be used to fix specific couplings in the vertices of the AMEGIC::Interaction_Model_Base.
    These calculations are performed by external spectrum generators that can be wrapped up through
    the base class MODEL::Spectrum_Generator_Base, a pointer to such an object is part of
    every model. Furthermore, all relevant parameters, including running couplings, are hosted here 
    and are accessible through searchable maps on the basis of word-strings.
  */
  /*!
    \var std::string MODEL_Base::m_dir
    The directory for the data file to be read-in.
  */
  /*!
    \var std::string MODEL_Base::m_file
    The name of the data file.
  */
  /*!
    \var std::string Model_Base::m_name
    The name of the model, e.g. SM, MSSM, ADD, ...
  */
  /*!
    \var ScalarNumbersMap * Model_Base::p_numbers
    A pointer to a searchable map of strings connected with integers defining the model.
    This - and all other - searchable maps are the main ingredients of a particular model
    to be communicated to the outside world.
  */
  /*!
    \var ScalarConstantsMap * Model_Base::p_constants
    A pointer to a searchable map of strings connected with double precision numbers defining 
    the model. These numbers might be, for instance, masses of the particles, fixed couplnigs, etc..
    This - and all other - searchable maps are the main ingredients of a particular model to be communicated 
    to the outside world. 
  */
  /*!
    \var ScalarFunctionsMap * Model_Base::p_functions
    A pointer to a searchable map of strings connected with function base objects defining 
    the model. The function bases might be objects like, e.g. running couplings and the like.
    This - and all other - searchable maps are the main ingredients of a particular model to be 
    communicated to the outside world.
  */
  /*!
    \var ComplexMatricesMap * Model_Base::p_matrices 
    A pointer to a searchable map of strings connected with complex matrices defining the model. Examples
    for such matrices include the CKM matrix, SUSY mixing matrices, etc. . This - and all other - searchable 
    maps are the main ingredients of a particular model to be communicated to the outside world. Possibly,
    more maps will have to be added in the future, for instance for the complex mass scheme this might
    be an option.
  */
  /*!
    \var ATOOLS::Data_Read * Model_Base::p_dataread
    A pointer to the Data_Read object used to read in all parameters needed to define the model.
  */
  /*!
    \var int Model_Base::m_spectrum
    The mode of the spectrum generator that might be used for the caluclation of particle masses and branching 
    ratios.
  */
  /*!
    \var Spectrum_Generator_Base * Model_Base::p_spectrumgenerator
    A pointer to the spectrum generator object to be used in the calculation of masses, branching widths
    and the like. Right now, only one such generator is actually implemented, namely LesHouches.
  */
  /*!
    \fn Model_Base::Model_Base(std::string,std::string)
    The constructor of the abstract Model_Base just initalizes all pointers (NULL) and the directory path 
    and filename of the input data file.
  */
  /*!
    \fn Model_Base::~Model_Base()
    In the destructor all maps are emptied, all members deleted.
  */
  /*!
    \fn virtual void Model_Base::FillSpectrum()
    A virtual method to read in the data file and set the parameters, this will be specified on a model-by-model base. 
  */
  /*!
    \fn virtual int Model_Base::ScalarNumber(const std::string)
    Through this method the integer constants can be accessed from outside. Within the map,
    it is looked for the string, if found, the corresponding number is returned. Otherwise,
    after an error message, zero is returned.
  */
  /*!
    \fn virtual double Model_Base::ScalarConstant(const std::string)
    Through this method the double precision constants can be accessed from outside. Within the map,
    it is looked for the string, if found, the corresponding number is returned. Otherwise,
    after an error message, zero is returned.
  */
  /*!
    \fn virtual ATOOLS::Function_Base * Model_Base::GetScalarFunction(const std::string)
    Through this method the pointers to the function base objects can be accessed from outside. Within 
    the map, it is looked for the string, if found, the corresponding function is returned. Otherwise,
    after an error message, a NULL-pointer is returned.
  */
  /*!
    \fn virtual double Model_Base::ScalarFunction(const std::string,double)
    Through this method the values of the function base objects for a certain value can be accessed from 
    outside. Within the map, it is looked for the string, if found, the corresponding value is evaluated
    within the function base and returned. Otherwise, after an error message, zero is returned.
  */
  /*!
    \fn virtual double Model_Base::ScalarFunction(const std::string)
    Through this method the default values of the function base objects can be accessed from 
    outside. Within the map, it is looked for the string, if found, the corresponding value is evaluated
    within the function base and returned. Otherwise, after an error message, zero is returned.
  */
  /*!
    \fn virtual ATOOLS::CMatrix Model_Base::ComplexMatrix(const std::string)
    Through this method the complex matrices can be accessed from outside. Within the map, it is looked 
    for the string, if found, the corresponding matrix is returned. Otherwise, after an error message, an 
    idendity matrix of rank 1 is returned.
  */
  /*!
    \fn virtual Complex Model_Base::ComplexMatrixElement(const std::string,const int,const int)
    Through this method the components of the complex matrices can be accessed from 
    outside. Within the map, it is looked for the string, if found, the corresponding component is returned.
    Otherwise, i.e. if either the matrix does not exist or the indices are out of range, after an error message, 
    zero is returned.
  */
}

#endif


