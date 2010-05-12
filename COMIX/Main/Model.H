#ifndef COMIX_Main_Model_H
#define COMIX_Main_Model_H

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Getter_Function.H"

#define CZERO Complex(0.0,0.0)

namespace COMIX {

  typedef std::map<std::string,Complex> Complex_Map;

  class Model {
  private:

    std::string m_name;

  protected:

    MODEL::Model_Base *p_model;

    std::string m_cplscheme, m_masscheme, m_widthscheme, m_yukscheme;

    Complex_Map m_consts;

    static void AddModel(std::vector<std::string> &mds,
			 const std::string &md);
    static void AddFlavour(std::vector<ATOOLS::Flavour> &fls,
			   const ATOOLS::Flavour &fl);

  public:

    // constructor
    Model(const std::string &name);

    virtual ~Model();

    // member functions
    virtual void Initialize(MODEL::Model_Base *const model,
			    const std::string &file);

    int    ScalarNumber(const std::string &name);
    double ScalarConstant(const std::string &name);

    ATOOLS::CMatrix        ComplexMatrix(const std::string &name);
    ATOOLS::Function_Base *ScalarFunction(const std::string &name);

    double  ScalarFunction(const std::string &name,const double &t);
    Complex ComplexMatrixElement(const std::string &name,
				 const int i,const int j);

    bool IncludesModel(const std::string &name) const;
    bool IncludesFlavour(const ATOOLS::Flavour &fl) const;

    virtual std::vector<std::string>     IncludedModels() const;
    virtual std::vector<ATOOLS::Flavour> IncludedFlavours() const;

    // inline functions
    inline std::string Name() const { return m_name; }

    inline MODEL::Model_Base *GetModel() const { return p_model; }

    inline Complex Constant(const std::string &name) 
    { return m_consts.find(name)!=m_consts.end()?m_consts[name]:CZERO; }
    
  };// end of class Model

  typedef ATOOLS::Getter_Function<Model,std::string,
				  std::less<std::string> > Model_Getter;

}// end of namespace COMIX

#endif
