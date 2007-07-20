#ifndef XS_Model_Base_H
#define XS_Model_Base_H

#include "Model_Base.H"
#include "MyComplex.H"

#define CZERO Complex(0.0,0.0)

namespace EXTRAXS {

  typedef std::map<std::string,Complex> Complex_Map;

  class XS_Model_Base {
  private:

    std::string m_name;

  protected:

    MODEL::Model_Base *p_model;

    std::string m_cplscheme, m_masscheme, m_widthscheme, m_yukscheme;

    Complex_Map m_consts;

  public:

    // constructor
    XS_Model_Base(const std::string &name);

    virtual ~XS_Model_Base();

    // member functions
    virtual void Initialize(MODEL::Model_Base *const model,
			    const std::string &file);

    int    ScalarNumber(const std::string _name);
    double ScalarConstant(const std::string _name);

    ATOOLS::CMatrix        ComplexMatrix(const std::string _name);
    ATOOLS::Function_Base *ScalarFunction(const std::string name);

    double  ScalarFunction(const std::string name,double t);
    Complex ComplexMatrixElement(const std::string name,
				 const int i,const int j);

    virtual bool IncludesModel(const std::string &name) const;

    virtual std::vector<std::string> IncludedModels() const;

    // inline functions
    inline std::string Name() const { return m_name; }

    inline Complex Constant(const std::string &name) 
    { return m_consts.find(name)!=m_consts.end()?m_consts[name]:CZERO; }
    
  };// end of class XS_Model_Base

}// end of namespace EXTRAXS

#endif
