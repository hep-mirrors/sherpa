#include "MODEL/Interaction_Models/Interaction_Model_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE MODEL::Interaction_Model_Base
#define PARAMETER_TYPE MODEL::Interaction_Model_Arguments
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/Message.H"

using namespace MODEL;

Interaction_Model_Base::Interaction_Model_Base(const std::string &code,MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  p_model(_model), m_code(code), m_cplscheme(_cplscheme), m_yukscheme(_yukscheme), 
  m_tensors(false), m_loops(false), m_agcs(false)
{ 
}

int Interaction_Model_Base::ScalarNumber(const std::string _name) {
  return p_model->ScalarNumber(_name);
}

double Interaction_Model_Base::ScalarConstant(const std::string _name) {
  return p_model->ScalarConstant(_name);
}

Complex Interaction_Model_Base::ComplexConstant(const std::string _name) {
  return p_model->ComplexConstant(_name);
}

ATOOLS::CMatrix Interaction_Model_Base::ComplexMatrix(const std::string _name) {
  return p_model->ComplexMatrix(_name);
}

Complex Interaction_Model_Base::ComplexMatrixElement(const std::string _name,const int _i,const int _j) {
  return p_model->ComplexMatrixElement(_name,_i,_j);
}

ATOOLS::Function_Base * Interaction_Model_Base::ScalarFunction(const std::string _name) {
  return p_model->GetScalarFunction(_name);
}

double Interaction_Model_Base::ScalarFunction(const std::string _name,double _t) {
  if (p_model->GetScalarFunction(_name)->Type()==std::string("Running Coupling")) {
    if (m_cplscheme==std::string("Running")) return p_model->ScalarFunction(_name,_t);
    if (m_cplscheme==std::string("Running_alpha_S")&&_name==std::string("alpha_S"))
	return p_model->ScalarFunction(_name,_t);
    if (m_cplscheme==std::string("Running_alpha_QED")&&_name==std::string("alpha_QED"))
	return p_model->ScalarFunction(_name,_t);
    return p_model->ScalarFunction(_name);
  }
  if (p_model->GetScalarFunction(_name)->Type()==std::string("Running Mass")) {
    if (m_yukscheme==std::string("Running")) return p_model->ScalarFunction(_name,_t);
    return p_model->ScalarFunction(_name);
  }
  return 0.0;
}

Interaction_Model_Base::~Interaction_Model_Base()
{
} 

std::string Interaction_Model_Base::Name()    
{ 
  return p_model->Name(); 
} 
