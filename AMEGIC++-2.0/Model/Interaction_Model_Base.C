#include "Interaction_Model_Base.H"
#include "Message.H"

using namespace AMEGIC;

Interaction_Model_Base::Interaction_Model_Base(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  p_model(_model), p_vertex(NULL), m_cplscheme(_cplscheme), m_yukscheme(_yukscheme)  { }

void Interaction_Model_Base::Init_Vertex() {
  AORGTOOLS::msg.Debugging()<<"Initialize new vertices !"<<std::endl;
  if (p_vertex!=NULL) delete p_vertex;
  p_vertex = new Vertex(this);
} 

int Interaction_Model_Base::ScalarNumber(const std::string _name) {
  return p_model->ScalarNumber(_name);
}

double Interaction_Model_Base::ScalarConstant(const std::string _name) {
  return p_model->ScalarConstant(_name);
}

AMATOOLS::CMatrix Interaction_Model_Base::ComplexMatrix(const std::string _name) {
  return p_model->ComplexMatrix(_name);
}

Complex Interaction_Model_Base::ComplexMatrixElement(const std::string _name,const int _i,const int _j) {
  return p_model->ComplexMatrixElement(_name,_i,_j);
}

AMATOOLS::Function_Base * Interaction_Model_Base::ScalarFunction(const std::string _name) {
  return p_model->GetScalarFunction(_name);
}

double Interaction_Model_Base::ScalarFunction(const std::string _name,double _t) {
  if (p_model->GetScalarFunction(_name)->Type()==string("Running Coupling")) {
    AORGTOOLS::msg.Out()<<"Match for Running Coupling : "<<_name<<endl;
    if (m_cplscheme==std::string("Running")) return p_model->ScalarFunction(_name,_t);
    return p_model->ScalarFunction(_name);
  }
  if (p_model->GetScalarFunction(_name)->Type()==string("Running Mass")) {
    AORGTOOLS::msg.Out()<<"Match for Running Mass : "<<_name<<endl;
    if (m_yukscheme==std::string("Running")) return p_model->ScalarFunction(_name,_t);
    return p_model->ScalarFunction(_name);
  }
  AORGTOOLS::msg.Out()<<"Miss : "<<_name<<endl;
}


