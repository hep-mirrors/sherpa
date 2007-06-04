#include "XS_Model_Base.H"

#include "Data_Read.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace EXTRAXS;

XS_Model_Base::XS_Model_Base(const std::string &name): 
  m_name(name), p_model(NULL) {}

XS_Model_Base::~XS_Model_Base()
{
}

void XS_Model_Base::Initialize(MODEL::Model_Base *const model,
			       const std::string &file)
{
  msg_Debugging()<<METHOD<<"(\""<<file<<"\"): {\n";
  p_model=model;
  Data_Read read(file);
  m_cplscheme=read.GetValue<std::string>("COUPLING_SCHEME","Running");
  m_yukscheme=read.GetValue<std::string>("YUKAWA_MASSES","Running");
  m_widthscheme=read.GetValue<std::string>("WIDTH_SCHEME","Fixed");
  msg_Debugging()<<"}\n";
}

bool XS_Model_Base::IncludesModel(const std::string &name) const
{
  return false;
}

std::vector<std::string> XS_Model_Base::IncludedModels() const
{
  return std::vector<std::string>();
}

int XS_Model_Base::ScalarNumber(const std::string _name) 
{
  return p_model->ScalarNumber(_name);
}

double XS_Model_Base::ScalarConstant(const std::string _name) 
{
  return p_model->ScalarConstant(_name);
}

ATOOLS::CMatrix XS_Model_Base::ComplexMatrix(const std::string _name) 
{
  return p_model->ComplexMatrix(_name);
}

Complex XS_Model_Base::ComplexMatrixElement
(const std::string _name,const int _i,const int _j) 
{
  return p_model->ComplexMatrixElement(_name,_i,_j);
}

ATOOLS::Function_Base * XS_Model_Base::ScalarFunction(const std::string _name)
{
  return p_model->GetScalarFunction(_name);
}

double XS_Model_Base::ScalarFunction(const std::string _name,double _t) 
{
  if (p_model->GetScalarFunction(_name)->Type()=="Running Coupling") {
    if (m_cplscheme=="Running") return p_model->ScalarFunction(_name,_t);
    if (m_cplscheme=="Running_alpha_S" && _name=="alpha_S")
      return p_model->ScalarFunction(_name,_t);
    if (m_cplscheme=="Running_alpha_QED" && _name=="alpha_QED")
      return p_model->ScalarFunction(_name,_t);
    return p_model->ScalarFunction(_name);
  }
  if (p_model->GetScalarFunction(_name)->Type()=="Running Mass") {
    if (m_yukscheme=="Running") return p_model->ScalarFunction(_name,_t);
    return p_model->ScalarFunction(_name);
  }
  return 0.0;
}
