#include "COMIX/Main/Model.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE std::string
#define OBJECT_TYPE COMIX::Model
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Smart_Pointer.C"
#include "ATOOLS/Org/Data_Reader.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace COMIX;

namespace ATOOLS { template class SP(Model); }

Model::Model(const std::string &name): 
  m_name(name), p_model(NULL) {}

Model::~Model()
{
}

void Model::Initialize(MODEL::Model_Base *const model,
			    const std::string &file)
{
  msg_Debugging()<<METHOD<<"(\""<<file<<"\"): {\n";
  p_model=model;
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputFile(file);
  m_widthscheme=read.GetValue<std::string>("WIDTH_SCHEME","Fixed");
  msg_Debugging()<<"}\n";
}

void Model::AddModel(std::vector<std::string> &mds,
			       const std::string &md)
{
  for (std::vector<std::string>::const_iterator mit(mds.begin());
       mit!=mds.end();++mit) if (*mit==md) return;
  mds.push_back(md);
}

void Model::AddFlavour(std::vector<ATOOLS::Flavour> &fls,
				 const ATOOLS::Flavour &fl) 
{
  if (!fl.IsOn()) return;
  for (std::vector<Flavour>::const_iterator fit(fls.begin());
       fit!=fls.end();++fit) if (*fit==fl) return;
  fls.push_back(fl);
}

bool Model::IncludesModel(const std::string &name) const
{
  std::vector<std::string> imd(IncludedModels());
  for (size_t i(0);i<imd.size();++i)
    if (imd[i]==name) return true;
  return false;
}

bool Model::IncludesFlavour(const ATOOLS::Flavour &fl) const
{
  std::vector<Flavour> ifl(IncludedFlavours());
  for (size_t i(0);i<ifl.size();++i)
    if (ifl[i]==fl) return true;
  return false;
}

std::vector<std::string> Model::IncludedModels() const
{
  return std::vector<std::string>();
}

std::vector<ATOOLS::Flavour> Model::IncludedFlavours() const
{
  return std::vector<Flavour>();
}

int Model::ScalarNumber(const std::string &name) 
{
  return p_model->GetInteractionModel()->ScalarNumber(name);
}

double Model::ScalarConstant(const std::string &name) 
{
  return p_model->ScalarConstant(name);
}

ATOOLS::CMatrix Model::ComplexMatrix(const std::string &name) 
{
  return p_model->GetInteractionModel()->ComplexMatrix(name);
}

Complex Model::ComplexMatrixElement
(const std::string &name,const int i,const int j) 
{
  return p_model->GetInteractionModel()->ComplexMatrixElement(name,i,j);
}

ATOOLS::Function_Base *Model::ScalarFunction
(const std::string &name)
{
  return p_model->GetInteractionModel()->ScalarFunction(name);
}

double Model::ScalarFunction(const std::string &name,const double &t) 
{
  return p_model->GetInteractionModel()->ScalarFunction(name,t);
}
