#include "PI_Interface.H"

#include "Run_Parameter.H"

using namespace PHASIC;

PI_Interface::PI_Interface(Phase_Space_Handler *const pshandler,
			   const std::string &key,const size_t dim):
  p_integrator(new ATOOLS::Primitive_Integrator()), p_pshandler(pshandler), 
  m_key(key), m_initialize(false), m_point(dim), m_integral(0.0),
  m_nopt(dim*200), m_nmax(dim*100000), m_ncells(dim*100)
{
  p_integrator->SetDimension(dim);
  p_integrator->SetFunction(this);
  p_integrator->SetVariableName("xs");
  p_integrator->SetUnitName("pb");
  p_integrator->SetScale(ATOOLS::rpa.Picobarn());
}

PI_Interface::~PI_Interface()
{
}

bool PI_Interface::Initialize()
{
  p_integrator->SetError(p_pshandler->Error());
  p_integrator->SetMode(0);
  p_integrator->SetNOpt(m_nopt);
  p_integrator->SetNMax(m_nmax);
  p_integrator->SetNCells(m_ncells);
  m_initialize=true;
  m_integral=p_integrator->Integrate(this);
  m_initialize=false;
  return true;
}

void PI_Interface::GeneratePoint()
{
  p_integrator->Point(m_point);
}

double PI_Interface::GenerateWeight()
{
  return m_integral*p_integrator->Weight(m_point);
}

double PI_Interface::operator()(const std::vector<double> &x) const
{
  PI_Interface *cur=(PI_Interface*)this;
  cur->m_point=x;
  double value=p_pshandler->
    Differential(p_pshandler->Active(),m_mode|psm::pi_call);
  return value;
}

bool PI_Interface::WriteOut(const std::string &path) const
{
  return p_integrator->WriteOut(path+m_key);
}

bool PI_Interface::ReadIn(const std::string &path)
{
  bool success=p_integrator->ReadIn(path+m_key);
  m_integral=p_integrator->Mean();
  return success;
}
