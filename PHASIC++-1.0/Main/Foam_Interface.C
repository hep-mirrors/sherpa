#include "Foam_Interface.H"

#include "Phase_Space_Handler.H"

using namespace PHASIC;

Foam_Interface::Foam_Interface(Phase_Space_Handler *const pshandler,
			       const std::string &key,const size_t dim):
  p_ran(NULL), p_foam(NULL), p_pshandler(pshandler), 
  m_key(key), m_initialize(false), m_values(dim), m_integral(0.0)
{
}

Foam_Interface::~Foam_Interface()
{
}

bool Foam_Interface::Initialize()
{
  if (p_foam!=NULL) return false;
  p_foam = new TFOAM(m_key.c_str());
  p_ran = new TPSEMAR();
  p_foam->SetkDim(m_values.size());
  p_foam->SetnCells(2000);
  m_initialize=true;
  p_foam->Initialize(p_ran,this);
  double relerr;
  p_foam->Finalize(m_integral,relerr);
  m_initialize=false;
  return true;
}

double Foam_Interface::MCEvent()
{
  return m_integral*p_foam->MCgenerate(&m_values.front());
}

double Foam_Interface::Density(int ndim,double *val)
{
  SetValues(val);
  return p_pshandler->Differential(p_pshandler->Active(),-1);
}

void Foam_Interface::SetValues(double *val)
{
  for (size_t i=0;i<m_values.size();++i) m_values[i]=val[i];
}

