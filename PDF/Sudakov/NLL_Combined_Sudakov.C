#include "NLL_Combined_Sudakov.H"
#include "Run_Parameter.H"

using namespace SHERPA;
using namespace ATOOLS;

NLL_Combined_Sudakov::NLL_Combined_Sudakov(int mode) 
{
  m_calcmode = (Sudakov::code)(mode&896);
  m_cutmode  = (Sudakov::code)(mode&7);
}

NLL_Combined_Sudakov::~NLL_Combined_Sudakov()
{
  std::set<NLL_Sudakov_Base*> deleted;
  for (size_t i=0;i<m_suds.size();++i) {
    if (deleted.find(m_suds[i])==deleted.end()) {
      delete m_suds[i]; 
      deleted.insert(m_suds[i]);
    }
  }
  m_suds.clear();
}

void NLL_Combined_Sudakov::Add(NLL_Single_Sudakov* sud) 
{
  m_suds.push_back(sud);
}

double NLL_Combined_Sudakov::Log(double Q, double q) 
{
  double sum(0.);
  for (size_t i=0;i<m_suds.size();++i) sum += m_suds[i]->Log(Q,q);
  return sum;
}

double NLL_Combined_Sudakov::operator()(double Q, double q) 
{
  if (Q<q) return 1.;
  return exp(-Log(Q,q));
}

double NLL_Combined_Sudakov::operator()(double Q) {
  return operator()(Q,0.);
}

