#include "LL_Sudakov.H"

#include "Run_Parameter.H"
#include "MyStrStream.H"

#define Nc 3
#define CA Nc
#define Nf 4
#define BETA0 (11.*CA-2.*Nf)/3.

using namespace SHERPA;

LL_Sudakov::LL_Sudakov(MODEL::Running_AlphaS *alphas):
  p_alphas(alphas)
{
  FixLambda2();
}

LL_Sudakov::~LL_Sudakov()
{
  for (size_t i=0; i<m_all_suds.size();++i) delete m_all_suds[i];
  m_all_suds.clear();
  m_sud_map.clear();
}

void LL_Sudakov::FixLambda2() 
{
  m_mu2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass());
  m_asmu=(*MODEL::as)(m_mu2);
  m_lambda=sqrt(m_mu2*exp(-4.*M_PI/(BETA0*m_asmu)));
}                 

void LL_Sudakov::AssignKeys(ATOOLS::Integration_Info *const info) 
{
  for (std::vector<NLL_Sudakov_Base*>::iterator sit=m_all_suds.begin();
       sit!=m_all_suds.end();++sit) {
    (*sit)->AssignKey(info);
  }
}

void LL_Sudakov::Initialize()
{
  int smode=Sudakov::numeric;
  LL_Single_Sudakov *ssud=NULL;
  for (int k=1;k<=Nf;++k) {
    ATOOLS::Flavour fl=ATOOLS::Flavour((ATOOLS::kf::code)k);
    ssud = new LL_Single_Sudakov(new PDF::LL_Branching(fl,p_alphas),smode);
    m_sud_map[fl]=ssud;
    m_sud_map[fl.Bar()]=ssud;
    m_all_suds.push_back(ssud);
  }
  ssud = new LL_Single_Sudakov(new PDF::LL_Branching(ATOOLS::Flavour(ATOOLS::kf::gluon),
						     p_alphas),smode);
  m_all_suds.push_back(ssud);
  m_sud_map[ATOOLS::Flavour(ATOOLS::kf::gluon)]=ssud;
}

NLL_Sudakov_Base &LL_Sudakov::Delta(const ATOOLS::Flavour &fl) 
{
  Sudakov_Map::const_iterator sit=m_sud_map.find(fl);
  if (sit!=m_sud_map.end()) {
    return *(sit->second);
  }
  ATOOLS::msg.Error()<<"LL_Sudakov::Delta("<<fl<<"): "<<ATOOLS::om::red
		     <<"Did not find sudakov !"<<ATOOLS::om::reset<<std::endl;
  return *m_sud_map[ATOOLS::Flavour(ATOOLS::kf::none)];
}

double LL_Sudakov::Delta(const ATOOLS::Flavour &fl,double Q,double q) 
{
  return Delta(fl)(Q,q);
}

void LL_Sudakov::SetOutPath(const std::string outpath) 
{
  m_outpath=outpath;
}


