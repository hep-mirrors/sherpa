#include "LL_Branching.H"

#include "QCD_Splitting_Functions.H"
#include "Running_AlphaS.H"
#include "MyStrStream.H"
#include "Exception.H"

using namespace PDF;

const double PDF::LL_Branching::s_nfmax=5.;
const double PDF::LL_Branching::s_nc=3.;
const double PDF::LL_Branching::s_cf=(3.*3.-1.)/(2.*3.);
const double PDF::LL_Branching::s_ca=3.;

LL_Branching::SF_Set PDF::LL_Branching::s_splittings;
LL_Branching::Initializer_Function PDF::LL_Branching::s_initializer; 

void LL_Branching::Initializer_Function::Initialize() const
{
  static bool initialized=false;
  if (initialized) return;
  for (int i=1;i<=s_nfmax;++i) {
    ATOOLS::Flavour cur((ATOOLS::kf::code)i);
    Insert(new APACIC::q_qg(cur));
    Insert(new APACIC::q_qg(cur.Bar()));
    Insert(new APACIC::q_gq(cur));
    Insert(new APACIC::q_gq(cur.Bar()));
    Insert(new APACIC::g_qq(cur));
    Insert(new APACIC::g_qq(cur.Bar()));
  }
  Insert(new APACIC::g_gg());
  initialized=true;
}

LL_Branching::Initializer_Function::~Initializer_Function() 
{
  while (s_splittings.size()>0) {
    delete *s_splittings.begin();
    s_splittings.erase(s_splittings.begin());
  }
}

LL_Branching::LL_Branching(const ATOOLS::Flavour flavour,
			   MODEL::Running_AlphaS *alphas):
  m_flavour(flavour),
  p_alphas(alphas)
{
  s_initializer.Initialize();
  m_splittings.push_back(new APACIC::Splitting_Group());
  m_splittings[0]->Add(*Find(m_flavour,m_flavour));
  if (m_flavour.IsGluon()) {
    m_splittings.push_back(new APACIC::Splitting_Group());
    for (int i=1;i<=(int)s_nfmax;++i) {
      m_splittings[1]->Add(*Find(m_flavour,(ATOOLS::kf::code)i));
    }
  }
  else if (!m_flavour.IsQuark()) {
    THROW(fatal_error,"No splitting function available.");
  }
  GenerateName();
}

LL_Branching::~LL_Branching() {}

void LL_Branching::GenerateName() 
{
  MyStrStream sstr;
  std::string temp;
  m_name=std::string("lls_");
  sstr.clear();
  sstr<<m_splittings[0]->GetFlA()<<"_as("
      <<ATOOLS::Flavour((ATOOLS::kf::code)24).PSMass()<<")-"
      <<p_alphas->AsMZ()<<"_o-"<<p_alphas->Order()<<".dat";
  sstr>>temp;
  m_name+=temp;
}

std::set<APACIC::Splitting_Function*>::iterator 
LL_Branching::Find(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
{
  std::set<APACIC::Splitting_Function*>::iterator sfit;
  for (sfit=s_splittings.begin();sfit!=s_splittings.end();++sfit) {
    if ((*sfit)->GetFlA()==a && (*sfit)->GetFlB()==b) {
      return sfit;
    }
  }
  return s_splittings.end();
}

void LL_Branching::Insert(APACIC::Splitting_Function *splitting) 
{
  bool found=false;
  std::set<APACIC::Splitting_Function*>::iterator sfit;
  for (sfit=s_splittings.begin();sfit!=s_splittings.end();++sfit) {
    if ((*sfit)->GetFlA()==splitting->GetFlA() &&
	(*sfit)->GetFlB()==splitting->GetFlB() &&
	(*sfit)->GetFlC()==splitting->GetFlC()) {
      found=true;
    }
  }
  if (!found) s_splittings.insert(splitting);
}

double LL_Branching::Gamma(double q, double Q)
{
  double splitting=0.;
  if (m_flavour.IsQuark()) {
    splitting=m_splittings[0]->Integral(0.,Q/(Q+q));
  }
  else if (m_flavour.IsGluon()) {
    splitting=MODEL::as->Nf(q*q)/s_nfmax*m_splittings[1]->Integral(0.,1.);
    splitting+=m_splittings[0]->Integral(q/(Q+q),Q/(Q+q));
  }
  return (*p_alphas)(q*q)/(M_PI*q)*splitting;
}

double LL_Branching::IntGamma(double q, double Q)
{
  return -1.; 
}


