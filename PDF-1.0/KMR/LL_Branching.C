#include "LL_Branching.H"

#include "QCD_Splitting_Kernels.H"
#include "Running_AlphaS.H"
#include "MyStrStream.H"
#include "Exception.H"
#include "MathTools.H"

using namespace PDF;
using namespace ATOOLS;

const int PDF::LL_Branching::s_nfmax(5);

LL_Branching::SF_Set PDF::LL_Branching::s_splittings;
LL_Branching::Initializer_Function PDF::LL_Branching::s_initializer; 

void LL_Branching::Initializer_Function::Initialize() const
{
  static bool initialized=false;
  if (initialized) return;
  for (int i=1;i<=s_nfmax;++i) {
    Flavour cur((kf::code)i);
    s_splittings.insert(new Q_QG(cur));
    s_splittings.insert(new Q_QG(cur.Bar()));
    s_splittings.insert(new Q_GQ(cur));
    s_splittings.insert(new Q_GQ(cur.Bar()));
    s_splittings.insert(new G_QQ(cur));
    s_splittings.insert(new G_QQ(cur.Bar()));
  }
  // account for factor 2 in P_{gg}(z)
  s_splittings.insert(new G_GG());
  s_splittings.insert(new G_GG());
  initialized=true;
}

LL_Branching::Initializer_Function::~Initializer_Function() 
{
  while (s_splittings.size()>0) {
    delete *s_splittings.begin();
    s_splittings.erase(s_splittings.begin());
  }
}

LL_Branching::LL_Branching(const Flavour flavour,
			   MODEL::Running_AlphaS *alphas):
  m_flavour(flavour),
  p_alphas(alphas)
{
  s_initializer.Initialize();
  m_splittings.push_back(*Find(m_flavour,m_flavour));
  if (m_flavour.IsGluon()) {
    for (int i=1;i<=s_nfmax;++i) {
      m_splittings.push_back(*Find(m_flavour,(kf::code)i));
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
  sstr<<m_splittings[0]->GetA()<<"_as("
      <<Flavour((kf::code)24).PSMass()<<")-"
      <<p_alphas->AsMZ()<<"_o-"<<p_alphas->Order()<<".dat";
  sstr>>temp;
  m_name+=temp;
}

std::set<Splitting_Kernel*>::iterator 
LL_Branching::Find(const Flavour &a,const Flavour &b)
{
  std::set<Splitting_Kernel*>::iterator sfit;
  for (sfit=s_splittings.begin();sfit!=s_splittings.end();++sfit)
    if ((*sfit)->GetA()==a && (*sfit)->GetB()==b) return sfit;
  return s_splittings.end();
}

double LL_Branching::Gamma(double q, double Q)
{
  double splitting(0.0);
  if (m_flavour.IsQuark()) {
    splitting=m_splittings[0]->Integral(0.0,Q/(Q+q));
  }
  else if (m_flavour.IsGluon()) {
    splitting=m_splittings[0]->Integral(q/(Q+q),Q/(Q+q));
    size_t nf(Min(s_nfmax,MODEL::as->Nf(q*q)));
    for (size_t i(1);i<=nf;++i) splitting+=m_splittings[i]->Integral(0.0,1.0);
  }
  else {
    THROW(fatal_error,"Invalid flavour.");
  }
  return (*p_alphas)(q*q)/(M_PI*q)*splitting;
}

double LL_Branching::IntGamma(double q, double Q)
{
  return -1.; 
}


