#include "LL_Branching.H"

#include "QCD_Splitting_Functions.H"
#include "MyStrStream.H"
#include "Exception.H"

#define NC 3.
#define CF (NC*NC-1.)/(2.*NC)
 
using namespace PDF;

const unsigned int PDF::LL_Branching::s_nf=3;

LL_Branching::SF_Set PDF::LL_Branching::s_splitting=LL_Branching::SF_Set();
LL_Branching PDF::LLB;

void LL_Branching::DeleteSplittings()
{
  while (s_splitting.size()>0) {
    delete *s_splitting.begin();
    s_splitting.erase(s_splitting.begin());
  }
}

LL_Branching::LL_Branching() {}

LL_Branching::LL_Branching(const ATOOLS::Flavour flavour,MODEL::Running_AlphaS *alphas):
  m_flavour(flavour),
  p_alphas(alphas)
{
  APACIC::Splitting_Function *newsplit;
  m_splitting.push_back(new APACIC::Splitting_Group());
  SF_Set::iterator sfit;
  if (m_flavour.IsQuark()) {
    if ((sfit=Find(m_flavour,ATOOLS::Flavour(ATOOLS::kf::gluon)))==s_splitting.end()) {
      Insert(new APACIC::q_gq(m_flavour));
    }
    if ((sfit=Find(m_flavour,m_flavour))==s_splitting.end()) {
      newsplit = new APACIC::q_qg(m_flavour);
      Insert(newsplit);
    }
    else {
      newsplit=*sfit;
    }
    m_splitting[0]->Add(newsplit);
  }
  else if (m_flavour.IsGluon()) {
    if ((sfit=Find(m_flavour,m_flavour))==s_splitting.end()) {
      newsplit = new APACIC::g_gg();
      Insert(newsplit);
    }
    else {
      newsplit=*sfit;
    }
    m_splitting[0]->Add(newsplit);
    m_splitting.push_back(new APACIC::Splitting_Group());
    for (int i=-(1+s_nf);i<(int)s_nf;++i) {
      if (i!=0) {
	if ((sfit=Find(m_flavour,(ATOOLS::kf::code)i+1))==s_splitting.end()) {
	  newsplit = new APACIC::g_qq((ATOOLS::kf::code)i+1);
	  Insert(newsplit);
	}
	else {
	  newsplit=*sfit;
	}
	m_splitting[1]->Add(newsplit);
	m_const.push_back(m_splitting[1]->Integral(0.,1.));
      }
    }
  }
  else {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"No splitting function available.",
			    "LL_Branching","LL_Branching"));
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
  sstr<<m_splitting[0]->GetFlA()<<"_as("
      <<ATOOLS::Flavour((ATOOLS::kf::code)24).PSMass()<<")-"
      <<p_alphas->AsMZ()<<"_o-"<<p_alphas->Order()<<".dat";
  sstr>>temp;
  m_name+=temp;
}

std::set<APACIC::Splitting_Function*>::iterator 
LL_Branching::Find(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
{
  std::set<APACIC::Splitting_Function*>::iterator sfit;
  for (sfit=s_splitting.begin();sfit!=s_splitting.end();++sfit) {
    if ((*sfit)->GetFlA()==a && (*sfit)->GetFlB()==b) {
      return sfit;
    }
  }
  return s_splitting.end();
}

void LL_Branching::Insert(APACIC::Splitting_Function *splitting) 
{
  bool found=false;
  std::set<APACIC::Splitting_Function*>::iterator sfit;
  for (sfit=s_splitting.begin();sfit!=s_splitting.end();++sfit) {
    if ((*sfit)->GetFlA()==splitting->GetFlA() &&
	(*sfit)->GetFlB()==splitting->GetFlB() &&
	(*sfit)->GetFlC()==splitting->GetFlC()) {
      found=true;
    }
  }
  if (!found) s_splitting.insert(splitting);
}

double LL_Branching::Gamma(double q, double Q)
{
  double splitting=0.;
  if (m_flavour.IsQuark()) {
    splitting=m_splitting[0]->Integral(0.,Q/(Q+q));
  }
  else if (m_flavour.IsGluon()) {
    //splitting=p_alphas->Nf(q*q)*m_const[0]/s_nf;
    splitting=m_const[0];
    splitting+=m_splitting[0]->Integral(q/(Q+q),Q/(Q+q));
  }
  return (*p_alphas)(q*q)/(M_PI*q)*splitting;
}

double LL_Branching::IntGamma(double q, double Q)
{
  return -1.; 
}

const LL_Branching::SF_Set &LL_Branching::AllSplittings() const
{
  return s_splitting; 
}

