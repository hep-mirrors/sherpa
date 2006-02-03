#include "NLL_Sudakov.H"
#include "Message.H"
#include "MathTools.H"
#include "Running_AlphaS.H"
#include <iomanip> 
#include <stdio.h>

#include "NLL_Single_Sudakov.H"
#include "NLL_Combined_Sudakov.H"
#include "NLL_Branching_Probabilities.H"
#include "NLL_JetRate.H"
#include "MyStrStream.H"

namespace SHERPA {
  const double   Nc    = 3;
  const double   CA    = Nc;
  const double   CF    = (Nc*Nc-1.)/(2.*Nc);
  const double   Nf    = 5;
  const double   TR    =  1./2.;
  const double   BETA0 = (11.*CA-2.*Nf)/3.;
  const double   BETA1 = (17.*CA*CA- 3.*CF*Nf-5.*CA*Nf)/3.;
}

using namespace SHERPA;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

NLL_Sudakov::NLL_Sudakov(int mode, double _tmax,double _tmin,
			 MODEL::Running_AlphaS * runas,double asfac): 
  m_mode(mode), p_runas(runas), m_as_factor(asfac)
{
  FixLambda2();
  PrepareMap();
}

NLL_Sudakov::~NLL_Sudakov()
{
  std::set<NLL_Sudakov_Base*> deleted;
  for (Sudakov_Map::iterator sud=m_sudakovs.begin();
       sud!=m_sudakovs.end();sud++) {
    if (sud->second!=NULL) { 
      if (deleted.find(sud->second)==deleted.end()) {
	delete sud->second; 
	deleted.insert(sud->second);
      }
      sud->second=NULL; 
    }
  }
  m_sudakovs.clear();
}

void NLL_Sudakov::PrepareMap() 
{
  int  smode=Sudakov::numeric;
  NLL_Dummy_Sudakov * dsud(new NLL_Dummy_Sudakov());
  m_sudakovs[Flavour(kf::none)] = dsud;
  NLL_Combined_Sudakov *csud(NULL);
  NLL_Single_Sudakov   *ssud(NULL);
  Flavour flav;
  bpm::code bpmode((bpm::code)(m_mode));
  if (bpmode & bpm::fs) {
    for (int k=1;k<=6;++k) {
      flav=Flavour(kf::code(k));
      ssud = new NLL_Single_Sudakov
	(new GammaQ_QG_Lambda
	 (bpmode,m_lambda,p_runas,flav.Mass(),m_as_factor),smode);
      m_sudakovs[flav]=ssud;
      m_sudakovs[flav.Bar()]=ssud;
    }
    csud = new NLL_Combined_Sudakov(smode);
    ssud = new NLL_Single_Sudakov
      (new GammaG_GG_Lambda(bpmode,m_lambda,p_runas,m_as_factor),smode);
    csud->Add(ssud);
    for (int k=1;k<=6;++k) {
      flav=Flavour(kf::code(k));
      ssud = new NLL_Single_Sudakov
	(new GammaG_QQ_Lambda
	 (bpmode,m_lambda,p_runas,flav.Mass(),m_as_factor),smode);
      csud->Add(ssud);
    }
    m_sudakovs[Flavour(kf::gluon)]=csud;
  }
  if (bpmode & bpm::is) {
    csud = new NLL_Combined_Sudakov(smode);
    ssud = new NLL_Single_Sudakov
      (new GammaQ_QG_Lambda(bpmode,m_lambda,p_runas,m_as_factor),smode);
    csud->Add(ssud);
//     ssud = new NLL_Single_Sudakov
//       (new GammaG_QQ_Lambda(bpmode,m_lambda,p_runas,m_as_factor),smode);
//     csud->Add(ssud);
    for (int k=1;k<=5;++k) {
      flav                   = Flavour(kf::code(k));
      m_sudakovs[flav]       = csud;
      m_sudakovs[flav.Bar()] = csud;
    }
    csud = new NLL_Combined_Sudakov(smode);
    ssud = new NLL_Single_Sudakov
      (new GammaG_GG_Lambda(bpmode,m_lambda,p_runas,m_as_factor),smode);
    csud->Add(ssud);
    csud->Add(ssud);
//     ssud = new NLL_Single_Sudakov
//       (new GammaQ_GQ_Lambda(bpmode,m_lambda,p_runas,m_as_factor),smode);
//     for (int k=1;k<=5;++k) csud->Add(ssud);
    m_sudakovs[Flavour(kf::gluon)] = csud;
  }
}

void NLL_Sudakov::FixLambda2() 
{
  m_mu2    = sqr(Flavour(kf::Z).Mass());
  m_asmu   = (*as)(m_mu2);
  m_lambda = sqrt( m_mu2 * exp(-4.*M_PI/(BETA0 * m_asmu)));
}                 

NLL_Sudakov_Base &  NLL_Sudakov::Delta(const ATOOLS::Flavour & fl) {
  Sudakov_Map::const_iterator sit=m_sudakovs.find(fl);
  if (sit!=m_sudakovs.end()) return *(sit->second);
  ATOOLS::msg.Out()<<"ERROR in  NLL_Sudakov::Delta : "<<std::endl
		   <<"   Did not find sudakov form factor for "<<fl<<", return default."<<std::endl;
  return *(m_sudakovs[ATOOLS::Flavour(ATOOLS::kf::none)]);
}
