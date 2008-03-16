#include "NLL_Sudakov.H"
#include "Message.H"
#include "MathTools.H"
#include "Running_AlphaS.H"
#include "Model_Base.H"
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

NLL_Sudakov::NLL_Sudakov(int mode,MODEL::Running_AlphaS *runas,
			 double asfac):
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
  m_sudakovs[Flavour(kf_none)] = dsud;
  NLL_Combined_Sudakov *csud(NULL);
  NLL_Single_Sudakov   *ssud(NULL);
  Flavour flav;
  bpm::code bpmode((bpm::code)(m_mode));
  //Quark Sudakovs
  for (int k=1;k<=6;++k) {
    flav=Flavour((kf_code)(k));
    ssud = new NLL_Single_Sudakov
      (new GammaQ_QG_Lambda
       (bpmode,m_lambda,p_runas,flav.PSMass(),m_as_factor),smode);
    m_sudakovs[flav]=ssud;
    m_sudakovs[flav.Bar()]=ssud;
  }
  
  if (MODEL::s_model->Name()==std::string("MSSM")) {
    //Gluino Sudakov
    flav=Flavour(kf_Gluino);
    GammaQ_QG_Lambda * GL = new GammaQ_QG_Lambda
      (bpmode,m_lambda,p_runas,flav.PSMass(),m_as_factor);
    GL->SetColFac(CA);
    ssud = new NLL_Single_Sudakov(GL,smode);
    m_sudakovs[flav]=ssud;
    
    //sQuark Sudakovs
    for (short int l=1;l<3;l++) {
      for (short int i=1;i<7;i++) {
	int fl = l*1000000 + i;
	flav = Flavour((kf_code)(fl));
	std::cout<<" flav is : "<<flav<<" mass is "<<flav.PSMass()<<std::endl;  
	ssud = new NLL_Single_Sudakov
	  (new GammasQ_sQG_Lambda
	   (bpmode,m_lambda,p_runas,flav.PSMass(),m_as_factor),smode);
	m_sudakovs[flav]=ssud;
	m_sudakovs[flav.Bar()]=ssud;
      }
    }
  }
  
  //Gluon Sudakov
  csud = new NLL_Combined_Sudakov(smode);
  ssud = new NLL_Single_Sudakov
    (new GammaG_GG_Lambda(bpmode,m_lambda,p_runas,m_as_factor),smode);
  csud->Add(ssud);
  for (int k=1;k<=6;++k) {
    flav=Flavour((kf_code)(k));
    ssud = new NLL_Single_Sudakov
      (new GammaG_QQ_Lambda
       (bpmode,m_lambda,p_runas,flav.PSMass(),m_as_factor),smode);
    csud->Add(ssud);
  }
  m_sudakovs[Flavour(kf_gluon)]=csud;
}

void NLL_Sudakov::FixLambda2() 
{
  m_mu2    = sqr(Flavour(kf_Z).Mass());
  m_asmu   = (*as)(m_mu2);
  m_lambda = sqrt( m_mu2 * exp(-4.*M_PI/(BETA0 * m_asmu)));
}                 

NLL_Sudakov_Base &  NLL_Sudakov::Delta(const ATOOLS::Flavour & fl) {
  Sudakov_Map::const_iterator sit=m_sudakovs.find(fl);
  if (sit!=m_sudakovs.end()) return *(sit->second);
  msg_Out()<<"ERROR in  NLL_Sudakov::Delta : "<<std::endl
		   <<"   Did not find sudakov form factor for "<<fl<<", return default."<<std::endl;
  return *(m_sudakovs[ATOOLS::Flavour(kf_none)]);
}
