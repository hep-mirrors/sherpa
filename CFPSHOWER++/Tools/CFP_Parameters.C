#include "CFPSHOWER++/Tools/CFP_Parameters.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;


CFP_Parameters * CFPSHOWER::cfp_pars = NULL;

CFP_Parameters::CFP_Parameters() {}

bool CFP_Parameters::Init(Default_Reader *const reader)
{
  m_switches["kinematics"]     = reader->Get<int>("CSS_KIN_SCHEME",1);
  m_switches["kfactor"]        = reader->Get<int>("CSS_KFACTOR_SCHEME",1);
  m_switches["couplings"]      = reader->Get<int>("CSS_COUPLING_SCHEME",1);
  m_switches["ME_corrections"] = reader->Get<int>("CSS_ME_CORRECTION",0);
  m_switches["max_emissions"]  = reader->Get<unsigned int>
    ("CSS_MAXEM",100000); //std::numeric_limits<unsigned int>::max());
  m_switches["max_particles"]  = reader->Get<unsigned int>
    ("CSS_MAXPART",100000); //std::numeric_limits<unsigned int>::max());
  m_parameters["recalc_fac"]   = reader->Get<double>("CSS_RECALC_FACTOR",4.0);
  m_parameters["PDF_min"]      = reader->Get<double>("CSS_PDF_MIN",1.0e-4);
  m_parameters["PDF_min_X"]    = reader->Get<double>("CSS_PDF_MIN_X",1.0e-2);
  m_parameters["pt2min(FS)"]   = ToType<double>(rpa->gen.Variable("CSS_FS_PT2MIN"));
  m_parameters["pt2min(IS)"]   = ToType<double>(rpa->gen.Variable("CSS_IS_PT2MIN"));
  m_parameters["k_alpha(FS)"]  = ToType<double>(rpa->gen.Variable("CSS_FS_AS_FAC"));
  m_parameters["k_alpha(IS)"]  = ToType<double>(rpa->gen.Variable("CSS_IS_AS_FAC"));
  m_parameters["k_muR"]        = ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  m_parameters["k_muF"]        = ToType<double>(rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR"));
  /*
  m_rcf=reader->Get<double>("CSS_RECALC_FACTOR",4.0);
  m_tcef=reader->Get<double>("CSS_TC_ENHANCE",1.0);
  m_maxrewem=reader->Get<unsigned int>
    ("REWEIGHT_MAXEM",std::numeric_limits<unsigned int>::max());
  m_rewtmin=reader->Get<double>("CSS_REWEIGHT_SCALE_CUTOFF", 5.0);
  m_oef=reader->Get<double>("CSS_OEF",3.0);
  */
  Output();
  return true;
}

int CFP_Parameters::operator[](string keyword) 
{
  map<string,int>::iterator siter = m_switches.find(keyword);
  if (siter!=m_switches.end()) return siter->second;
  msg_Error()<<"Error in CFP_Switches("<<keyword<<"), "
	     <<"in total "<<m_switches.size()<<" switches.\n"
	     <<"   Keyword not found. Return 0 and hope for the best.\n";
  exit(1);
  return 0;
}

double CFP_Parameters::operator()(string keyword) 
{
  map<string,double>::iterator piter = m_parameters.find(keyword);
  if (piter!=m_parameters.end()) return piter->second;
  msg_Error()<<"Error in CFP_Parameters("<<keyword<<") "
	     <<"in "<<m_parameters.size()<<".\n"
	     <<"   Keyword not found. Return 0 and hope for the best.\n";
  exit(1);
  return 0.;
}

void CFP_Parameters::Output() {
  msg_Out()<<"*********************************************\n"
	   <<METHOD<<" has "<<m_switches.size()<<" switches and "
	   <<m_parameters.size()<<" parameters:\n";
  for (map<string,int>::iterator siter=m_switches.begin();
       siter!=m_switches.end();siter++) {
    msg_Out()<<"  *** "<<siter->first<<" = "<<siter->second<<" (switch)\n";
  }
  for (map<string,double>::iterator piter=m_parameters.begin();
       piter!=m_parameters.end();piter++) {
    msg_Out()<<"  *** "<<piter->first<<" = "<<piter->second<<" (parameter)\n";
  }
  msg_Out()<<"*********************************************\n";
}
