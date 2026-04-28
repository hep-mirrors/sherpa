#include "BEAM/Ions/Ion_Parameters.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"


using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Ion_Parameters* BEAM::ionpars = nullptr;


void Ion_Parameters::Initialise() {
  auto s = Settings::GetMainSettings()["IONS"];
  m_parametermap["K_Fermi"] =
    s["K_Fermi"].SetDefault(0.2).Get<double>();
  m_parametermap["Rho_0"] =
    s["Rho_0"].SetDefault(0.17).Get<double>();
  m_parametermap["R_min^2(Pauli)"] =
    sqr(s["R_min(Pauli)"].SetDefault(1.5).Get<double>());
  m_parametermap["Sigma_Psi"] =
    s["Sigma_Psi"].SetDefault(1.03).Get<double>();
  m_parametermap["V^(loc)_2"] =
    s["V^(loc)_2"].SetDefault(-0.1).Get<double>();
  m_parametermap["V^(loc)_3"] =
    s["V^(loc)_3"].SetDefault(0.0).Get<double>();
  m_parametermap["R^(loc)_max"] =
    s["R^(loc)_max"].SetDefault(20.).Get<double>();
  m_parametermap["nu^(loc)"] =
    s["nu^(loc)"].SetDefault(2.).Get<double>();
  m_parametermap["a^(Yuk)"] =
    s["a^(Yuk)"].SetDefault(0.8).Get<double>();
  m_parametermap["V^(Yuk)"] =
    s["V^(Yuk)"].SetDefault(0.0).Get<double>();
  m_parametermap["R^(Yuk)_max"] =
    s["R^(Yuk)_max"].SetDefault(20.).Get<double>();
  m_parametermap["delta_t"] =
    s["delta_t"].SetDefault(0.1).Get<double>();
  m_integermap["Evolver"] =
    s["Evolver"].SetDefault(0).Get<int>();
  m_integermap["R_Form"] =
    s["R_Form"].SetDefault(1).Get<int>();
  m_integermap["P_Form"] =
    s["P_Form"].SetDefault(1).Get<int>();
  m_integermap["N_steps"] =
    s["N_steps"].SetDefault(1000).Get<int>();
}

const double Ion_Parameters::Get(string keyword) {
  map<string,double>::const_iterator piter = m_parametermap.find(keyword);
  if (piter!=m_parametermap.end()) return piter->second;
  msg_Error()<<"Error in Ion_Parameters::Get("<<keyword<<") "
	     <<"in "<<m_parametermap.size()<<".\n"
	     <<"   Keyword not found. Return 0 and hope for the best.\n";
  return 0.;
}

const int Ion_Parameters::GetInt(string keyword) {
  map<string,int>::const_iterator piter = m_integermap.find(keyword);
  if (piter!=m_integermap.end()) return piter->second;
  msg_Error()<<"Error in Ion_Parameters::Get("<<keyword<<") "
	     <<"in "<<m_integermap.size()<<".\n"
	     <<"   Keyword not found. Return -1 and hope for the best.\n";
  for (map<string,int>::const_iterator piter = m_integermap.begin();
       piter!=m_integermap.end();piter++)
    msg_Out()<<piter->first<<": "<<piter->second<<"\n";
  return -1;
}
