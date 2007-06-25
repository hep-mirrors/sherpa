#include "Interaction_Model_Handler.H"

#include "Interaction_Model_QCD.H"
#include "Interaction_Model_EW.H"
#include "Interaction_Model_SM.H"
#include "Interaction_Model_THDM.H"
#include "Interaction_Model_MSSM.H"
#include "Interaction_Model_ADD.H"
#include "Interaction_Model_SMEHC.H"
#include "Interaction_Model_MSSM_EHC.H"
#include "Interaction_Model_SM_ZPrime.H"
#include "Interaction_Model_SM_AGC.H"
#include "Interaction_Model_SM_Phantom_U1.H"
#include "Interaction_Model_4GenLep.H"

#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_Handler::Interaction_Model_Handler(MODEL::Model_Base * model) :
  p_model(model)
{ }

Interaction_Model_Base * Interaction_Model_Handler::GetModel(std::string modeltype,
							     std::string cplscheme,
							     std::string yukscheme)
{
  if (modeltype==std::string("pure_QCD")||modeltype==std::string("QCD")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::pure_QCD);
    return new Interaction_Model_QCD(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("pure_EW")||modeltype==std::string("EW")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::pure_EW);
    return new Interaction_Model_EW(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("SM"))  {
    rpa.gen.SetModelType(ATOOLS::Model_Type::SM);
    return new Interaction_Model_SM(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("THDM")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::THDM);
    return new Interaction_Model_THDM(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("MSSM")) {
      rpa.gen.SetModelType(ATOOLS::Model_Type::MSSM);
      return new Interaction_Model_MSSM(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("MSSM+EHC")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::MSSM_EHC);
    return new Interaction_Model_MSSM_EHC(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("ADD")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::ADD);
    return new Interaction_Model_ADD(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("SM+EHC")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::SMEHC);
    return new Interaction_Model_SMEHC(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("SM+ZPrime")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::SM_ZPrime);
    return new Interaction_Model_SM_ZPrime(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("SM+AGC")) {
    rpa.gen.SetModelType(ATOOLS::Model_Type::SM_AGC);
    return new Interaction_Model_SM_AGC(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("SM+Phantom_U1"))  {
    rpa.gen.SetModelType(ATOOLS::Model_Type::SM_Phantom_U1);
    return new Interaction_Model_SM_Phantom_U1(p_model,cplscheme,yukscheme);
  }
  if (modeltype==std::string("FOURTH_GEN_LEPTONS"))  {
    rpa.gen.SetModelType(ATOOLS::Model_Type::FourthGen_Lep);
    return new Interaction_Model_4GenLep(p_model,cplscheme,yukscheme);
  }

  msg.Error()<<"Error in Interaction_Model_Handler::GetModel("<<modeltype<<") : "<<endl
	     <<"   Model not found. Initialize Standard Model."<<endl;
  rpa.gen.SetModelType(ATOOLS::Model_Type::SM);
  return new Interaction_Model_SM(p_model,cplscheme,yukscheme);
}
