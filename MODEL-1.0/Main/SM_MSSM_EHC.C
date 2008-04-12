#include "SM_MSSM_EHC.H"
#include "Standard_Model.H"
#include "THDM.H"
#include "MSSM.H"
#include "Effective_Higgs_Coupling.H"
#include "Message.H"
#include "Data_Reader.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(SM_EHC_Getter,"SM+EHC",Model_Base,Model_Arguments);

Model_Base *SM_EHC_Getter::operator()(const Model_Arguments &args) const
{
  return new SM_EHC(args.m_path,args.m_file,args.m_elementary);
}

void SM_EHC_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + Effective Higgs Couplings to Photons/Gluons\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- FINITE_TOP_MASS (0(1) neglect(consider) top mass effects in loop)\n"
     <<std::setw(width+7)<<" "<<"- HIGGS_PP_EFF (Higgs-Photon-Photon coupling constant)\n"
     <<std::setw(width+4)<<" "<<"}\n";
}

SM_EHC::SM_EHC(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model \\w EHC from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+EHC");
  p_numbers   = new ScalarNumbersMap();
  p_constants = new ScalarConstantsMap();
  p_functions = new ScalarFunctionsMap();
  p_matrices  = new ComplexMatricesMap();
 
  Standard_Model * sm = new Standard_Model(m_dir,m_file,false);
  p_numbers   = sm->ExtractScalarNumbers();
  p_constants = sm->ExtractScalarConstants();
  p_functions = sm->ExtractScalarFunctions();
  p_matrices  = sm->ExtractComplexMatrices();

  delete sm;
  
  ParticleInit();
  FillSpectrum();
}

void SM_EHC::ParticleInit() {
  //add pseudo gluon
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[kf_shgluon] = new Particle_Info(kf_shgluon,0.,0.0,0,0,8,2,-1,1,1,0,"shgluon","shgluon",1);
  
  //particle data read by SM
}

void SM_EHC::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
  
  //Effective coupling for Higgs-Gluon-Gluon / Higgs-3 Gluon /Higgs-4 Gluon vertices 
  double eh=2./3.;
  if (p_dataread->GetValue<int>("FINITE_TOP_MASS",0)==1) {
    double hm=Flavour(kf_h0).Mass();
    Effective_Higgs_Coupling ehc(hm);
    eh = ehc.GetFermionContribution(Flavour(kf_t).Mass());
  }
  double hpp = p_dataread->GetValue<double>("HIGGS_PP_EFF",0.);
  p_constants->insert(std::make_pair(std::string("HIGGS_PP_EFF"),hpp));
  p_constants->insert(std::make_pair(std::string("Higgs_gg_fac"),eh));
}


DECLARE_GETTER(MSSM_EHC_Getter,"MSSM+EHC",Model_Base,Model_Arguments);

Model_Base *MSSM_EHC_Getter::operator()(const Model_Arguments &args) const
{
  return new MSSM_EHC(args.m_path,args.m_file,args.m_elementary);
}

void MSSM_EHC_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MSSM + Effective Higgs Coupling to Gluons\n"
    <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the MSSM parameters\n"
     <<std::setw(width+7)<<" "<<"- FINITE_TOP_MASS (0(1) neglect(consider) top mass effects in loop)\n"
     <<std::setw(width+4)<<" "<<"}\n";
}

MSSM_EHC::MSSM_EHC(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  if (m_elementary)
    msg_Info()<<"Initialize the MSSM \\w EHCs from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("MSSM+EHC");
  p_numbers   = new ScalarNumbersMap();
  p_constants = new ScalarConstantsMap();
  p_functions = new ScalarFunctionsMap();
  p_matrices  = new ComplexMatricesMap();
 
  MSSM * mssm = new MSSM(m_dir,m_file,false);
  p_numbers   = mssm->ExtractScalarNumbers();
  p_constants = mssm->ExtractScalarConstants();
  p_functions = mssm->ExtractScalarFunctions();
  p_matrices  = mssm->ExtractComplexMatrices();

  delete mssm;
  
  ParticleInit();
  FillSpectrum();
}

void MSSM_EHC::ParticleInit() {
  //add pseudo gluon
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[kf_shgluon] = new Particle_Info(kf_shgluon,0.,0.0,0,0,8,2,-1,1,1,0,"shgluon","shgluon",1);
  //particle data read by SM
}


void MSSM_EHC::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
  
  //Effective coupling for Higgs-Gluon-Gluon / Higgs-3 Gluon /Higgs-4 Gluon vertices 
  double eh=2./3.;
  if (p_dataread->GetValue<int>("FINITE_TOP_MASS",0)==1) {
    double hm=Flavour(kf_h0).Mass();
    Effective_Higgs_Coupling ehc(hm);
    eh = ehc.GetFermionContribution(Flavour(kf_t).Mass());
  }
  p_constants->insert(std::make_pair(std::string("Higgs_gg_fac"),eh));
}


