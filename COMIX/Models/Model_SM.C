#include "COMIX/Main/Model.H"
#include "ATOOLS/Org/Exception.H"

namespace COMIX {

  class Model_SM: public Model {
  public:

    // constructor
    Model_SM();
    
    // member functions
    void Initialize(MODEL::Model_Base *const model,
		    const std::string &file);

    std::vector<std::string>     IncludedModels() const;
    std::vector<ATOOLS::Flavour> IncludedFlavours() const;

  };// end of class Model_SM

}// end of namespace COMIX

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

#define M_I Complex(0.0,1.0)

using namespace COMIX;
using namespace MODEL;
using namespace ATOOLS;

Model_SM::Model_SM(): Model("SM") {}

void Model_SM::Initialize(MODEL::Model_Base *const model,
			       const std::string &file)
{
  Model::Initialize(model,file);
  if (model->GetScalarNumbers()->count("EW_SCHEME")>0 &&
      model->ScalarNumber("EW_SCHEME")!=0) {
    THROW(not_implemented, "Comix only implements EW_SCHEME=0 so far.");
  }
  double ecms2(sqr(rpa->gen.Ecms()));
  Complex mw(ScalarConstant("MW")), gw(ScalarConstant("GammaW"));
  Complex mz(ScalarConstant("MZ")), gz(ScalarConstant("GammaZ"));
  m_consts["m_{W,PS}^2"]=mw*mw-M_I*gw*mw;
  m_consts["m_{Z,PS}^2"]=mz*mz-M_I*gz*mz;
  m_consts["v_{EW}"]=ScalarConstant("vev");
  if (m_widthscheme!="CMS") {
    m_consts["\\sin\\theta_W"]=sqrt(ScalarConstant("sin2_thetaW"));
    m_consts["\\cos\\theta_W"]=sqrt(1.0-csqr(m_consts["\\sin\\theta_W"]));
    m_consts["m_W"]=sqrt(m_consts["m_W^2"]=mw*mw);
    m_consts["m_Z"]=sqrt(m_consts["m_Z^2"]=mz*mz);
  }
  else {
    msg_Info()<<METHOD<<"(): Initializing complex mass scheme.\n";
    m_consts["m_W"]=sqrt(m_consts["m_W^2"]=m_consts["m_{W,PS}^2"]);
    m_consts["m_Z"]=sqrt(m_consts["m_Z^2"]=m_consts["m_{Z,PS}^2"]);
    m_consts["\\cos\\theta_W"]=sqrt(m_consts["m_W^2"]/m_consts["m_Z^2"]);
    m_consts["\\sin\\theta_W"]=sqrt(1.0-csqr(m_consts["\\cos\\theta_W"]));
  }
  m_consts["g_1"]=sqrt(4.*M_PI*ScalarFunction("alpha_QED",ecms2));
  m_consts["g_2"]=m_consts["g_1"]/m_consts["\\sin\\theta_W"];
  m_consts["g_3"]=sqrt(4.*M_PI*ScalarFunction("alpha_S",ecms2));  
  msg_Debugging()<<METHOD<<"(): Initialized SM parameters {\n";
  msg_Debugging()<<"  m_W = "<<m_consts["m_W"]<<" ( m = "<<mw<<", w = "<<gw
		 <<" ), m_W^2 = "<<m_consts["m_W^2"]<<"\n";
  msg_Debugging()<<"  m_Z = "<<m_consts["m_Z"]<<" ( m = "<<mz<<", w = "<<gz
		 <<" ), m_Z^2 = "<<m_consts["m_Z^2"]<<"\n";
  msg_Debugging()<<"  \\sin\\theta_W = "<<m_consts["\\sin\\theta_W"]
		 <<" ( "<<sqrt(ScalarConstant("sin2_thetaW"))<<" )\n";
  msg_Debugging()<<"  \\cos\\theta_W = "<<m_consts["\\cos\\theta_W"]
		 <<" ( "<<sqrt(1.0-ScalarConstant("sin2_thetaW"))<<" )\n";
  msg_Debugging()<<"  v = "<<m_consts["v_{EW}"]<<"\n";
  msg_Debugging()<<"  g_1 = "<<m_consts["g_1"]<<" => 1/\\alpha_qed = "
		 <<(4.0*M_PI)/sqr(m_consts["g_1"])<<"\n";
  msg_Debugging()<<"  g_2 = "<<m_consts["g_2"]<<"\n";
  msg_Debugging()<<"  g_3 = "<<m_consts["g_3"]<<" => \\alpha_s = "
		 <<sqr(m_consts["g_3"])/(4.0*M_PI)<<"\n";
  msg_Debugging()<<"}\n";
}

std::vector<std::string> Model_SM::IncludedModels() const
{
  std::vector<std::string> models;
  AddModel(models,"QED");
  AddModel(models,"EW");
  AddModel(models,"QCD");
  return models;
}

std::vector<ATOOLS::Flavour> Model_SM::IncludedFlavours() const
{
  std::vector<ATOOLS::Flavour> fls;
  for (int i(1);i<=6;++i) AddFlavour(fls,Flavour((kf_code)i));
  for (int i(11);i<=16;++i) AddFlavour(fls,Flavour((kf_code)i));
  AddFlavour(fls,Flavour(kf_photon));
  AddFlavour(fls,Flavour(kf_Wplus));
  AddFlavour(fls,Flavour(kf_Z));
  AddFlavour(fls,Flavour(kf_h0));
  AddFlavour(fls,Flavour(kf_gluon));
  if (Flavour(kf_Wplus).IsOn()) {
    AddFlavour(fls,Flavour(kf_Z_qgc));
    if (Flavour(kf_photon).IsOn() ||
	Flavour(kf_Z).IsOn()) AddFlavour(fls,Flavour(kf_Wplus_qgc));
  }
  if (Flavour(kf_h0).IsOn()) AddFlavour(fls,Flavour(kf_h0_qsc));
  if (Flavour(kf_gluon).IsOn()) AddFlavour(fls,Flavour(kf_gluon_qgc));
  return fls;
}

DECLARE_GETTER(SM_Getter,"SM",Model,std::string);

Model *SM_Getter::operator()(const std::string &key) const
{ 
  return new Model_SM();
}

void SM_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"standard model"; 
}

