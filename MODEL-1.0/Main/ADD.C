#include "ADD.H"
#include "Message.H"
#include "Standard_Model.H"
#include "Isajet_Fortran_Interface.H"

using namespace MODEL;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;

ADD::ADD(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg.Events()<<"Initialize the ADD from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("ADD");

  Model_Base * SM = new Standard_Model(m_dir,m_file);
  p_numbers   = new ScalarNumbersMap(*(SM->GetScalarNumbers()));
  p_constants = new ScalarConstantsMap(*(SM->GetScalarConstants()));
  p_functions = new ScalarFunctionsMap(*(SM->GetScalarFunctions()));
  p_matrices  = new ComplexMatricesMap(*(SM->GetComplexMatrices()));

  delete SM;

  ReadInFile();
}

void ADD::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);

  
  p_constants->insert(std::make_pair(std::string("G_Newton"), 
				     p_dataread->GetValue<double>("G_NEWTON",6.707E-39)));
  p_numbers->insert(std::make_pair(std::string("ED"), 
				     p_dataread->GetValue<int>("N_ED",2)));
  p_constants->insert(std::make_pair(std::string("M_s"), 
				     p_dataread->GetValue<double>("M_S",0.)));
  p_constants->insert(std::make_pair(std::string("Radius"), 
				     p_dataread->GetValue<double>("RADIUS",0.)));
  p_numbers->insert(std::make_pair(std::string("KK_mode"), 
				     p_dataread->GetValue<int>("KK_CONVENTION",1)));
  

  double rad = ScalarConstant(std::string("Radius"));
  int    ed  = ScalarNumber(std::string("ED"));
  double gn  = ScalarConstant(std::string("G_Newton"));
  double ms  = ScalarConstant(std::string("M_s"));

  //Calculation of Gamma(ed/2)
  double gam;
  if(ed%2==0) gam=1.;
  else gam=sqrt(M_PI);
  for(int i=2-ed%2;i<ed;i+=2)gam*=0.5*i;
  
  //If Radius is set but not the scale M_s, M_s is calculated
  if (IsZero(ms) && rad > 0.) {
    ms=pow(gam*pow(4.*M_PI,.5*ed)/pow(rad,1.*ed)/gn,1./(2.+ed));
    (*p_constants)[std::string("M_s")] = ms;
  }

  (*p_constants)[std::string("Radius")] = pow(gam*pow(4.*M_PI,.5*ed)/pow(ms,2.+(double(ed)))/gn,1./(double(ed)));
  p_constants->insert(std::make_pair(std::string("kappa"),sqrt(8.*M_PI*gn))); 
  p_constants->insert(std::make_pair(std::string("omega"),sqrt(4.*(-1.+ed)/(3.*(2.+ed)))));
  p_constants->insert(std::make_pair(std::string("M2_s"),sqr(ms)));
}


