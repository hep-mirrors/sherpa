#include "ADD.H"
#include "Message.H"
#include "Standard_Model.H"
#ifdef USING__ISAJET
#include "Isajet_Fortran_Interface.H"
#else
#include "Spectrum_Generator_Base.H"
#endif

using namespace MODEL;
using namespace ATOOLS;


ADD::ADD(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the ADD from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("ADD");

  p_sm = new Standard_Model(m_dir,m_file);
  p_numbers   = new ScalarNumbersMap(*(p_sm->GetScalarNumbers()));
  p_constants = new ScalarConstantsMap(*(p_sm->GetScalarConstants()));
  p_functions = new ScalarFunctionsMap(*(p_sm->GetScalarFunctions()));
  p_matrices  = new ComplexMatricesMap(*(p_sm->GetComplexMatrices()));

  ReadInFile();
}

ADD::~ADD()
{
  delete p_sm;
}

void ADD::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);

  
  p_constants->insert(std::make_pair(std::string("G_Newton"), 
				     p_dataread->GetValue<double>("G_NEWTON",6.707E-39)));
  p_numbers->insert(std::make_pair(std::string("ED"), 
				     p_dataread->GetValue<int>("N_ED",2)));
  p_constants->insert(std::make_pair(std::string("M_s"), 
				     p_dataread->GetValue<double>("M_S",0.)));
  p_constants->insert(std::make_pair(std::string("M_cut"), 
				     p_dataread->GetValue<double>("M_CUT",ScalarConstant(std::string("M_s")))));
  p_constants->insert(std::make_pair(std::string("Radius"), 
				     p_dataread->GetValue<double>("RADIUS",0.)));
  p_numbers->insert(std::make_pair(std::string("KK_mode"), 
				     p_dataread->GetValue<int>("KK_CONVENTION",1)));


  int    mode = ScalarNumber(std::string("KK_mode"));
  double rad = ScalarConstant(std::string("Radius"));
  int    ed  = ScalarNumber(std::string("ED"));
  double gn  = ScalarConstant(std::string("G_Newton"));
  double ms  = ScalarConstant(std::string("M_s"));

  switch(mode){
  case 1:case 2:                            //HLZ
    //Calculation of Gamma(ed/2)
    double gam;
    if(ed%2==0) gam=1.;
    else gam=sqrt(M_PI);
    for(int i=2-ed%2;i<ed;i+=2)gam*=0.5*i;
    
    //If Radius is set but not the scale M_s, M_s is calculated
    if (IsZero(ms) && rad > 0.) {
      ms=pow(gam*pow(4.*M_PI,.5*ed)/pow(2.*M_PI*rad,1.*ed)/gn,1./(2.+ed));
      (*p_constants)[std::string("M_s")] = ms;
    }
    
    (*p_constants)[std::string("Radius")] = 
      pow(gam*pow(4.*M_PI,.5*ed)/pow(ms,2.+(double(ed)))/gn/2.,1./(double(ed)))/2/M_PI;
    break;
  case 5:                                   //GRW
    //If Radius is set but not the scale M_s, M_s is calculated
    if (IsZero(ms) && rad > 0.) {
      ms=pow(8.*M_PI*pow(rad,1.*ed)*gn,-1./(2.+ed));
      (*p_constants)[std::string("M_s")] = ms;
    }
    
    (*p_constants)[std::string("Radius")] = pow(8.*M_PI*pow(ms,2.+(double(ed)))*gn,-1./(double(ed)));
  }

  p_constants->insert(std::make_pair(std::string("kappa"),sqrt(8.*M_PI*gn))); 
  p_constants->insert(std::make_pair(std::string("omega"),sqrt(4.*(-1.+ed)/(3.*(2.+ed)))));
  p_constants->insert(std::make_pair(std::string("M2_s"),sqr(ms)));
}


