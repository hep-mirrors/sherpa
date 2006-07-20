#include "MSSM.H"
#include "Message.H"
#include "Standard_Model.H"
#include "LesHouches_Interface.H"
#ifdef USING__ISAJET
#include "Isajet_Fortran_Interface.H"
#else
#include "Spectrum_Generator_Base.H"
#endif

using namespace MODEL;
using namespace ATOOLS;


MSSM::MSSM(std::string _dir,std::string _file) :
  Model_Base(_dir,_file)
{
  msg_Info()<<"Initialize the MSSM from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("MSSM");

  Model_Base * SM = new Standard_Model(m_dir,m_file);
  p_numbers   = new ScalarNumbersMap(*(SM->GetScalarNumbers()));
  p_constants = new ScalarConstantsMap(*(SM->GetScalarConstants()));
  p_functions = new ScalarFunctionsMap(*(SM->GetScalarFunctions()));
  p_matrices  = new ComplexMatricesMap(*(SM->GetComplexMatrices()));

  p_constants->insert(std::make_pair(std::string("mT"),    
				     SM->ScalarConstant("Yukawa_t")));

  delete SM;

  ReadInFile();
}

void MSSM::ReadInFile() {
  p_dataread = new Data_Read(m_dir+m_file);

  m_scales     = m_unification = 0;
  m_benchmark  = p_dataread->GetValue<std::string>("BENCHMARK",std::string(""));
  m_spectrum   = p_dataread->GetValue<int>("GENERATOR_ON",1);

  if (m_benchmark.length()!=0) InitializeBenchmarkPoint();
  else {
    m_scenario = p_dataread->GetValue<std::string>("SUSY_SCENARIO",std::string("mSUGRA"));

    if (m_scenario==std::string("AMSB")) {
      p_constants->insert(std::make_pair(std::string("m0"),    
					 p_dataread->GetValue<double>("M0",0.)));
      p_constants->insert(std::make_pair(std::string("m32"),    
					 p_dataread->GetValue<double>("M32",0.)));
      p_constants->insert(std::make_pair(std::string("tan(beta)"),    
					 p_dataread->GetValue<double>("TAN(BETA)",0.)));
      p_constants->insert(std::make_pair(std::string("A0"),    
					 p_dataread->GetValue<double>("A0",0.)));
      p_numbers->insert(std::make_pair(std::string("sign(mu)"),    
				       p_dataread->GetValue<int>("SIGN(MU)",1)));
    }
    else if (m_scenario==std::string("mGMSB")) {
      p_constants->insert(std::make_pair(std::string("m_mes"),    
					 p_dataread->GetValue<double>("MESSENGER_SCALE",2000.)));
      p_constants->insert(std::make_pair(std::string("Lambda_m"),    
					 p_dataread->GetValue<double>("LAMBDA_M",0.)));
      p_constants->insert(std::make_pair(std::string("tan(beta)"),    
					 p_dataread->GetValue<double>("TAN(BETA)",0.)));
      p_constants->insert(std::make_pair(std::string("c_grav"),    
					 p_dataread->GetValue<double>("GRAVITINO_RATIO",1.)));
      p_numbers->insert(std::make_pair(std::string("sign(mu)"),    
				       p_dataread->GetValue<int>("SIGN(MU)",1)));
      p_numbers->insert(std::make_pair(std::string("n_mes"),    
				       p_dataread->GetValue<int>("NUMBER_MESSENGERS",1)));
    }
    else if (m_scenario==std::string("non-minimal GMSB")) {
      p_constants->insert(std::make_pair(std::string("m_mes"),    
					 p_dataread->GetValue<double>("MESSENGER_SCALE",2000.)));
      p_constants->insert(std::make_pair(std::string("Lambda_m"),    
					 p_dataread->GetValue<double>("LAMBDA_M",0.)));
      p_constants->insert(std::make_pair(std::string("tan(beta)"),    
					 p_dataread->GetValue<double>("TAN(BETA)",0.)));
      p_constants->insert(std::make_pair(std::string("c_grav"),    
					 p_dataread->GetValue<double>("GRAVITINO_RATIO",1.)));
      p_constants->insert(std::make_pair(std::string("c_gauge"),    
					 p_dataread->GetValue<double>("GAUGINO_RATIO",1.)));
      p_constants->insert(std::make_pair(std::string("Delta_Hu"),    
					 p_dataread->GetValue<double>("DELTA_M(Hu)",1.)));
      p_constants->insert(std::make_pair(std::string("Delta_Hd"),    
					 p_dataread->GetValue<double>("DELTA_M(Hd)",1.)));
      p_constants->insert(std::make_pair(std::string("Delta_Y"),    
					 p_dataread->GetValue<double>("DELTA_M(Y)",1.)));
      p_constants->insert(std::make_pair(std::string("n5_1"),    
					 p_dataread->GetValue<double>("N5(U1)",1.)));
      p_constants->insert(std::make_pair(std::string("n5_2"),    
					 p_dataread->GetValue<double>("N5(SU2)",1.)));
      p_constants->insert(std::make_pair(std::string("n5_3"),    
					 p_dataread->GetValue<double>("N5(SU3)",1.)));
      p_numbers->insert(std::make_pair(std::string("sign(mu)"),    
				       p_dataread->GetValue<int>("SIGN(MU)",1)));
      p_numbers->insert(std::make_pair(std::string("n_mes"),    
				       p_dataread->GetValue<int>("NUMBER_MESSENGERS",1)));
    }
    else if (m_scenario==std::string("non-universal SUGRA")) {
      p_constants->insert(std::make_pair(std::string("m0"),    
					 p_dataread->GetValue<double>("M0",0.)));
      p_constants->insert(std::make_pair(std::string("m12"),    
					 p_dataread->GetValue<double>("M12",0.)));
      p_constants->insert(std::make_pair(std::string("tan(beta)"),    
					 p_dataread->GetValue<double>("TAN(BETA)",0.)));
      p_constants->insert(std::make_pair(std::string("A0"),    
					 p_dataread->GetValue<double>("A0",0.)));
      p_constants->insert(std::make_pair(std::string("m_U1"),    
					 p_dataread->GetValue<double>("M(U1)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_SU2"),    
					 p_dataread->GetValue<double>("M(SU2)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_SU3"),    
					 p_dataread->GetValue<double>("M(SU3)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_e1R"),    
					 p_dataread->GetValue<double>("M(e1R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_l1"),    
					 p_dataread->GetValue<double>("M(l1)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_u1R"),    
					 p_dataread->GetValue<double>("M(u1R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_d1R"),    
					 p_dataread->GetValue<double>("M(d1R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_q1"),    
					 p_dataread->GetValue<double>("M(q1)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_e2R"),    
					 p_dataread->GetValue<double>("M(e2R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_l2"),    
					 p_dataread->GetValue<double>("M(l2)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_u2R"),    
					 p_dataread->GetValue<double>("M(u2R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_d2R"),    
					 p_dataread->GetValue<double>("M(d2R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_q2"),    
					 p_dataread->GetValue<double>("M(q2)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_e3R"),    
					 p_dataread->GetValue<double>("M(e3R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_l3"),    
					 p_dataread->GetValue<double>("M(l3)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_u3R"),    
					 p_dataread->GetValue<double>("M(u3R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_d3R"),    
					 p_dataread->GetValue<double>("M(d3R)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_q3"),    
					 p_dataread->GetValue<double>("M(q3)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("a_tau"),    
					 p_dataread->GetValue<double>("A(tau)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("a_b"),    
					 p_dataread->GetValue<double>("A(b)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("a_t"),    
					 p_dataread->GetValue<double>("A(t)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("H_d"),    
					 p_dataread->GetValue<double>("H(d)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("H_u"),    
					 p_dataread->GetValue<double>("H(u)_GUT",1.e20)));
      p_constants->insert(std::make_pair(std::string("m_SUSY"),    
					 p_dataread->GetValue<double>("SCALE(SUSY)",0.)));
      p_numbers->insert(std::make_pair(std::string("sign(mu)"),    
				       p_dataread->GetValue<int>("SIGN(MU)",1)));
      //if (m_SUSYScale>0.) m_scales = 1;
    }
    else if (m_scenario==std::string("SUGRA with enforced unification")) {
      m_unification = 1;
      p_constants->insert(std::make_pair(std::string("m0"),    
					 p_dataread->GetValue<double>("M0",0.)));
      p_constants->insert(std::make_pair(std::string("m12"),    
					 p_dataread->GetValue<double>("M12",0.)));
      p_constants->insert(std::make_pair(std::string("tan(beta)"),    
					 p_dataread->GetValue<double>("TAN(BETA)",0.)));
      p_constants->insert(std::make_pair(std::string("A0"),    
					 p_dataread->GetValue<double>("A0",0.)));
      p_numbers->insert(std::make_pair(std::string("sign(mu)"),    
				       p_dataread->GetValue<int>("SIGN(MU)",1)));
    }
    else if (m_scenario==std::string("mSUGRA")) {
      p_constants->insert(std::make_pair(std::string("m0"),    
					 p_dataread->GetValue<double>("M0",0.)));
      p_constants->insert(std::make_pair(std::string("m12"),    
					 p_dataread->GetValue<double>("M12",0.)));
      p_constants->insert(std::make_pair(std::string("tan(beta)"),    
					 p_dataread->GetValue<double>("TAN(BETA)",0.)));
      p_constants->insert(std::make_pair(std::string("A0"),    
					 p_dataread->GetValue<double>("A0",0.)));
      p_numbers->insert(std::make_pair(std::string("sign(mu)"),    
				       p_dataread->GetValue<int>("SIGN(MU)",1)));
    }
  }
}


bool MSSM::RunSpectrumGenerator() {
  if (m_spectrum) {
    m_generator = p_dataread->GetValue<std::string>("SUSY_GENERATOR",std::string("LesHouches"));
#ifdef USING__ISAJET
    if (m_generator==std::string("Isajet")) {
      p_spectrumgenerator = new ISAJET::Isajet_Fortran_Interface(p_dataread,this);
      p_spectrumgenerator->Run(std::string(m_scenario));
      p_spectrumgenerator->FillMasses();
      //p_spectrumgenerator->FillDecays();
      return 1;
    }
#endif
    if (m_generator==std::string("LesHouches")) {
      p_spectrumgenerator = new LesHouches_Interface(p_dataread,this,m_dir);
      p_spectrumgenerator->Run(std::string(m_scenario));
      p_spectrumgenerator->FillMasses();
      //p_spectrumgenerator->FillDecays();
      return 1;
    }
    
    msg.Error()<<"Error in MSSM::RunSpectrumGenerator."<<std::endl
	       <<"   Unknown spectrum generator : "<<m_generator<<" use internal solution."<<std::endl;
    return 0;
  }
  return 1;
}
