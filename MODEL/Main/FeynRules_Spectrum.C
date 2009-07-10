#include "MODEL/Main/FeynRules_Spectrum.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Shell_Tools.H"

//additional stuff for sprintf etc
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;


FeynRules_Spectrum::FeynRules_Spectrum(Data_Reader * _dataread,
				       Model_Base * _model, string _dir) :
  Spectrum_Generator_Base(_dataread,_model),m_dir(_dir) { 

  m_identfile    = p_dataread->GetValue<string>("FR_IDENTFILE",string("ident_card.dat"));
  m_paramdeffile = p_dataread->GetValue<string>("FR_PARAMDEF",string("param_definition.dat"));
  m_paramfile    = p_dataread->GetValue<string>("FR_PARAMCARD",string("param_card.dat"));
  string intfile = p_dataread->GetValue<string>("FR_INTERACTIONS",string("interactions.dat"));
  rpa.gen.SetVariable("INTERACTION_DATA_FILE",intfile);
}

FeynRules_Spectrum::~FeynRules_Spectrum() 
{ 
  exh->RemoveTerminatorObject(this);
}

void FeynRules_Spectrum::PrepareTerminate()
{
  string path(rpa.gen.Variable("SHERPA_STATUS_PATH")+"/");
  if (path=="/") return;
  CopyFile(m_dir+"/"+m_identfile,path+m_identfile);
  CopyFile(m_dir+"/"+m_paramdeffile,path+m_paramdeffile);
  CopyFile(m_dir+"/"+m_paramfile,path+m_paramfile);
}

void FeynRules_Spectrum::Run() {

  msg_Info()<<"========== Generate FeynRules Model ========== "<<endl;
  
  SetUpTools();
  SetExternalParameters();
  FillAlgebra();
  SetInternalParameters();
  
  msg_Info()<<"================================================ "<<endl;

  delete p_reader;
  exh->AddTerminatorObject(this);
}

void FeynRules_Spectrum::SetUpTools() {
  
  p_reader = new ATOOLS::Data_Reader(" ",";","#","=");
  
  p_reader->AddWordSeparator("\t");
  p_reader->SetAddCommandLine(false);
  p_reader->SetIgnoreCase(true);
  p_reader->AddComment("!");
  p_reader->SetInputPath(m_dir);
  
  p_algebra = p_reader->Interpreter();

}

void FeynRules_Spectrum::PrepareReader(string filename) {

  p_reader->SetInputFile(filename);
  p_reader->AddFileEnd(string("Block"));
  
  p_reader->OpenInFile();
  
  if (p_reader->InFileMode()==3) 
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,
			    "The FeynRules input file "+m_dir+filename+" could not be read!",
			    "FeynRules_Spectrum","PrepareReader"));
}

void FeynRules_Spectrum::SetExternalParameters() {

  PrepareReader(m_identfile);
  
  msg_Info()<<endl<<" Reading file "<<m_identfile<<endl;

  vector<vector<string> > vdef;
  p_reader->RereadInFile();
  if (p_reader->MatrixFromFile(vdef,"")) {
  
    msg_Info()<<"   Fill external parameters : "<<endl;
    for (size_t i=0;i<vdef.size();++i) {
      if (vdef[i][0]!=string("DECAY")) {
	PrepareReader(m_paramfile);
	vector<vector<double> > vparam;
	p_reader->SetFileBegin(string("Block ")+vdef[i][0]);
	p_reader->RereadInFile();
	if (p_reader->MatrixFromFile(vparam,"")) {
	  if (vdef[i].size()==4) {
	    for (size_t j=0;j<vparam.size();++j) {
	      char str[12];
	      sprintf (str,"%i",int(vparam[j][0]));
	      if (str==vdef[i][1]) {
		msg_Info()<<"     "<<vdef[i][2]<<" = "<<vparam[j][1]<<endl;  
		if (vdef[i][2]!=string("aS")) {
		  p_model->GetScalarConstants()->insert(make_pair(vdef[i][2],vparam[j][1]));
		  break;
		}
		else {
		  p_model->GetScalarConstants()->insert(make_pair(string("aS(MZ)"),vparam[j][1]));
		  break;
		}
	      }
	    }
	  }
	  if (vdef[i].size()==5) {
	    for (size_t j=0;j<vparam.size();++j) {
	      char str_a[12];
	      sprintf (str_a,"%i",int(vparam[j][0]));
	      char str_b[12];
	      sprintf (str_b,"%i",int(vparam[j][1]));
	      if (str_a==vdef[i][1] && str_b==vdef[i][2] ) {
		msg_Info()<<"     "<<vdef[i][3]<<" = "<<vparam[j][2]<<endl;  
		p_model->GetScalarConstants()->insert(make_pair(vdef[i][3],vparam[j][2]));
		break;
	      }
	    }
	  }
	  if (vdef[i].size()>5) {
	    msg_Error()<<METHOD<<" can't interprete input ..."<<std::endl; 
	  }
	}
      }
    }
  }
  else msg_Error()<<"FeynRules input file seems to be empty ... "
		  <<m_dir+m_identfile<<endl;
  
  PrepareReader(m_paramfile);
  
  msg_Info()<<endl<<" Reading file "
		<<m_paramfile<<endl;
  
  //read-in masses
  msg_Info()<<endl<<"   Reading mass data: "<<endl;
  p_reader->SetFileBegin(string("Block MASS"));
  p_reader->AddFileEnd(string("DECAY"));
  p_reader->RereadInFile();
  vector<vector<double> >   vm;
  if(p_reader->MatrixFromFile(vm,"")) {
    Flavour flav;    
    for (size_t i=0;i<vm.size();++i) {
      flav.FromHepEvt(int(vm[i][0]));
      flav.SetMass(dabs(vm[i][1]));
      if (vm[i][1]<0) flav.SetMassSign(-1);
      msg_Info()<<"     Set mass of "<<flav<<" to "<<vm[i][1]<<endl;
    }
  }
  else msg_Error()<<endl<<"Could not read Block MASS from FeynRules file: "
		  <<m_dir+m_paramfile<<endl;


  //append SM parameters in Sherpa notation
  //Fermi constant
  double GF = p_model->ScalarConstant(string("Gf"));
  p_model->GetScalarConstants()->insert(make_pair(string("GF"),GF));
  //strong coupling
  double asMZ = p_model->ScalarConstant(string("aS(MZ)"));
  p_model->GetScalarConstants()->insert(make_pair(string("alpha_S(MZ)"),asMZ));
  int    order_alphaS	= p_dataread->GetValue<int>("ORDER_ALPHAS",0);
  double alphaS_default = p_dataread->GetValue<double>("ALPHAS(default)",asMZ);
  double MZ = Flavour(kf_Z).Mass();
  as = new Running_AlphaS(asMZ,sqr(MZ),order_alphaS);
  as->SetDefault(alphaS_default);
  p_model->GetScalarFunctions()->insert(make_pair(string("alpha_S"),as));
  //insert aS @ cpl scale
  double as_cpl = p_model->ScalarFunction(std::string("alpha_S"),rpa.gen.CplScale());
  p_model->GetScalarConstants()->insert(make_pair(string("aS"),as_cpl));
  //QED coupling and EW parameters
  p_model->GetScalarConstants()->insert(make_pair(string("alpha_QED(0)"),1./137.03599976));
  double aQEDMZ = 1./p_model->ScalarConstant(string("aEWM1"));
  aqed = new Running_AlphaQED(aQEDMZ,sqr(MZ));
  aqed->SetDefault(aQEDMZ);
  p_model->GetScalarFunctions()->insert(make_pair(string("alpha_QED"),aqed));
  double sin2TW = 1.-sqr(Flavour(kf_Wplus).Mass()/Flavour(kf_Z).Mass());
  p_model->GetScalarConstants()->insert(make_pair(std::string("sin2_thetaW"),sin2TW));
  
  msg_Info()<<"   set alphaS(MZ) to "
	   <<p_model->ScalarFunction(string("alpha_S"),sqr(MZ))<<" "<<sqr(MZ)<<endl;
  msg_Info()<<"   set alphaQED(MZ) to "
	   <<p_model->ScalarFunction(string("alpha_QED"),sqr(MZ))<<" "<<sqr(MZ)<<endl;
  
  //read-in decay widths
  msg_Info()<<endl<<"   Reading decay data: "<<endl;
  p_reader->SetFileEnd("dummy");
  vector<vector<string> >   vw;
  
  //
  p_reader->SetFileBegin(string("DECAY"));
  p_reader->RereadInFile();
  //
  p_reader->MatrixFromFile(vw,"");
  if (vw.size()>0) vw.front().insert(vw.front().begin(),"DECAY");
  
  for (size_t k=0;k<vw.size();++k) {
    if (vw[k].size()==3) {
      if (vw[k][0]=="DECAY") {
       Flavour flav;    
       flav.FromHepEvt(ATOOLS::ToType<int>(vw[k][1]));
       flav.SetWidth(ATOOLS::ToType<double>(vw[k][2]));
       msg_Info()<<"     Set width of "<<flav<<" to "<<flav.Width()<<endl;
      }
    }
  }  
}

void FeynRules_Spectrum::FillAlgebra() {
 
  //fill external parameter into algebra
  msg_Info()<<"\n   Add algebra relation : "<<endl;
     map<string,double>::iterator mit = p_model->GetScalarConstants()->begin();
  for (;mit!=p_model->GetScalarConstants()->end();++mit) {
    msg_Info()<<"     "<<mit->first<<" = "<<mit->second<<" vs. "<<ToString(mit->second)<<endl; 
    p_algebra->AddTag(mit->first,ToString(mit->second));
  }
}


void FeynRules_Spectrum::SetInternalParameters() {

  PrepareReader(m_paramdeffile);
  p_reader->SetFileBegin("!");
  p_reader->AddComment("!");
  
  msg_Info()<<endl<<" Reading file "
		<<m_paramdeffile<<"\n"<<endl;
  
  vector<vector<string> > vdef;
  p_reader->RereadInFile();
  msg_Info()<<"   Derived parameters : "<<endl;
  if (p_reader->MatrixFromFile(vdef,"")) {
    for (size_t i=0;i<vdef.size();++i) {
      //add derived parameters to ScalarConstants & Algebra
      if (vdef[i][2]=="R") {
	double val = ToType<double>(p_algebra->Interprete(vdef[i][1]));
	//std::cout<<" interprete R : "<<vdef[i][1]<<" -> "<<val<<std::endl; 
	p_algebra->AddTag(vdef[i][0],ToString(val));
	p_model->GetScalarConstants()->insert(make_pair(vdef[i][0],val));
	msg_Info()<<"     "<<vdef[i][0]<<" = "<<val<<endl;
      }
      if (vdef[i][2]=="C") {
	Complex val = ToType<Complex>(p_algebra->Interprete(vdef[i][1]));
	//std::cout<<" interprete C : "<<vdef[i][1]<<" -> "<<val<<std::endl; 
	p_algebra->AddTag(vdef[i][0],ToString(val));
	p_model->GetComplexConstants()->insert(make_pair(vdef[i][0],val));
	msg_Info()<<"     "<<vdef[i][0]<<" = "<<val<<endl;
      }
    }
  }
  else {
    msg_Error()<<" Error reading "<<m_paramdeffile<<
      " in FeynRules_Spectrum::SetInternalParameters() "<<endl;  
  }
}
