#include "MODEL/Interaction_Models/Interaction_Model_FeynRules.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include <stdio.h>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_FeynRules_Getter,"FeynRules",
	       Interaction_Model_Base,Interaction_Model_Arguments);


Interaction_Model_Base *Interaction_Model_FeynRules_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_FeynRules
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_FeynRules_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The FeynRules Model"; 
}


Interaction_Model_FeynRules::Interaction_Model_FeynRules(MODEL::Model_Base * _model,
							 std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("FeynRules",_model,_cplscheme,_yukscheme)
{
  m_interactionfile = rpa.gen.Variable("INTERACTION_DATA_FILE");
  
  //setup the data reader
  p_reader = new ATOOLS::Data_Reader(" ",";","#","=");
  p_reader->AddWordSeparator("\t");
  p_reader->SetAddCommandLine(false);
  p_reader->SetIgnoreCase(true);
  p_reader->AddComment("!");
  //p_reader->SetInputPath(_model->m_dir);
  p_reader->SetInputFile(m_interactionfile);
  
  //equip the algebra interpreter
  p_algebra = p_reader->Interpreter();
  msg_Tracking()<<"\n   Add algebra relation : "<<endl;
  //scalar constants
  map<string,double>::iterator mit = _model->GetScalarConstants()->begin();
  for (;mit!=_model->GetScalarConstants()->end();++mit) {
    msg_Tracking()<<"     "<<mit->first<<" = "<<mit->second<<" vs. "<<ToString(mit->second)<<endl; 
    p_algebra->AddTag(mit->first,ToString(mit->second));
  }
  //complex constants
  map<string,Complex>::iterator cmit = _model->GetComplexConstants()->begin();
  for (;cmit!=_model->GetComplexConstants()->end();++cmit) {
    msg_Tracking()<<"     "<<cmit->first<<" = "<<cmit->second<<" vs. "<<ToString(cmit->second)<<endl; 
    p_algebra->AddTag(cmit->first,ToString(cmit->second));
  }

}

//three-point interactions
void Interaction_Model_FeynRules::c_FFV(std::vector<Single_Vertex>& vertex,int & vanz) {
  //
  p_reader->SetFileBegin(string("!"));
  p_reader->RereadInFile();
  //
  vector<vector<string> > vv;
  p_reader->MatrixFromFile(vv,"");
  bool valid = true;
  for (size_t i=0;i<vv.size();++i) {
    if (vv[i].size()>0) {
      if (vv[i][0]=="VERTEX") {
	valid = true;
	if (vv[i].size()>4) {
	  valid=false; 
	  continue;
	}
	vertex[vanz].nleg=3;
	for (size_t j=0;j<3;j++) {
	  int kf = ATOOLS::ToType<int>(vv[i][j+1]);
	  Flavour flav = kf>0?Flavour(kf):Flavour(-kf).Bar();
	  vertex[vanz].in[j] = flav;
	  if (!flav.IsOn()) valid=false;
	}
      }
      if (vv[i][0]=="1" && valid) {
	Complex cpl(0.,0.);
	//std::cout<<" interprete RC: "<<vv[i][1]<<std::endl;
	cpl = ToType<Complex>(p_algebra->Interprete(vv[i][1]));
	//std::cout<<" interpreted RC: "<<vv[i][1]<<" -> "<<cpl<<std::endl; 
	Kabbala kcpl = Kabbala(vv[i][1],cpl);
	vertex[vanz].cpl[0]  = kcpl;
	vertex[vanz].Str    += kcpl.String()+"*P_R";
      }
      if (vv[i][0]=="2" && valid) {
	Complex cpl(0.,0.);
	//std::cout<<" interprete LC: "<<vv[i][1]<<std::endl;
	cpl = ToType<Complex>(p_algebra->Interprete(vv[i][1]));
	//std::cout<<" interpreted LC: "<<vv[i][1]<<" -> "<<cpl<<std::endl; 
	Kabbala kcpl = Kabbala(vv[i][1],cpl);
	vertex[vanz].cpl[1]  = kcpl;
	vertex[vanz].Str    += kcpl.String()+"*P_L";
      }
      if (vv[i][0]=="3" && valid) {
	switch (vv[i][1][0]) {
	case 'N' : {
	  vertex[vanz].Color.push_back(Color_Function(cf::None));
	  break;
	}
	case 'F' : {
	  int arg1   = int(vv[i][1][2])-49;
	  int arg2   = int(vv[i][1][4])-49;
	  int arg3   = int(vv[i][1][6])-49;
	  char sarg1[2],sarg2[2],sarg3[2];
	  sprintf(sarg1,"%i",arg1);
	  sprintf(sarg2,"%i",arg2);
	  sprintf(sarg3,"%i",arg3);
	  //
	  vertex[vanz].Color.push_back(Color_Function(cf::F,arg1,arg2,arg3,sarg1[0],sarg2[0],sarg3[0]));
	  break;
	}
	case 'T' : {
	  int arg1   = int(vv[i][1][2])-49;
	  int arg2   = int(vv[i][1][4])-49;
	  int arg3   = int(vv[i][1][6])-49;
	  char sarg1[2],sarg2[2],sarg3[2];
	  sprintf(sarg1,"%i",arg1);
	  sprintf(sarg2,"%i",arg2);
	  sprintf(sarg3,"%i",arg3);
	  //
	  vertex[vanz].Color.push_back(Color_Function(cf::T,arg1,arg2,arg3,sarg1[0],sarg2[0],sarg3[0]));
	  break;
	}
	case 'D' : {
	  int arg1   = int(vv[i][1][2])-49;
	  int arg2   = int(vv[i][1][4])-49;
	  char sarg1[2],sarg2[2];
	  sprintf(sarg1,"%i",arg1);
	  sprintf(sarg2,"%i",arg2);
	  //
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(arg1,arg2);     
	  vertex[vanz].Color.back().SetStringArg(sarg1[0],sarg2[0]);
	  break;
	}
	case 'G' : {
	  int arg1   = int(vv[i][1][2])-49;
	  int arg2   = int(vv[i][1][4])-49;
	  char sarg1[2],sarg2[2];
	  sprintf(sarg1,"%i",arg1);
	  sprintf(sarg2,"%i",arg2);
	  //
	  vertex[vanz].Color.push_back(Color_Function(cf::G));     
	  vertex[vanz].Color.back().SetParticleArg(arg1,arg2);     
	  vertex[vanz].Color.back().SetStringArg(sarg1[0],sarg2[0]);
	  break;
	}
	default : {
	  valid = false;
	  msg_Error()<<" Can't handle colour structure "<<vv[i][1]<<endl;  	
	}
	}
      }
      if (vv[i][0]=="4"  && valid) {
	if (vv[i][1]=="FFS") {
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	}
	else if (vv[i][1]=="FFV") {
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	  //correct cpl sign
	  Kabbala kcpl0 = Kabbala("(-1.)*("+vertex[vanz].cpl[0].String()+")",vertex[vanz].cpl[0].Value()*(-1.));
	  Kabbala kcpl1 = Kabbala("(-1.)*("+vertex[vanz].cpl[1].String()+")",vertex[vanz].cpl[1].Value()*(-1.));
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	}
	else if (vv[i][1]=="Gauge3") {
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge3",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2);     
	}
	else if(vv[i][1]=="SSS") {
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	}
	else if(vv[i][1]=="VVS") {
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
	}
	else if(vv[i][1]=="SSV") {
	  //correct cpl sign
	  Kabbala kcpl0 = Kabbala("(-1.)*("+vertex[vanz].cpl[0].String()+")",vertex[vanz].cpl[0].Value()*(-1.));
	  Kabbala kcpl1 = Kabbala("(-1.)*("+vertex[vanz].cpl[1].String()+")",vertex[vanz].cpl[1].Value()*(-1.));
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);
	}
	else {
	  msg_Error()<<" Lorentz structure missing: 3-particle interaction "<<vv[i][1]<<endl; 
	  valid = false;
	}
	//fill in actual vertex
	if (valid && vertex[vanz].nleg!=4) {
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());
	  vanz++;
	  valid = true;
	}
      }
    }
  }
}


void Interaction_Model_FeynRules::c_FFS(std::vector<Single_Vertex>& vertex,int & vanz) {}
void Interaction_Model_FeynRules::c_VVV(std::vector<Single_Vertex>& vertex,int & vanz) {}
void Interaction_Model_FeynRules::c_VVS(std::vector<Single_Vertex>& vertex,int & vanz) {}
void Interaction_Model_FeynRules::c_SSV(std::vector<Single_Vertex>& vertex,int & vanz) {}
void Interaction_Model_FeynRules::c_SSS(std::vector<Single_Vertex>& vertex,int & vanz) {}

//four-point interactions
void Interaction_Model_FeynRules::c_VVVV(std::vector<Single_Vertex>& vertex,int & vanz){

  //
  p_reader->SetFileBegin(string("!"));
  p_reader->RereadInFile();
  //
  vector<vector<string> > vv;
  p_reader->MatrixFromFile(vv,"");
  bool valid = true;
  for (size_t i=0;i<vv.size();++i) {
    if (vv[i].size()>0) {
      if (vv[i][0]=="VERTEX") {
	valid = true;
	if (vv[i].size()==5) {
	  vertex[vanz].nleg=4;
	  for (size_t j=0;j<4;j++) {
	    int kf = ATOOLS::ToType<int>(vv[i][j+1]);
	    Flavour flav = kf>0?Flavour(kf):Flavour(-kf).Bar();
	    vertex[vanz].in[j] = flav;
	    if (!flav.IsOn()) valid=false;
	  }
	}
	else {
	  valid = false;
	  if (vv[i].size()>5) {
	    msg_Error()<<" Interaction_Model_FeynRules :";
	    msg_Error()<<"Five- and more point interactions not supported! "<<std::endl; 
	  }
	  continue;
	}
      }
      if (vv[i][0]=="1" && valid) {
	Complex cpl(0.,0.);
	cpl = ToType<Complex>(p_algebra->Interprete(vv[i][1]));
	//std::cout<<" interpreted RC: "<<vv[i][1]<<" -> "<<cpl<<std::endl; 
	Kabbala kcpl = Kabbala(vv[i][1],cpl);
	vertex[vanz].cpl[0]  = kcpl;
	vertex[vanz].Str    += kcpl.String()+"*P_R";
      }
      if (vv[i][0]=="2" && valid) {
	Complex cpl(0.,0.);
	cpl = ToType<Complex>(p_algebra->Interprete(vv[i][1]));
	//std::cout<<" interpreted LC: "<<vv[i][1]<<" -> "<<cpl<<std::endl; 
	Kabbala kcpl = Kabbala(vv[i][1],cpl);
	vertex[vanz].cpl[1]  = kcpl;
	vertex[vanz].Str    += kcpl.String()+"*P_L";
      }
      if (vv[i][0]=="3" && valid) {
	switch (vv[i][1][0]) {
	case 'N' : {
	  vertex[vanz].Color.push_back(Color_Function(cf::None));
	  break;
	}
	case 'F' : {
	  //needs to be generalised
	  if (vv[i][1].size()==8) {
	    int arg1   = int(vv[i][1][2])-49;
	    int arg2   = int(vv[i][1][4])-49;
	    int arg3   = int(vv[i][1][6])-49;
	    char sarg1[2],sarg2[2],sarg3[2];
	    sprintf(sarg1,"%i",arg1);
	    sprintf(sarg2,"%i",arg2);
	    sprintf(sarg3,"%i",arg3);
	    //
	    vertex[vanz].Color.push_back(Color_Function(cf::F,arg1,arg2,arg3,sarg1[0],sarg2[0],sarg3[0]));
	  }
	  else vertex[vanz].Color.push_back(Color_Function(cf::F,0,1,2,'0','1','2'));
	}
	case 'T' : {
	  if (vv[i][1].size()==8) {
	    int arg1   = int(vv[i][1][2])-49;
	    int arg2   = int(vv[i][1][4])-49;
	    int arg3   = int(vv[i][1][6])-49;
	    char sarg1[2],sarg2[2],sarg3[2];
	    sprintf(sarg1,"%i",arg1);
	    sprintf(sarg2,"%i",arg2);
	    sprintf(sarg3,"%i",arg3);
	    //
	    vertex[vanz].Color.push_back(Color_Function(cf::T,arg1,arg2,arg3,sarg1[0],sarg2[0],sarg3[0]));
	  }
	  if (vv[i][1].size()==10) {
	    int arg1   = int(vv[i][1][2])-49;
	    int arg2   = int(vv[i][1][4])-49;
	    int arg3   = int(vv[i][1][6])-49;
	    int arg4   = int(vv[i][1][8])-49;
	    char sarg1[2],sarg2[2],sarg3[2],sarg4[2];
	    sprintf(sarg1,"%i",arg1);
	    sprintf(sarg2,"%i",arg2);
	    sprintf(sarg3,"%i",arg3);
	    sprintf(sarg4,"%i",arg4);
	    //
	    vertex[vanz].Color.clear();
	    vertex[vanz].Color.resize(2);
	    vertex[vanz].Lorentz.resize(2); 
	    //
	    vertex[vanz].Color[0] = Color_Function(cf::T,arg1,arg3,4,sarg1[0],sarg3[0],'4',
						   new Color_Function(cf::T,arg2,4,arg4,sarg2[0],'4',sarg4[0]));
	    vertex[vanz].Lorentz[0]= LF_Getter::GetObject("VVSS",LF_Key());     
	    vertex[vanz].Lorentz[0]->SetParticleArg(0,3);     
	    //
	    vertex[vanz].Color[1] = Color_Function(cf::T,arg2,arg3,4,sarg2[0],sarg3[0],'4',
						   new Color_Function(cf::T,arg1,4,arg4,sarg1[0],'4',sarg4[0]));
	    vertex[vanz].Lorentz[1]= LF_Getter::GetObject("VVSS",LF_Key());     
	    vertex[vanz].Lorentz[1]->SetParticleArg(0,3);     
	  }
	  break;
	}
	case 'D' : {
	  int arg1   = int(vv[i][1][2])-49;
	  int arg2   = int(vv[i][1][4])-49;
	  char sarg1[2],sarg2[2];
	  sprintf(sarg1,"%i",arg1);
	  sprintf(sarg2,"%i",arg2);
	  //
	  vertex[vanz].Color.push_back(Color_Function(cf::D,arg1,arg2,-1,sarg1[0],sarg2[0]));
	  break;
	}
	case 'G' : {
	  int arg1   = int(vv[i][1][2])-49;
	  int arg2   = int(vv[i][1][4])-49;
	  char sarg1[2],sarg2[2];
	  sprintf(sarg1,"%i",arg1);
	  sprintf(sarg2,"%i",arg2);
	  //
	  vertex[vanz].Color.push_back(Color_Function(cf::G,arg1,arg2,-1,sarg1[0],sarg2[0]));
	  break;
	}
	default : {
	  valid = false;
	  msg_Error()<<" Can't handle colour structure "<<vv[i][1]<<endl;  	
	}
	}
      }
      if (vv[i][0]=="4" && valid) {
	if (vv[i][1]=="Gauge4") {
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge4",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2,3);     
	  
	  //correct coupling
	  Kabbala kcpl0 = Kabbala("i*("+vertex[vanz].cpl[0].String()+")",vertex[vanz].cpl[0].Value()*Complex(0.,1.));
	  Kabbala kcpl1 = Kabbala("i*("+vertex[vanz].cpl[1].String()+")",vertex[vanz].cpl[1].Value()*Complex(0.,1.));
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;
	}
      	else if (vv[i][1]=="Gluon4") {  
	  vertex[vanz].Color.clear();
	  vertex[vanz].Color.resize(3);
	  vertex[vanz].Lorentz.resize(3); 
	  vertex[vanz].Color[0]   = Color_Function(cf::F,0,2,4,'0','2','4',
							new Color_Function(cf::F,1,3,4,'1','3','4'));
	  vertex[vanz].Lorentz[0] = LF_Getter::GetObject("Gluon4",LF_Key());
	  vertex[vanz].Lorentz[0]->SetParticleArg(0,1,2,3);     
	  
	  vertex[vanz].Color[1]        = Color_Function(cf::F,0,3,4,'0','3','4',
							new Color_Function(cf::F,1,2,4,'1','2','4'));
	  vertex[vanz].Lorentz[1] = LF_Getter::GetObject("Gluon4",LF_Key());
	  vertex[vanz].Lorentz[1]->SetParticleArg(0,1,3,2);     
	  
	  vertex[vanz].Color[2]   = Color_Function(cf::F,0,1,4,'0','1','4',
							new Color_Function(cf::F,3,2,4,'3','2','4')); 
	  vertex[vanz].Lorentz[2] = LF_Getter::GetObject("Gluon4",LF_Key());     
	  vertex[vanz].Lorentz[2]->SetParticleArg(0,3,1,2);     
	  //correct coupling
	  Kabbala kcpl0 = Kabbala(vertex[vanz].cpl[0].String()+"*i",vertex[vanz].cpl[0].Value()*Complex(0.,1.));
	  Kabbala kcpl1 = Kabbala(vertex[vanz].cpl[1].String()+"*i",vertex[vanz].cpl[1].Value()*Complex(0.,1.));
	  vertex[vanz].cpl[0] = kcpl0;
	  vertex[vanz].cpl[1] = kcpl1;

  	}
	else if(vv[i][1]=="SSSS") {
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));
	}
	else if(vv[i][1]=="VVSS") {
	  if (vertex[vanz].Color.size()==1) {
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
	    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
	  }
	}
	else {
	  msg_Error()<<" Lorentz structure missing: 4-particle interaction "<<vv[i][1]<<endl; 
	  valid = false;
	}
	//fill in actual vertex
	if (valid && vertex[vanz].nleg==4) {
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());
	  vanz++;
	  valid = true;
	}
      }
    }
  }
}

void Interaction_Model_FeynRules::c_SSVV(std::vector<Single_Vertex>& vertex,int & vanz){}
void Interaction_Model_FeynRules::c_SSSS(std::vector<Single_Vertex>& vertex,int & vanz){}
void Interaction_Model_FeynRules::c_FFT(std::vector<Single_Vertex>& vertex,int& vanz)  {}
void Interaction_Model_FeynRules::c_VVT(std::vector<Single_Vertex>& vertex,int& vanz)  {}
void Interaction_Model_FeynRules::c_SST(std::vector<Single_Vertex>& vertex,int& vanz)  {}
void Interaction_Model_FeynRules::c_VVVT(std::vector<Single_Vertex>& vertex,int& vanz) {}
void Interaction_Model_FeynRules::c_FFVT(std::vector<Single_Vertex>& vertex,int& vanz) {}
void Interaction_Model_FeynRules::c_SSST(std::vector<Single_Vertex>& vertex,int& vanz) {}

Interaction_Model_FeynRules::~Interaction_Model_FeynRules()
{
  delete p_reader;
}
