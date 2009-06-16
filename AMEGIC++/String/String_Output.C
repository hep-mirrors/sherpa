#include "AMEGIC++/String/String_Output.H"

#include <stdio.h>  
#include <iostream>
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

#define use_templates

String_Output::String_Output(const string &_path,int _maxgraph,int _maxhel, int mode) 
  : slib(mode), path(_path), maxgraph(_maxgraph), maxhel(_maxhel), m_mode(mode)
{
  pathID = path;
  pID    = path;

  for (short int i=path.length()-1;i>=0;i--) {
    if (path[i]=='/') {
      pID    = string("V")+path.substr(i+1);
      //      pathID = path.substr(i+1);
      break;
    }
  }
  //kill +- in ID
  string help;
  short int i;
  for (;;) {
    i = pID.find("+");
    if (i==-1) i = pID.find("-");
    if (i==-1) break;
    help = pID.substr(0,i) + pID.substr(i+1);
    pID = help;
  }

  path=rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+path;
}

void String_Output::Output(sknot*** sk,String_Tree* stree,
			   Virtual_String_Generator* sgen,Helicity* hel)
{
  //Create names
    
  string headername = path+string("/V.H");
  if (slib.IsFile(headername)==1){
    return;
  }

  string cfilename  = path+string("/V");

  slib.InitMakefile(pathID);
  
  ofstream header;
  header.open(headername.c_str());
  Make_Header(header,sgen);

  int maxlines  = 200;
  int tolerance = 50;

  Zform(header,maxlines,tolerance,sgen,stree);

  Cform(header,maxlines,tolerance,sgen,sk,stree,cfilename,hel);

  header<<"};"<<endl;
  header<<"}"<<endl<<endl;
  header<<"#endif"<<endl;

  header.close();
  Add_To_Set_Values();
}

void String_Output::Cform(ofstream& header,int maxlines,int tolerance,
			  Virtual_String_Generator* sgen,
			  sknot*** sk,String_Tree* stree,
			  const string& cfilename,Helicity* hel)
{
  ofstream cfile;
  cfile.open((cfilename+string(".C")).c_str());

  string Makefile = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+pathID+string("/Makefile");
  slib.AddToMakefile(Makefile,pathID,string("V"));

  int lines   = 0;
  int fnumber = 0;
  
  cfile<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
  cfile<<"using namespace AMEGIC;"<<endl;
  cfile<<"using namespace ATOOLS;"<<endl;
  cfile<<"using namespace std;"<<endl<<endl;

  cfile<<"extern "<<'"'<<"C"<<'"'<<" Values* Getter_"<<pID<<"(Basic_Sfuncs* bs) {"<<endl;
  cfile<<"  return new "<<pID<<"(bs);"<<endl;
  cfile<<"}"<<endl<<endl;

  cfile<<pID<<"::"<<pID<<"(Basic_Sfuncs* _BS) :"<<endl;
  cfile<<"     Basic_Func(0,_BS)";  
  if (sgen->UsesFunction(4)) cfile<<","<<endl<<"     Basic_Yfunc(0,_BS)";  
  if (sgen->UsesFunction(1)) cfile<<","<<endl<<"     Basic_Zfunc(0,_BS)";  
  if (sgen->UsesFunction(0)) cfile<<","<<endl<<"     Basic_Xfunc(0,_BS)";  
  if (sgen->UsesFunction(3)||sgen->UsesFunction(9)) cfile<<","<<endl<<"     Basic_Vfunc(0,_BS)";  
  if (sgen->UsesFunction(5)) cfile<<","<<endl<<"     Basic_Pfunc(0,_BS)"; 
  if (sgen->UsesFunction(10)) cfile<<","<<endl<<"     Basic_Epsilonfunc(0,_BS)";
  if (sgen->UsesFunction(7)) cfile<<","<<endl<<"     Basic_MassTermfunc(0,_BS)";
  if (sgen->UsesFunction(11)) cfile<<","<<endl<<"     Unitarityfunc(0,_BS)";
  cfile<<endl<<"{"<<endl;
  cfile<<"  f = new int["<<sgen->GetFlavours()->size()<<"];"<<endl;
  cfile<<"  c = new Complex["<<sgen->NumberOfCouplings()<<"];"<<endl;
  cfile<<"  Z = new Complex["<<sgen->ZXMaxNumber()<<"];"<<endl;
  cfile<<"  M = new Complex*["<<maxhel<<"];"<<endl;
  cfile<<"  for(int i=0;i<"<<maxhel<<";i++) M[i] = new Complex["<<maxgraph<<"];"<<endl;
  cfile<<"  cl = new int["<<maxhel<<"];"<<endl;
  cfile<<"}"<<endl<<endl;
  cfile<<pID<<"::"<<"~"<<pID<<"()"<<endl; 
  cfile<<"{"<<endl;
  cfile<<"  if (Z)  delete[] Z;"<<endl;
  cfile<<"  if (f)  delete[] f;"<<endl;
  cfile<<"  if (c)  delete[] c;"<<endl;
  cfile<<"  if (cl) delete[] cl;"<<endl;
  cfile<<"  if (M) {"<<endl;
  cfile<<"    for(int i=0;i<"<<maxhel<<";i++) delete[] M[i];"<<endl;
  cfile<<"    delete[] M;"<<endl; 
  cfile<<"  }"<<endl;
  cfile<<"}"<<endl<<endl;
  lines+=25;

  cfile<<"Complex "<<pID<<"::Evaluate"<<"(int m,int n)"<<endl;
  cfile<<"{"<<endl;
  cfile<<"  if (cl[n]) return M[n][m];"<<endl;
  cfile<<"  switch (n) {"<<endl;

  lines += 4;

  for (int ihel=0;ihel<maxhel;ihel++) {
    if (hel->On(ihel)!=0) {
      header<<"  void Calculate_M"<<ihel<<"();"<<endl;
      cfile<<"    case "<<ihel<<": Calculate_M"<<ihel<<"(); break;"<<endl;
      lines++;
    }
    if (hel->On(ihel)==0) 
      if (hel->Partner(ihel)!=-1) {
	cfile<<"    case "<<ihel<<": return Evaluate(m,"<<hel->Partner(ihel)<<");"<<endl;
	lines++;
      }
  }  
  cfile<<"  }"<<endl;
  cfile<<"  cl[n]=1;"<<endl;
  cfile<<"  return M[n][m];"<<endl;
  cfile<<"}"<<endl<<endl;
  lines += 4;
  
  string str;
  int divnum;
  for (int ihel=0;ihel<maxhel;ihel++) {
    if (hel->On(ihel)!=0) {
      if ((maxlines-tolerance)<lines) {
	//cut
	lines = 0;
	cfile.close();
	// new file
	fnumber++;
	char numb[5];
	sprintf(numb,"%i",fnumber);	
	cfile.open((cfilename+string("_")+string(numb)+string(".C")).c_str());
	slib.AddToMakefile(Makefile,pathID,string("V_")+string(numb));
	cfile<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
	cfile<<"using namespace AMEGIC;"<<endl;  
	cfile<<"using namespace ATOOLS;"<<endl<<endl;
      }

      divnum = 0;
      cfile<<"void "<<pID<<"::Calculate_M"<<ihel<<"()"<<endl;
      cfile<<"{"<<endl;
      lines += 2;

      for (int igraph=0;igraph<maxgraph;igraph++) {
	if (sk[igraph][ihel]!=0) 
	  str = stree->Tree2String(sk[igraph][ihel],0);
	else
	  str = string("");
	if (str!=string("")) {
	  cfile<<"  M["<<ihel<<"]["<<igraph<<"] = ";
	  lines += Line_Form(cfile,str)+1;
	}
      
	if (((maxlines+tolerance)<lines) && (ihel!=maxhel-1)) {
	  lines = 0;
	  divnum++;
	  //close the old one
	  cfile<<"  Calculate_M"<<ihel<<"_"<<divnum<<"();"<<endl;
	  cfile<<"}"<<endl<<endl;
	  cfile.close();
	  // new file
	  fnumber++;
	  char numb[5];
	  sprintf(numb,"%i",fnumber);	
	  cfile.open((cfilename+string("_")+string(numb)+string(".C")).c_str());
	  slib.AddToMakefile(Makefile,pathID,string("V_")+string(numb));
	  header<<"  void Calculate_M"<<ihel<<"_"<<divnum<<"();"<<endl;
	  cfile<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
	  cfile<<"using namespace AMEGIC;"<<endl;  
	  cfile<<"using namespace ATOOLS;"<<endl<<endl;
	  cfile<<"void "<<pID<<"::Calculate_M"<<ihel<<"_"<<divnum<<"()"<<endl;
	  cfile<<"{"<<endl;
	  lines += 2;
	}
      }
      cfile<<"}"<<endl<<endl;
    }
  }
  cfile.close();
}

void String_Output::Zform(ofstream& header,int maxlines,int tolerance,
			  Virtual_String_Generator* sgen,
			  String_Tree* stree)
{  
  string Zname = path+string("/V_Z");
  
  ofstream zf,szf;
  zf.open((Zname+string(".C")).c_str());

  string Makefile = rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/")+pathID+string("/Makefile");
  slib.AddToMakefile(Makefile,pathID,string("V_Z"));

  int lines = 0;
  zf<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
  zf<<"using namespace AMEGIC;"<<endl;  
  zf<<"using namespace ATOOLS;"<<endl;
  zf<<"using namespace std;"<<endl<<endl;  

  //Flavours and Couplings
  zf<<"void "<<pID<<"::SetCouplFlav(vector<Complex>& coupl)"<<endl;
  zf<<"{"<<endl;
  for (size_t i=0;i<sgen->GetFlavours()->size();i++) {
    zf<<"  f["<<i<<"] = "<<(*sgen->GetFlavours())[i]<<";"<<endl;
  }
  zf<<endl;
  zf<<"  for (int i=0;i<"<<sgen->NumberOfCouplings()<<";i++) "<<"c[i] = coupl[i];"<<endl;
  zf<<"  for (int i=0;i<"<<maxhel<<";i++)"<<endl;
  zf<<"    for (int j=0;j<"<<maxgraph<<";j++) M[i][j] = Complex(0.,0.);"<<endl<<endl;
  zf<<"  Z[0] = Complex(0.,0.);"<<endl;
  zf<<"}"<<endl<<endl;

  zf<<"void "<<pID<<"::Calculate()"<<endl;
  zf<<"{"<<endl;
  zf<<"  for(int i=0;i<"<<maxhel<<";i++) cl[i] = 0;"<<endl<<endl;

  bool mvz = (sgen->ZXMaxNumber()>maxlines+tolerance);

  ofstream* pz=&zf;
  int divnum = 0;
  if (mvz) { 
    divnum = 1;
    szf.open((Zname+string("_")+ToString(divnum)+string(".C")).c_str());
    slib.AddToMakefile(Makefile,pathID,string("V_Z_")+ToString(divnum));
    szf<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
    szf<<"using namespace AMEGIC;"<<endl;  
    szf<<"using namespace ATOOLS;"<<endl<<endl;
    szf<<"void "<<pID.c_str()<<"::Calculate_"<<divnum
       <<"()"<<endl;
    szf<<"{"<<endl;

    header<<"  void Calculate_"<<divnum<<"();"<<endl;
    zf<<"  Calculate_"<<divnum<<"();"<<endl;

    pz=&szf;
  }

  ZXlist* zx;
  int hit;
  Complex norm;
  for (int i=1;i<sgen->ZXMaxNumber();i++) {
    zx = sgen->GetZXl(i);
    if (zx->on) {
      lines++;
      (*pz)<<"  Z["<<i<<"] = ";
      int* arg = zx->arg;
      switch (zx->zlist) {
      case 0: 
#ifdef use_templates
	(*pz)<<"XT<"<<arg[1]<<","<<arg[4]<<">";
	(*pz)<<"("<<arg[0]<<","<<arg[2]<<","<<arg[3];
	(*pz)<<",c["<<arg[5]<<"],c["<<arg[6]<<"]);"<<endl;
#else
	(*pz)<<"Xcalc("<<arg[0]<<","<<arg[1]<<","<<arg[2]<<","<<arg[3]<<","<<arg[4];
	(*pz)<<",c["<<arg[5]<<"],c["<<arg[6]<<"]);"<<endl;
#endif
	break;
      case 1:
#ifdef use_templates
	(*pz)<<"ZT";
	if (!sgen->Massless(i)) (*pz)<<"M";
	(*pz)<<"<"<<arg[1]<<","<<arg[3]<<","<<arg[5]<<","<<arg[7]<<">";
	(*pz)<<"("<<arg[0]<<","<<arg[2]<<","<<arg[4]<<","<<arg[6];
	(*pz)<<",c["<<arg[8]<<"],c["<<arg[9]<<"],c["<<arg[10]<<"],c["<<arg[11]<<"]);"<<endl;
#else
	(*pz)<<"Zcalc("<<arg[0]<<","<<arg[1]<<","<<arg[2]<<","<<arg[3];
	(*pz)<<","<<arg[4]<<","<<arg[5]<<","<<arg[6]<<","<<arg[7];
	(*pz)<<",c["<<arg[8]<<"],c["<<arg[9]<<"],c["<<arg[10]<<"],c["<<arg[11]<<"]);"<<endl;
#endif
	break;
      case 2: 
	norm = zx->value.Value();
	hit = 0;
	//couplings
	for (short int j=0;j<sgen->NumberOfCouplings();j++) {
	  if ( ATOOLS::IsEqual(norm,sgen->GetCoupling(j)) ||
	       ATOOLS::IsEqual(norm,-sgen->GetCoupling(j)) ) {
	    hit = 1;
	    if (ATOOLS::IsEqual(norm,-sgen->GetCoupling(j))) (*pz)<<"-";
	    (*pz)<<"c["<<j<<"];"<<endl;
	    break;
	  }
	}
	if (hit) break;
	//masses
	if (real(norm)<0) {
	  (*pz)<<"-";
	  norm = -norm;
	}
	if (ATOOLS::IsEqual(norm,1./sqr(Flavour(kf_Z).Mass()))) { 
	  hit = 1;
	  (*pz)<<"Complex(1./sqr(Flavour(kf_Z).Mass()),0.);"<<endl;
	  break;
	}
	//new
	if (ATOOLS::IsEqual(norm,sqr(Flavour(kf_Z).Mass()))) { 
	  hit = 1;
	  (*pz)<<"Complex(sqr(Flavour(kf_Z).Mass()),0.);"<<endl;
	  break;
	}
	//new
	if (ATOOLS::IsEqual(norm,1./sqr(Flavour(kf_Wplus).Mass()))) { 
	  hit = 1;
	  (*pz)<<"Complex(1./sqr(Flavour(kf_Wplus).Mass()),0.);"<<endl;
	  break;
	}
	//new
	if (ATOOLS::IsEqual(norm,1./sqr(Flavour(kf_h0).Mass()))) { 
	  hit = 1;
	  (*pz)<<"Complex(1./sqr(Flavour(kf_h0).Mass()),0.);"<<endl;
	  break;
	}
	if (ATOOLS::IsEqual(norm,1./sqr(sqr(Flavour(kf_Z).Mass())))) { 
	  hit = 1;
	  (*pz)<<"Complex(1./sqr(sqr(Flavour(kf_Z).Mass())),0.);"<<endl;
	  break;
	}	  
	// double masses
	if (ATOOLS::IsEqual(norm,1./sqr(Flavour(kf_Z).Mass()*Flavour(kf_Wplus).Mass()))) { 
	  hit = 1;
	  (*pz)<<"Complex(1./sqr(Flavour(kf_Z).Mass()*Flavour(kf_Wplus).Mass()),0.);"<<endl;
	  break;	
	}
	if (ATOOLS::IsEqual(norm,1./sqr(Flavour(kf_Wplus).Mass()*Flavour(kf_Wplus).Mass()))) { 
	  hit = 1;
	  (*pz)<<"Complex(1./sqr(Flavour(kf_Wplus).Mass()*Flavour(kf_Wplus).Mass()),0.);"<<endl;
	  break;
	}
	if (ATOOLS::IsEqual(norm,0.5)) {
          hit = 1;
          (*pz)<<"Complex(0.5,0.);"<<endl;
          break;
	}
	if (ATOOLS::IsEqual(norm,-0.5)) {
          hit = 1;
          (*pz)<<"Complex(-0.5,0.);"<<endl;
	  break;        
	}
	if (ATOOLS::IsEqual(norm,1./3.)) {
          hit = 1;
          (*pz)<<"Complex(1./3.,0.);"<<endl;
          break;
	}
	if (ATOOLS::IsEqual(norm,1.)) { 
	  hit = 1;
	  (*pz)<<"Complex(1.,0.);"<<endl;
	    break;
	}
	if (ATOOLS::IsEqual(norm,2.)) { 
	  hit = 1;
	  (*pz)<<"Complex(2.,0.);"<<endl;
	  break;
	}
	if (ATOOLS::IsEqual(norm,Complex(0.,1.))) { 
	  hit = 1;
	  (*pz)<<"Complex(0.,1.);"<<endl;
	  break;
	}
	if (ATOOLS::IsEqual(norm,Complex(0.,-1.))) { 
	  hit = 1;
	  (*pz)<<"Complex(0.,-1.);"<<endl;
	  break;
	}
	if (ATOOLS::IsEqual(norm,Complex(0.,-1./4.))) { 
	  hit = 1;
	  (*pz)<<"Complex(0.,-1./4.);"<<endl;
	  break;
	}
	if (hit==0) {
	  msg_Error()<<"No match for E-function:"<<zx->value.Value()<<endl;
	  abort();
	}
	break;
      case 3: 
	(*pz)<<"Vcalc("<<arg[0]<<","<<arg[1]<<");"<<endl;
	break;
      case 4: 
#ifdef use_templates
	(*pz)<<"YT<"<<arg[1]<<","<<arg[3]<<">";
	(*pz)<<"("<<arg[0]<<","<<arg[2];
	(*pz)<<",c["<<arg[4]<<"],c["<<arg[5]<<"]);"<<endl;
#else
	(*pz)<<"Ycalc("<<arg[0]<<","<<arg[1]<<","<<arg[2]<<","<<arg[3];
	(*pz)<<",c["<<arg[4]<<"],c["<<arg[5]<<"]);"<<endl;
#endif
	break;
      case 5: 
	(*pz)<<"Pcalc(f["<<arg[0]<<"],"<<arg[1]<<");"<<endl;
	break;				       
      case 6:
	if (zx->sk!=0) {
	  lines += Line_Form((*pz),stree->Tree2String(zx->sk,0));
	  (*pz)<<endl;
	}
	else 
	  (*pz)<<"Complex(0.,0.);"<<endl;
	break;
      case 7: 
	(*pz)<<"MassTermCalc("<<arg[1]<<",f["<<arg[0]<<"]);"<<endl;
	break;
      case 8: 
	(*pz)<<"(1./Complex(sqr(Flavour("<<arg[0]<<").Mass()),"
	     <<"-Flavour("<<arg[0]<<").Mass()*Flavour("<<arg[0]<<").Width()));"<<endl;
	break;
      case 9: 
 	(*pz)<<"Vcplxcalc("<<arg[0]<<","<<arg[1]<<");"<<endl;
 	break;
      case 10: 
 	(*pz)<<"EpsCalc<"<<arg[4]<<">("<<arg[0]<<","<<arg[1]<<","<<arg[2]<<","<<arg[3]<<");"<<endl;
 	break;
      case 11: 
 	(*pz)<<"Ucalc();"<<endl;
 	break;
      }
    }
    if (mvz && ((maxlines+tolerance)<lines) && (i!=sgen->ZXMaxNumber()-1)) {
      lines = 0;
      divnum++;
      //close the old one
      (*pz)<<"}"<<endl;
      (*pz).close();
      // new file
      (*pz).open((Zname+string("_")+ToString(divnum)+string(".C")).c_str());
      slib.AddToMakefile(Makefile,pathID,string("V_Z_")+ToString(divnum));
      (*pz)<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
      (*pz)<<"using namespace AMEGIC;"<<endl;  
      (*pz)<<"using namespace ATOOLS;"<<endl<<endl;
      (*pz)<<"void "<<pID.c_str()<<"::Calculate_"<<divnum
	<<"()"<<endl;
      (*pz)<<"{"<<endl;
      header<<"  void Calculate_"<<divnum<<"();"<<endl;
      zf<<"  Calculate_"<<divnum<<"();"<<endl;
      lines += 4;
    }
  }
  zf<<"}"<<endl;
  zf.close();
  if (mvz) {
    (*pz)<<"}"<<endl;
    (*pz).close();
  }
}

void String_Output::Make_Header(ofstream &header,Virtual_String_Generator* sgen)
{
  header<<"//Header file for process "<<pID.c_str()<<endl<<endl;
  header<<"#ifndef "<<pID<<"_on"<<endl;
  header<<"#define "<<pID<<"_on"<<endl;  
  header<<"#include "<<'"'<<"AMEGIC++/String/Values.H"<<'"'<<endl<<endl;

  header<<"extern "<<'"'<<"C"<<'"'<<" AMEGIC::Values* Getter_"<<pID<<"(AMEGIC::Basic_Sfuncs* bs);"<<endl<<endl;

  header<<"namespace AMEGIC {"<<endl<<endl;  
  header<<"class "<<pID.c_str()<<" : public Values";
  if (sgen->UsesFunction(4)) header<<","<<endl<<"  public Basic_Yfunc"; 
  if (sgen->UsesFunction(1)) header<<","<<endl<<"  public Basic_Zfunc"; 
  if (sgen->UsesFunction(0)) header<<","<<endl<<"  public Basic_Xfunc"; 
  if (sgen->UsesFunction(3)||sgen->UsesFunction(9)) header<<","<<endl<<"  public Basic_Vfunc"; 
  if (sgen->UsesFunction(5)) header<<","<<endl<<"  public Basic_Pfunc"; 
  if (sgen->UsesFunction(10)) header<<","<<endl<<"  public Basic_Epsilonfunc";
  if (sgen->UsesFunction(7)) header<<","<<endl<<"  public Basic_MassTermfunc";
  if (sgen->UsesFunction(11)) header<<","<<endl<<"  public Unitarityfunc";
  header<<" {"<<endl; 

  header<<"  Complex*  Z;"<<endl; 
  header<<"  int*      f;"<<endl; 
  header<<"  Complex*  c;"<<endl; 
  header<<"  Complex** M;"<<endl; 
  header<<"  int*      cl;"<<endl; 
  header<<"public:"<<endl;
  header<<"  "<<pID<<"(Basic_Sfuncs* _BS);"<<endl; 
  header<<"  ~"<<pID<<"();"<<endl;
  header<<"  void SetCouplFlav(std::vector<Complex>&);"<<endl;
  header<<"  int NumberOfCouplings() { return "<<sgen->NumberOfCouplings()<<"; }"<<endl;
  header<<"  Complex Evaluate(int,int);"<<endl;
  header<<"  void    Calculate();"<<endl;
}

int String_Output::Line_Form(ofstream& file,const string &str)
{
  int counter = 0;
  int lines   = 0;
  for (size_t j=0;j<str.length();j++) {
    if (counter>70) {
      int hit = 0;
      switch (str[j]) {
      case '+':hit = 1;break;
      case '*':hit = 1;break;
      case '-':hit = 1;break;
      }
      if (hit) {
	file<<endl<<"           ";
	lines++;
	counter = 0;
      }
      file<<str[j];
    }
    else {
      counter++;
      file<<str[j];
    }
  }
  file<<";"<<endl;

  return lines;
}

void String_Output::Add_To_Set_Values()
{
  //manipulate Set_Values

  // only include in Set_Values.C if neccessary
  if (m_mode==1) return;


  ifstream from;
  ofstream to;

  from.open((rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/Set_Values.C")).c_str());
  to.open((rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/Set_Values.C.tmp")).c_str());

  int hit = 0;

  string buffer;

  //include into first line
  to<<"#include "<<'"'<<(pathID+string("/V.H")).c_str()<<'"'<<endl;  

  for(;from;) {
    getline(from,buffer);
    
    if (buffer.find(pID)!=string::npos) break;

    if (buffer.find(string("return 0"))!=string::npos || 
	buffer.find(string("libname"))!=string::npos) {
      hit = 1;
      to<<"#ifdef "<<pID<<"_on"<<endl;
      to<<"  if (pID==string("<<'"'<<pID<<'"'<<")) return (new "<<pID<<"(BS));"<<endl;
      to<<"#endif"<<endl;
    }
    to<<buffer<<endl;
  }
  from.close();
  to.close();

  if (hit) 
    slib.Copy(rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/Set_Values.C.tmp"),rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/Set_Values.C"));
  else 
    remove((rpa.gen.Variable("SHERPA_CPP_PATH")+string("/Process/Set_Values.C.tmp")).c_str());
}

