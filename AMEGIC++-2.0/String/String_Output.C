#include "String_Output.H"

#ifdef Kabbala_on

#include <stdio.h>  
#include <iostream>
//#include <stdlib.h>
#include "Message.H"
#include "prof.hh"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

#define use_templates

String_Output::String_Output(const string &_path,int _maxgraph,int _maxhel) 
  : path(_path), maxgraph(_maxgraph), maxhel(_maxhel)  
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

  path=string("Process/")+path;

  msg.Debugging()<<"String_Output::String_Output("<<pathID<<"/"<<pID<<")"<<endl;
  msg.Debugging()<<"                       path= "<<path<<endl;

}

void String_Output::Output(sknot*** sk,String_Tree* stree,
			   Virtual_String_Generator* sgen,Helicity* hel)
{
  //Create names
    
  string headername = path+string("/V.H");
  if (slib.IsFile(headername)==1){
    msg.Debugging()<<"String_Output::Output   Header exists : "<<headername<<endl;
    return;
  }

  string cfilename  = path+string("/V");

  slib.InitMakefile(pathID);
  
  ofstream header;
  header.open(headername.c_str());
  Make_Header(header);

  int maxlines  = 200;
  int tolerance = 50;

  Zform(header,maxlines,tolerance,sgen,stree);

  Cform(header,maxlines,tolerance,sk,stree,cfilename,hel);

  header<<"};"<<endl;
  header<<"}"<<endl<<endl;
  header<<"#endif"<<endl;

  header.close();
  msg.Tracking()<<headername.c_str()<<" saved."<<endl;
  Add_To_Set_Values();
  /*
  AORGTOOLS::msg.Out()<<"New Libraries produced ! Please enter \"make install\" and rerun !"<<endl;
  ofstream file("profile.out");
  print_profile( file );
  abort();
  */
}

void String_Output::Cform(ofstream& header,int maxlines,int tolerance,
			  sknot*** sk,String_Tree* stree,
			  const string& cfilename,Helicity* hel)
{
  ofstream cfile;
  cfile.open((cfilename+string(".C")).c_str());

  string Makefile = string("Process/")+pathID+string("/Makefile");
  slib.AddToMakefile(Makefile,pathID,string("V"));

  int lines   = 0;
  int fnumber = 0;
  
  cfile<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
  cfile<<"using namespace AMEGIC;"<<endl;
  cfile<<"using namespace APHYTOOLS;"<<endl;
  cfile<<"using namespace AMATOOLS;"<<endl;
  cfile<<"using namespace std;"<<endl<<endl;
  cfile<<"Complex "<<pID<<"::Evaluate"<<"(int& m,int& n)"<<endl;
  cfile<<"{"<<endl;
  cfile<<"  switch (m) {"<<endl;

  lines += 11;

  for (int i=0;i<maxgraph;i++) {
    cfile<<"    case "<<i<<": return M"<<i<<"(n);"<<endl;
    lines++;
  }  
  cfile<<"  }"<<endl;
  cfile<<"  return Complex(0.,0.);"<<endl;
  cfile<<"}"<<endl<<endl;
  lines += 3;
  
  string str;
  int divnum;
  for (short int igraph=0;igraph<maxgraph;igraph++) {
    divnum = 0;
    header<<"  Complex M"<<igraph<<"(const int&);"<<endl;
    cfile<<"Complex "<<pID<<"::M"<<igraph<<"(const int& n)"<<endl;
    cfile<<"{"<<endl;
    cfile<<"  Complex M;"<<endl; 
    cfile<<"  switch (n) {"<<endl;
    lines += 4;
    for (short int ihel=0;ihel<maxhel;ihel++) {
      if (sk[igraph][ihel]!=0) 
	str = stree->Tree2String(sk[igraph][ihel],0);
      else
	str = string("");
      if (str!=string("")) {
	if (hel->On(ihel)==0) {
	  if (hel->Partner(ihel)!=-1) {
	    cfile<<"    case "<<ihel<<": M = M"<<igraph<<"("<<hel->Partner(ihel)<<");break;"<<endl;
	    lines++;
	  }
	}
	else {
	  cfile<<"    case "<<ihel<<": M = ";
	  lines += Line_Form(cfile,str);
	  cfile<<"break;"<<endl;
	  lines++;
	}
      }
      if (((maxlines+tolerance)<lines) && (ihel!=maxhel-1)) {
	lines = 0;
	divnum++;
	//close the old one
	cfile<<"    default: M = M"<<igraph<<"_"<<divnum<<"(n);"<<endl;
	cfile<<"  }"<<endl;
	cfile<<"  return M;"<<endl;
	cfile<<"}"<<endl<<endl;
	cfile.close();
	// new file
	fnumber++;
	char numb[5];
	sprintf(numb,"%i",fnumber);	
	cfile.open((cfilename+string("_")+string(numb)+string(".C")).c_str());
	slib.AddToMakefile(Makefile,pathID,string("V_")+string(numb));
	header<<"  Complex M"<<igraph<<"_"<<divnum<<"(const int&);"<<endl;
	cfile<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
	cfile<<"using namespace AMEGIC;"<<endl;  
	cfile<<"using namespace APHYTOOLS;"<<endl;
	cfile<<"using namespace AMATOOLS;"<<endl;
	cfile<<"using namespace std;"<<endl<<endl;  

	cfile<<"Complex "<<pID<<"::M"<<igraph<<"_"<<divnum<<"(const int& n)"<<endl;
	cfile<<"{"<<endl;
	cfile<<"  Complex M;"<<endl; 
	cfile<<"  switch (n) {"<<endl;
	lines += 6;
      }
    }
    cfile<<"    default: M = Complex(0.,0.);"<<endl;
    cfile<<"  }"<<endl;
    cfile<<"  return M;"<<endl;
    cfile<<"}"<<endl<<endl;
    if (((maxlines-tolerance)<lines) && igraph!=maxgraph-1) {
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
      cfile<<"using namespace APHYTOOLS;"<<endl;
      cfile<<"using namespace AMATOOLS;"<<endl;
      cfile<<"using namespace std;"<<endl<<endl;  
    }
  }

  msg.Tracking()<<cfilename<<"_"<<fnumber<<".C saved."<<endl;
  cfile.close();
}

void String_Output::Zform(ofstream& header,int &maxlines,int &tolerance,
			  Virtual_String_Generator* sgen,
			  String_Tree* stree)
{  
  string Zname = path+string("/V_Z");
  
  ofstream zf;
  zf.open((Zname+string(".C")).c_str());

  string Makefile = string("Process/")+pathID+string("/Makefile");
  slib.AddToMakefile(Makefile,pathID,string("V_Z"));

  int lines = 0;
  zf<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
  zf<<"using namespace AMEGIC;"<<endl;  
  zf<<"using namespace APHYTOOLS;"<<endl;
  zf<<"using namespace AMATOOLS;"<<endl;
  zf<<"using namespace std;"<<endl<<endl;  

  //Flavours and Couplings
  header<<"  void SetCouplFlav();"<<endl;
  zf<<"void "<<pID<<"::SetCouplFlav()"<<endl;
  zf<<"{"<<endl;
  zf<<"  if (f==NULL) f = new int["<<sgen->GetFlavours()->size()<<"];"<<endl<<endl;
  lines += 3;
  for (short int i=0;i<sgen->GetFlavours()->size();i++) {
    zf<<"  f["<<i<<"] = "<<(*sgen->GetFlavours())[i]<<";"<<endl;
    lines++;
  }
  zf<<"}"<<endl<<endl;
  lines +=2 ;

  zf<<"void "<<pID<<"::Calculate(vector<Complex>& c)"<<endl;
  zf<<"{"<<endl;
  int hit;
  Complex norm;
  zf<<"  if (Z==NULL) Z = new Complex["<<sgen->ZX_Max_Number()<<"];"<<endl<<endl;
  zf<<"  Z[0] = Complex(0.,0.);"<<endl;
  lines += 7;

  int divnum = 0;
  ZXlist* zx;
  for (int i=1;i<sgen->ZX_Max_Number();i++) {
    zx = sgen->Get_ZXl(i);
    if (zx->on) {
      lines++;
      zf<<"  Z["<<i<<"] = ";
      int* arg = zx->arg;
      switch (zx->zlist) {
      case 0: 
#ifdef use_templates
	zf<<"XT<"<<arg[1]<<","<<arg[4]<<">";
	zf<<"("<<arg[0]<<","<<arg[2]<<","<<arg[3];
	zf<<",c["<<arg[5]<<"],c["<<arg[6]<<"]);"<<endl;
#else
	zf<<"Xcalc("<<arg[0]<<","<<arg[1]<<","<<arg[2]<<","<<arg[3]<<","<<arg[4];
	zf<<",c["<<arg[5]<<"],c["<<arg[6]<<"]);"<<endl;
#endif
	break;
      case 1:
#ifdef use_templates
	zf<<"ZT";
	if (!sgen->Massless(i)) zf<<"M";
	zf<<"<"<<arg[1]<<","<<arg[3]<<","<<arg[5]<<","<<arg[7]<<">";
	zf<<"("<<arg[0]<<","<<arg[2]<<","<<arg[4]<<","<<arg[6];
	zf<<",c["<<arg[8]<<"],c["<<arg[9]<<"],c["<<arg[10]<<"],c["<<arg[11]<<"]);"<<endl;
#else
	zf<<"Zcalc("<<arg[0]<<","<<arg[1]<<","<<arg[2]<<","<<arg[3];
	zf<<","<<arg[4]<<","<<arg[5]<<","<<arg[6]<<","<<arg[7];
	zf<<",c["<<arg[8]<<"],c["<<arg[9]<<"],c["<<arg[10]<<"],c["<<arg[11]<<"]);"<<endl;
#endif
	break;
      case 2: 
	// E = 1/M(Z,W,H)^(2,4); E = coupl
	norm = zx->value.Value();
	hit = 0;
	//couplings
	for (short int j=0;j<sgen->Coupl_Max_Number();j++) {
	  if ( AMATOOLS::IsEqual(norm,sgen->Get_Coupl(j)) ||
	       AMATOOLS::IsEqual(norm,-sgen->Get_Coupl(j)) ) {
	    hit = 1;
	    if (AMATOOLS::IsEqual(norm,-sgen->Get_Coupl(j))) zf<<"-";
	    zf<<"c["<<j<<"];"<<endl;
	    break;
	  }
	}
	if (hit) break;
	//masses
	if (real(norm)<0) {
	  zf<<"-";
	  norm = -norm;
	}
	if (AMATOOLS::IsEqual(norm,1./sqr(Flavour(kf::Z).mass()))) { 
	  hit = 1;
	  zf<<"Complex(1./sqr(Flavour(kf::Z).mass()),0.);"<<endl;
	  break;
	}
	//new
	if (AMATOOLS::IsEqual(norm,1./(Complex(sqr(Flavour(kf::Z).mass()),
			      -Flavour(kf::Z).mass()*Flavour(kf::Z).width())))) { 
	    hit = 1;
	    zf<<"(1./Complex(sqr(Flavour(kf::Z).mass()),"
	      <<"-Flavour(kf::Z).mass()*Flavour(kf::Z).width()));"<<endl;
	    break;
	}
	if (AMATOOLS::IsEqual(norm,1./sqr(Flavour(kf::W).mass()))) { 
	  hit = 1;
	  zf<<"Complex(1./sqr(Flavour(kf::W).mass()),0.);"<<endl;
	  break;
	}
	//new
	if (AMATOOLS::IsEqual(norm,1./(Complex(sqr(Flavour(kf::W).mass()),
			      -Flavour(kf::W).mass()*Flavour(kf::W).width())))) { 
	    hit = 1;
	    zf<<"(1./Complex(sqr(Flavour(kf::W).mass()),"
	      <<"-Flavour(kf::W).mass()*Flavour(kf::W).width()));"<<endl;
	    break;
	}
	if (AMATOOLS::IsEqual(norm,1./sqr(Flavour(kf::h).mass()))) { 
	  hit = 1;
	  zf<<"Complex(1./sqr(Flavour(kf::h).mass()),0.);"<<endl;
	  break;
	}
	if (AMATOOLS::IsEqual(norm,1./sqr(sqr(Flavour(kf::Z).mass())))) { 
	  hit = 1;
	  zf<<"Complex(1./sqr(sqr(Flavour(kf::Z).mass())),0.);"<<endl;
	  break;
	}	  
	// double masses
	if (AMATOOLS::IsEqual(norm,1./sqr(Flavour(kf::Z).mass()*Flavour(kf::W).mass()))) { 
	  hit = 1;
	  zf<<"Complex(1./sqr(Flavour(kf::Z).mass()*Flavour(kf::W).mass()),0.);"<<endl;
	  break;	
	}
	if (AMATOOLS::IsEqual(norm,1./sqr(Flavour(kf::W).mass()*Flavour(kf::W).mass()))) { 
	  hit = 1;
	  zf<<"Complex(1./sqr(Flavour(kf::W).mass()*Flavour(kf::W).mass()),0.);"<<endl;
	  break;
	}
	if (AMATOOLS::IsEqual(norm,0.5)) {
          hit = 1;
          zf<<"Complex(0.5,0.);"<<endl;
          break;
	}
	if (AMATOOLS::IsEqual(norm,-0.5)) {
          hit = 1;
          zf<<"Complex(-0.5,0.);"<<endl;
	  break;        
	}
	if (AMATOOLS::IsEqual(norm,1.)) { 
	  hit = 1;
	  zf<<"Complex(1.,0.);"<<endl;
	    break;
	}
	if (AMATOOLS::IsEqual(norm,2.)) { 
	  hit = 1;
	  zf<<"Complex(2.,0.);"<<endl;
	  break;
	}
	if (AMATOOLS::IsEqual(norm,Complex(0.,1.))) { 
	  hit = 1;
	  zf<<"Complex(0.,1.);"<<endl;
	  break;
	}
	if (AMATOOLS::IsEqual(norm,Complex(0.,-1.))) { 
	  hit = 1;
	  zf<<"Complex(0.,-1.);"<<endl;
	  break;
	}
	if (AMATOOLS::IsEqual(norm,Complex(0.,-1./4.))) { 
	  hit = 1;
	  zf<<"Complex(0.,-1./4.);"<<endl;
	  break;
	}
	if (hit==0) msg.Error()<<"No match for E-function:"<<zx->value.Value()<<endl;
	break;
      case 3: 
	zf<<"Vcalc("<<arg[0]<<","<<arg[1]<<");"<<endl;
	break;
      case 4: 
	zf<<"Ycalc("<<arg[0]<<","<<arg[1]<<","<<arg[2]<<","<<arg[3];
	zf<<",c["<<arg[4]<<"],c["<<arg[5]<<"]);"<<endl;
	break;
      case 5: 
	zf<<"Pcalc(f["<<arg[0]<<"],"<<arg[1]<<");"<<endl;
	break;				       
      case 6:
	if (zx->sk!=0) {
	  lines += Line_Form(zf,stree->Tree2String(zx->sk,0));
	  zf<<endl;
	}
	else 
	  zf<<"Complex(0.,0.);"<<endl;
	break;
      case 7: 
	zf<<"MassTermCalc("<<arg[1]<<",f["<<arg[0]<<"]);"<<endl;
	break;
	
      case 8: 
	zf<<"VNCcalc("<<arg[0]<<","<<arg[1]<<");"<<endl;
 	break;
      case 9: 
 	zf<<"Vcplxcalc("<<arg[0]<<","<<arg[1]<<");"<<endl;
 	break;
      }
    }
    if (((maxlines+tolerance)<lines) && (i!=sgen->ZX_Max_Number()-1)) {
      lines = 0;
      divnum++;
      //close the old one
      zf<<"  Calculate_"<<divnum<<"(c);"<<endl;
      zf<<"}"<<endl;
      zf.close();
      // new file
      char numb[5];
      sprintf(numb,"%i",divnum);
      zf.open((Zname+string("_")+string(numb)+string(".C")).c_str());
      slib.AddToMakefile(Makefile,pathID,string("V_Z_")+string(numb));
      zf<<"#include "<<'"'<<"V.H"<<'"'<<endl<<endl;  
      zf<<"using namespace AMEGIC;"<<endl;  
      zf<<"using namespace APHYTOOLS;"<<endl;
      zf<<"using namespace AMATOOLS;"<<endl;
      zf<<"using namespace std;"<<endl<<endl;  
      zf<<"void "<<pID.c_str()<<"::Calculate_"<<divnum
	<<"(vector<Complex>& c)"<<endl;
      zf<<"{"<<endl;
      header<<"  void Calculate_"<<divnum<<"(std::vector<Complex>&);"<<endl;
      lines += 4;
    }
  }
  zf<<"}"<<endl;
  zf.close();
  msg.Tracking()<<Zname<<"_"<<divnum<<".C saved."<<endl;
}

void String_Output::Make_Header(ofstream &header)
{
  header<<"//Header file for process "<<pID.c_str()<<endl<<endl;
  header<<"#ifndef "<<pID.c_str()<<"_on"<<endl;
  header<<"#define "<<pID.c_str()<<"_on"<<endl;  
  header<<"#include "<<'"'<<"Values.H"<<'"'<<endl<<endl;
  header<<"namespace AMEGIC {"<<endl<<endl;  
  header<<"class "<<pID.c_str()<<" : public Values,"<<endl;
  header<<"  public Basic_Yfunc,"<<endl; 
  header<<"  public Basic_Zfunc,"<<endl; 
  header<<"  public Basic_Xfunc,"<<endl; 
  header<<"  public Basic_Mfunc,"<<endl; 
  header<<"  public Basic_Vfunc,"<<endl; 
  header<<"  public Basic_Pfunc,"<<endl; 
  header<<"  public Basic_MassTermfunc {"<<endl; 

  header<<"  Complex* Z;"<<endl; 
  header<<"  int*     f;"<<endl; 
  header<<"  Complex* c;"<<endl; 
  header<<"public:"<<endl;
  header<<"  "<<pID<<"(Basic_Sfuncs* _BS) :"<<endl; 
  header<<"     Basic_Func(0,_BS),"<<endl;  
  header<<"     Basic_Yfunc(0,_BS),"<<endl;  
  header<<"     Basic_Zfunc(0,_BS),"<<endl;  
  header<<"     Basic_Xfunc(0,_BS),"<<endl;  
  header<<"     Basic_Mfunc(0,_BS),"<<endl;  
  header<<"     Basic_Vfunc(0,_BS),"<<endl;  
  header<<"     Basic_Pfunc(0,_BS),"<<endl; 
  header<<"     Basic_MassTermfunc(0,_BS) {Z=NULL;f=NULL;c=NULL;}"<<endl;
  header<<"  ~"<<pID<<"() {"<<endl;
  header<<"    if (Z) delete[] Z;"<<endl;
  header<<"    if (f) delete[] f;"<<endl;
  header<<"    if (c) delete[] c;"<<endl;
  header<<"  }"<<endl;
  header<<"  Complex Evaluate(int&,int&);"<<endl;
  header<<"  void    Calculate(std::vector<Complex>&);"<<endl;
}

int String_Output::Line_Form(ofstream& file,const string &str)
{
  int counter = 0;
  int lines   = 0;
  for (short int j=0;j<str.length();j++) {
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

  ifstream from;
  ofstream to;

  from.open("Process/Set_Values.C");
  to.open("Process/Set_Values.C.tmp");

  int hit = 0;

  char buffer[512];

  //include into first line
  to<<"#include "<<'"'<<(pathID+string("/V.H")).c_str()<<'"'<<endl;  

  for(;from;) {
    from.getline(buffer,512);
    
    if (string(buffer).find(pID)!=-1) break;

    if (string(buffer).find(string("return 0"))!=-1) {
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
    slib.Copy(string("Process/Set_Values.C.tmp"),string("Process/Set_Values.C"));
  else 
    remove("Process/Set_Values.C.tmp");
}

#endif
