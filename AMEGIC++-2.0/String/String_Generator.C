#include <iostream>
#include "String_Generator.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "MyStrStream.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace std;
#ifdef __GNUC__
#include <stdio.h>
#endif

String_Generator::~String_Generator() {}

void String_Generator::Reset()
{
  couplings.clear();
  flavours.clear();
  zxl.clear();

  ZXlist zero;
  // Z0 == 0
#ifdef Kabbala_on
  zero.value = Kabbala(string("Z[0]"),Complex(0.,0.));
#else
  zero.value = Complex(0.,0.);
#endif    
  zero.zlist = -1;

  zxl.push_back(zero);
}

void String_Generator::Print()
{
  if (!(rpa.gen.Debugging())) return;
  for (long int i=0;i<zxl.size();i++) {
    msg.Out()<<i<<". Zfunction: Type="<<zxl[i].zlist<<";On="<<zxl[i].on<<";Value="<<zxl[i].value.String(); 
    if (zxl[i].narg>0) msg.Out()<<";Arg[0] = "<<zxl[i].arg[0];
    msg.Out()<<endl;
  }
}

int String_Generator::Z_Count() 
{
  int count = 0;
  for (long int i=1;i<zxl.size();i++) {
    if ((zxl[i].zlist==1) && (zxl[i].on)) count++;
  }
  return count;
}

int String_Generator::X_Count() 
{
  int count = 0;
  for (long int i=1;i<zxl.size();i++) {
    if ((zxl[i].zlist==0) && (zxl[i].on)) count++;
  } 
  return count;
}

int String_Generator::E_Count() 
{
  int count = 0;
  for (long int i=1;i<zxl.size();i++) {
    if ((zxl[i].zlist==2) && (zxl[i].on)) count++;
  }
  return count;
}

int String_Generator::ZX_Count() 
{
  int count = 0;
  for (long int i=1;i<zxl.size();i++) {
    if (zxl[i].on) count++;
  } 
  return count;
}

int String_Generator::ZXY_Number(int type,int narg,int* arg,int ncoupl,int* coupl)
{
  for (long int i=1;i<zxl.size();i++) {
    if (zxl[i].zlist==type) {
      int hit = 1;
      for (short int j=0;j<narg;j++) {
	if (arg[j]!=zxl[i].arg[j]) {
	  hit = 0;
	  break;
	}
      }
      if (hit) {
	for (short int j=narg;j<narg+ncoupl;j++) {
	  if (coupl[j-narg]!=zxl[i].arg[j]) {
	    hit = 0;
	    break;
	  }
	}
      }
      if (hit) return i;
    }
  }
  return -1;
}

int String_Generator::Get_Cnumber(Complex coupl)
{
  for (short int i=0;i<couplings.size();i++) {
    if (AMATOOLS::IsEqual(coupl,couplings[i])) return i;
  }
  couplings.push_back(coupl);
  return couplings.size()-1;
}

int String_Generator::Get_Fnumber(int fl)
{
  for (short int i=0;i<flavours.size();i++) {
    if (fl==flavours[i]) return i;
  }
  flavours.push_back(fl);
  return flavours.size()-1;
}

Kabbala String_Generator::Number(int n,Complex value)
{
#ifdef Kabbala_on
  char help[10];
  sprintf(help,"Z[%i]",n);
  return Kabbala(string(help),value);
#else
  return value;
#endif      
}

int String_Generator::GetNumber(int type,Complex value)
{
  if (AMATOOLS::IsEqual(zxl[0].value.Value(),value)) return 0;

  for (long int i=1;i<zxl.size();i++) {
    if (zxl[i].zlist==type) {
      if (AMATOOLS::IsEqual(zxl[i].value.Value(),value)) return i;
    }
  }
  return zxl.size();
}


Kabbala String_Generator::Get_CZnumber(Complex value,string str)
{
  int numb = GetNumber(6,value);

  if (numb!=zxl.size()) return zxl[numb].value;

  ZXlist newz;

  newz.zlist = 6;

  String_Tree st;
  newz.sk = st.String2Tree(str);
  st.Delete(newz.sk,string("Z[0]"));
  
  if (newz.sk->op==0) {
    if (newz.sk->Str()==string("0")) return zxl[0].value;
  }

  st.DeleteMinus(newz.sk);

  String_Tree st2;
  newz.sk = st2.String2Tree(st.Tree2String(newz.sk,0));

  st2.Cluster(newz.sk,0);
  st2.DeleteMinus(newz.sk);
  //cout<<"Vor Delete: "<<st2.Tree2String(newz.sk,0)<<endl;

  st2.Delete(newz.sk,string("Z[0]"));

  if (newz.sk->op==0) {
    if (newz.sk->Str()==string("0")) return zxl[0].value;
  }
  /*
  list<sknot*> endpoint;
  stree.GetEnd(newz.sk,endpoint);
  for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it)
    (*it)->value = Get_Kabbala((*it)->Str());

  Complex vorher = stree.Evaluate(newz.sk);
  stree.Simplify(newz.sk);
  Complex nachher = stree.Evaluate(newz.sk);

  if (!AMATOOLS::IsEqual(vorher,nachher)) {
    cout<<"Error in String_Generator::Get_CZnumber(): "<<vorher<<";"<<nachher<<endl;
    //correcting
    newz.sk = st2.String2Tree(str);
    st2.Delete(newz.sk,string("Z[0]"));
  }
  */ 
  string newstr = st2.Tree2String(newz.sk,0);
  if ( newstr.find(string("+"))==-1 &&
       newstr.find(string("-"))==-1 &&
       newstr.find(string("*"))==-1 ) return Kabbala(newstr,value);

  //cout<<"Taken: "<<newstr<<endl;

  newz.sk = stree.String2Tree(newstr);
  newz.value = Number(zxl.size(),value);
  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Znumber(int* arg,Complex* coupl,Complex value) 
{
  int numb = GetNumber(1,value);
  if (numb!=zxl.size()) return zxl[numb].value;
  
  //new Zfunc  

  ZXlist newz;

  newz.zlist  = 1;
  newz.narg   = 12;
  newz.value  = Number(numb,value);
  newz.arg    = new int[12];
  for (short int i=0;i<8;i++)  newz.arg[i] = arg[i];
  for (short int i=8;i<12;i++) newz.arg[i] = Get_Cnumber(coupl[i-8]);

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Xnumber(int* arg,Complex* coupl,Complex value) 
{
  int numb = GetNumber(0,value);
  if (numb!=zxl.size()) return zxl[numb].value;

  //new Zfunc (X)
  ZXlist newz;

  newz.zlist = 0;
  newz.narg  = 7;
  newz.arg   = new int[7];

  newz.value = Number(numb,value);
  for (short int i=0;i<5;i++) newz.arg[i] = arg[i];
  for (short int i=5;i<7;i++) newz.arg[i] = Get_Cnumber(coupl[i-5]);

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Ynumber(int* arg,Complex* coupl,Complex value) 
{
  int numb = GetNumber(4,value);
  if (numb!=zxl.size()) return zxl[numb].value;

  //new Zfunc (Y)
  ZXlist newz;

  newz.zlist  = 4;
  newz.narg   = 6;
  newz.arg    = new int[6];
  newz.value = Number(numb,value);
  for (short int i=0;i<4;i++) newz.arg[i] = arg[i];
  for (short int i=4;i<6;i++) newz.arg[i] = Get_Cnumber(coupl[i-4]);

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Enumber(Complex value) 
{
  if (AMATOOLS::IsZero(value)) return zxl[0].value;

  int numb = GetNumber(2,value);
  if (numb!=zxl.size()) return zxl[numb].value;

  //new Zfunc E
  ZXlist newz;

  newz.zlist  = 2;
  newz.value  = Number(numb,value);

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Pnumber(Pfunc* pl,int numb) 
{
  for (long int i=0;i<zxl.size();i++) {
    if (zxl[i].zlist==5) {
      if( AMATOOLS::IsEqual(zxl[i].value.Value(),pl->value) &&
	  (flavours[int(zxl[i].arg[0])]==(pl->fl).Kfcode()) )
	return zxl[i].value;
    }
  }
  //new Zfunc P -> propagator

  ZXlist newz;

  newz.zlist  = 5;
  newz.value  = Number(zxl.size(),pl->value);
  newz.narg   = 2;
  newz.arg    = new int[2];
  newz.arg[0] = Get_Fnumber((pl->fl).Kfcode());
  newz.arg[1] = numb;

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Massnumber(int numb,Flavour fl,Complex value) 
{
  for (long int i=0;i<zxl.size();i++) {
    if (zxl[i].zlist==7) {
      if( AMATOOLS::IsEqual(zxl[i].value.Value(),value) &&
	  ( (flavours[int(zxl[i].arg[0])]==fl.Kfcode()) ||
	    (flavours[int(zxl[i].arg[0])]==-fl.Kfcode() && fl.IsAnti()) )
	  ) 
	return zxl[i].value;
    }
  }
  //new Zfunc P -> propagator
  ZXlist newz;

  newz.zlist = 7;
  newz.value = Number(zxl.size(),value);
  newz.narg  = 2;
  newz.arg   = new int[2];

  if (fl.IsAnti()) newz.arg[0] = Get_Fnumber(-fl.Kfcode());
              else newz.arg[0] = Get_Fnumber(fl.Kfcode());

  newz.arg[1] = numb;

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Snumber(int& a1,int& a2,Complex value) 
{
  for (long int i=0;i<zxl.size();i++) {
    if (zxl[i].zlist==3) {
      if ((zxl[i].arg[0]==a1) && (zxl[i].arg[1]==a2)) return zxl[i].value;
      if ((zxl[i].arg[1]==a1) && (zxl[i].arg[0]==a2)) return zxl[i].value;
    }
  }
  //new Zfunc S
  ZXlist newz;

  newz.zlist  = 3;
  newz.value  = Number(zxl.size(),value);
  newz.narg   = 2;
  newz.arg    = new int[2];
  newz.arg[0] = a1;
  newz.arg[1] = a2;

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_VNCnumber(int& a1,int& a2,Complex value) 
{
  for (long int i=0;i<zxl.size();i++) {
    if (zxl[i].zlist==3) {
      if ((zxl[i].arg[0]==a1) && (zxl[i].arg[1]==a2)) return zxl[i].value;
      if ((zxl[i].arg[1]==a1) && (zxl[i].arg[0]==a2)) return -zxl[i].value;
    }
  }
  //new Zfunc VNC
  ZXlist newz;

  newz.zlist  = 8;
  newz.value  = Number(zxl.size(),value);
  newz.narg   = 2;
  newz.arg    = new int[2];
  newz.arg[0] = a1;
  newz.arg[1] = a2;

  zxl.push_back(newz);

  return newz.value;
}

Kabbala String_Generator::Get_Scplxnumber(const int& a1,const int& a2,Complex value) 
{
  if (AMATOOLS::IsZero(value)) return zxl[0].value;
  for (long int i=0;i<zxl.size();i++) {
    if (zxl[i].zlist==9) {
      if ((zxl[i].arg[0]==a1) && (zxl[i].arg[1]==a2)) return zxl[i].value;
      if ((zxl[i].arg[1]==a1) && (zxl[i].arg[0]==a2)) return zxl[i].value;
    }
  }
  //AORGTOOLS::msg.Error()<<"*****Get_Scplxnumber: "<<a1<<","<<a2<<";"<<value<<endl;  //new Zfunc S
  ZXlist newz;

  newz.zlist  = 9;
  newz.value  = Number(zxl.size(),value);
  newz.narg   = 2;
  newz.arg    = new int[2];
  newz.arg[0] = a1;
  newz.arg[1] = a2;

  zxl.push_back(newz);

  return newz.value;
}

void String_Generator::Calculate(Values* val) 
{  
#ifdef Kabbala_on  
  if (val!=0) {
    val->Calculate(couplings);
    return;
  }

  for (long int i=1;i<zxl.size();i++) {
    if (zxl[i].on) {
      //cout<<"Calc Number "<<i<<";type = "<<zxl[i].zlist<<endl;
      //cout<<"    "<<zxl[i].value.String()<<" = "<<zxl[i].value.Value()<<" "<<endl;
      
      int* arg = zxl[i].arg;
      switch (zxl[i].zlist) {
      case 0: zxl[i].value = 
		Kabbala(zxl[i].value.String(),
			Xcalc(arg[0],arg[1],arg[2],arg[3],arg[4],
			      couplings[arg[5]],couplings[arg[6]]));
      break;
      case 1: zxl[i].value = 
		Kabbala(zxl[i].value.String(),
			Zcalc(arg[0],arg[1],arg[2],arg[3],arg[4],
			      arg[5],arg[6],arg[7],
			      couplings[arg[8]],couplings[arg[9]],
			      couplings[arg[10]],couplings[arg[11]]));
      break;
      case 2: break; //do nothing 				     
      case 3: zxl[i].value = 
		Kabbala(zxl[i].value.String(),Vcalc(arg[0],arg[1]));
      break;
      case 4: zxl[i].value = 
		Kabbala(zxl[i].value.String(),
			Ycalc(arg[0],arg[1],arg[2],arg[3],
			      couplings[arg[4]],couplings[arg[5]]));
      break;
      case 5: zxl[i].value = 
		Kabbala(zxl[i].value.String(), 
			Pcalc(flavours[arg[0]],arg[1]));
      break;      
      case 6: zxl[i].value = 
	  Kabbala(zxl[i].value.String(), 
		  stree.Evaluate(zxl[i].sk));
      break;	
      case 7: 
	zxl[i].value = 
		Kabbala(zxl[i].value.String(), 
			MassTermCalc(arg[1],flavours[arg[0]]));
      break;      
      case 8: zxl[i].value = 
		Kabbala(zxl[i].value.String(),VNCcalc(arg[0],arg[1]));
	break;
      case 9: zxl[i].value = 
 		Kabbala(zxl[i].value.String(),Vcplxcalc(arg[0],arg[1]));
      break;
      default:msg.Error()<<"Unknown Z-Type: "<<zxl[i].zlist<<endl;
      }
    }
  }
#endif
}

int String_Generator::Massless(int i)
{  
#ifdef Kabbala_on  
  int* arg = zxl[i].arg;
  switch (zxl[i].zlist) {
  case 1: return Zmassless(arg[0],arg[1],arg[2],arg[3],arg[4],
			   arg[5],arg[6],arg[7],
			   couplings[arg[8]],couplings[arg[9]],
			   couplings[arg[10]],couplings[arg[11]]);
  break;
  }
#endif
  return 0;
}

void String_Generator::SetOn(const string& str)
{
  if (str==string("0")) return;

  if (str==string("1") || str==string("1.") || str==string("0.5") || str==string("2") || str==string("2.")) return;

  string tststring = str.substr(2);
  tststring = tststring.substr(0,tststring.length()-1);
  MyStrStream   msstr;
  //  std::strstream msstr;  
  msstr<<tststring;
  int number;
  msstr>>number;
    
  if (zxl[number].value.String()==str) {
    zxl[number].on = 1;
    return;
  }

  msg.Error()<<"Error in String_Generator::SetOn()!"<<endl;
}

Kabbala* String_Generator::Get_Kabbala(const string& str)
{
#ifdef Kabbala_on  
  if (str==string("0")) return &zxl[0].value;

  if (str!=string("1") && str!=string("1.") && str!=string("0.5") && str!=string("2") && str!=string("2.")) {
    string tststring = str.substr(2);
    tststring = tststring.substr(0,tststring.length()-1);
    MyStrStream msstr;
    // std::stringstream msstr;
    msstr<<tststring;
    int number;
    msstr>>number;
    
    if (zxl[number].value.String()==str) {
      //zxl[number].on = 1;
      return &zxl[number].value; 
    }
    else {
      msg.Error()<<"Error in String_Generator::Get_Kabbala() : "<<str<<endl;
      abort();
    }
  }

  for (long int i=0;i<zxl.size();i++) {
    if (zxl[i].zlist==2) {
      if (str==string("1") && zxl[i].value.Value()==Complex(1.,0.)) {
	zxl[i].on = 1;
	return &zxl[i].value; 
      }
      if (str==string("1.") && zxl[i].value.Value()==Complex(1.,0.)) {
        zxl[i].on = 1;
        return &zxl[i].value;
      }
      if (str==string("0.5") && zxl[i].value.Value()==Complex(1./2.,0.)) {
	zxl[i].on = 1;
	return &zxl[i].value;
      } 
      if (str==string("2") && zxl[i].value.Value()==Complex(2.,0.)) {
	zxl[i].on = 1;
	return &zxl[i].value;
      } 
      if (str==string("2.") && zxl[i].value.Value()==Complex(2.,0.)) {
	zxl[i].on = 1;
	return &zxl[i].value;
      } 
    }
  }
  
  if (str==string("1") || str==string("1.")) {
    //new Zfunc E
    ZXlist newz;

    newz.zlist  = 2;
    newz.value  = Number(zxl.size(),Complex(1.,0.));
    newz.on     = 1;

    zxl.push_back(newz);

    return &zxl[zxl.size()-1].value;
  }

  if (str==string("2") || str==string("2.")) {
    //new Zfunc E
    ZXlist newz;

    newz.zlist  = 2;
    newz.value  = Number(zxl.size(),Complex(2.,0.));
    newz.on     = 1;

    zxl.push_back(newz);

    return &zxl[zxl.size()-1].value;
  }

  if (str==string("0.5")) {
    //new Zfunc E
    ZXlist newz;

    newz.zlist  = 2;
    newz.value  = Number(zxl.size(),Complex(1./2.,0.));
    newz.on     = 1;

    zxl.push_back(newz);

    return &zxl[zxl.size()-1].value;
  }
#endif
  msg.Error()<<"Error: No Zvalue for String: "<<str<<endl;
  return 0;
}



