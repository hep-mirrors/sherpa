#include <iostream>
#include "AMEGIC++/String/String_Generator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;
#ifdef __GNUC__
#include <stdio.h>
#endif


ZXlist& ZXlist::operator=(const ZXlist& copy) {
  if (this!=&copy) {
    narg  = copy.narg;
    value = copy.value;
    zlist = copy.zlist;
    on    = copy.on;
    sk    = copy.sk;
    if (arg!=0) delete[] arg;
    arg = new int[narg];
    for (short int i=0;i<narg;i++) arg[i] = copy.arg[i];
  }
  return *this;
}

String_Generator::String_Generator(Basic_Sfuncs* _BS) : 
  Basic_Func(this,_BS), Basic_Yfunc(this,_BS), Basic_Zfunc(this,_BS), 
  Basic_Xfunc(this,_BS), Basic_Mfunc(this,_BS), Basic_Vfunc(this,_BS), 
  Basic_Pfunc(this,_BS), Basic_MassTermfunc(this,_BS), 
  Basic_Epsilonfunc(this,_BS), Unitarityfunc(this,_BS),
  p_zxlsave(NULL), p_couplingssave(NULL), p_flavourssave(NULL), m_copied(0)
{ 
  m_zuse.resize(13);
  for (size_t i=0;i<13;i++) m_zuse[i]=0;
  p_zxl       = new vector<ZXlist>;
  p_couplings = new vector<Complex>;
  p_flavours  = new vector<long int>;
  Reset();
}

String_Generator::~String_Generator() 
{
  if (p_zxl)           { delete p_zxl;           p_zxl           = NULL; }
  if (p_couplings)     { delete p_couplings;     p_couplings     = NULL; }
  if (p_flavours)      { delete p_flavours;      p_flavours      = NULL; }
  if (p_zxlsave)       { delete p_zxlsave;       p_zxlsave       = NULL; }
  if (p_couplingssave) { delete p_couplingssave; p_couplingssave = NULL; }
  if (p_flavourssave)  { delete p_flavourssave;  p_flavourssave  = NULL; }
}

void String_Generator::Reset(int f)
{
  if (f) (*p_couplings).clear();
  (*p_flavours).clear();
  (*p_zxl).clear();

  ZXlist zero;
  zero.value = Kabbala(string("Z[0]"),Complex(0.,0.));
  zero.zlist = -1;

  (*p_zxl).push_back(zero);
}

void String_Generator::ReplaceZXlist(Virtual_String_Generator* _sgen)
{
  if (!m_copied) {
    p_zxlsave       = p_zxl;
    //p_couplingssave = p_couplings;
    p_flavourssave  = p_flavours;
  }
  p_zxl           = _sgen->GetZXlist();
  //p_couplings     = coupl;
  p_flavours      = _sgen->GetFlavours();

  m_copied = 1;
}

void String_Generator::ReStore()
{
  if (m_copied) {
    p_zxl           = p_zxlsave;
    //p_couplings     = p_couplingssave;
    p_flavours      = p_flavourssave;
    p_zxlsave       = NULL;
    //p_couplingssave = NULL;
    p_flavourssave  = NULL;

    m_copied        = 0;
    return;
  }
  msg_Error()<<"Error in String_Generator::ReStore():"<<endl
	     <<"   No information stored to be restored ! Abort."<<endl;
  abort();
}


void String_Generator::Print()
{
  if (!(msg_LevelIsDebugging())) return;
  for (size_t i=0;i<(*p_zxl).size();i++) {
    msg_Out()<<i<<". Zfunction: Type="<<(*p_zxl)[i].zlist<<";On="<<(*p_zxl)[i].on<<";Value="<<(*p_zxl)[i].value.String(); 
    if ((*p_zxl)[i].narg>0) msg_Out()<<";Arg[0] = "<<(*p_zxl)[i].arg[0];
    msg_Out()<<endl;
  }
}

int String_Generator::ZCount() 
{
  int count = 0;
  for (size_t i=1;i<(*p_zxl).size();i++) {
    if (((*p_zxl)[i].zlist==1) && ((*p_zxl)[i].on)) count++;
  }
  return count;
}

int String_Generator::XCount() 
{
  int count = 0;
  for (size_t i=1;i<(*p_zxl).size();i++) {
    if (((*p_zxl)[i].zlist==0) && ((*p_zxl)[i].on)) count++;
  } 
  return count;
}

int String_Generator::ECount() 
{
  int count = 0;
  for (size_t i=1;i<(*p_zxl).size();i++) {
    if (((*p_zxl)[i].zlist==2) && ((*p_zxl)[i].on)) count++;
  }
  return count;
}

int String_Generator::ZXCount() 
{
  int count = 0;
  for (size_t i=1;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].on) count++;
  } 
  return count;
}

int String_Generator::ZXYNumber(int type,int narg,int* arg,int ncoupl,int* coupl)
{
  for (size_t i=1;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==type) {
      int hit = 1;
      for (int j=0;j<narg;j++) {
	if (arg[j]!=(*p_zxl)[i].arg[j]) {
	  hit = 0;
	  break;
	}
      }
      if (hit) {
	for (short int j=narg;j<narg+ncoupl;j++) {
	  if (coupl[j-narg]!=(*p_zxl)[i].arg[j]) {
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

int String_Generator::GetCnumber(Complex coupl)
{
  for (size_t i=0;i<(*p_couplings).size();i++) {
    if (ATOOLS::IsEqual(coupl,(*p_couplings)[i])) return i;
  }
  (*p_couplings).push_back(coupl);
  return (*p_couplings).size()-1;
}

int String_Generator::GetFnumber(long int fl)
{
  for (size_t i=0;i<(*p_flavours).size();i++) {
    if (fl==(*p_flavours)[i]) return i;
  }
  (*p_flavours).push_back(fl);
  return (*p_flavours).size()-1;
}

Kabbala String_Generator::Number(int n,Complex value)
{
  char help[10];
  sprintf(help,"Z[%i]",n);
  return Kabbala(string(help),value);
}

int String_Generator::GetNumber(int type,Complex value)
{
  if (ATOOLS::IsEqual((*p_zxl)[0].value.Value(),value)) return 0;

  for (size_t i=1;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==type) {
      if (ATOOLS::IsEqual((*p_zxl)[i].value.Value(),value)) return i;
    }
  }
  return (*p_zxl).size();
}


Kabbala String_Generator::GetCZnumber(Complex value,string str)
{
  int numb = GetNumber(6,value);

  if (numb!=(int)(*p_zxl).size()) return (*p_zxl)[numb].value;

  ZXlist newz;

  newz.zlist = 6;

  m_sthelp.Reset();
  newz.sk = m_sthelp.String2Tree(str);
  m_sthelp.Delete(newz.sk,string("Z[0]"));
  
  if (newz.sk->op==0) {
    if (newz.sk->Str()==string("0")) return (*p_zxl)[0].value;
  }
  m_sthelp.DeleteMinus(newz.sk);

  if (newz.sk->op==0) {
    if (newz.sk->Str()==string("0")) return (*p_zxl)[0].value;
  }
  string newstr = m_sthelp.Tree2String(newz.sk,0);
  if ( newstr.find(string("+"))==string::npos &&
       newstr.find(string("-"))==string::npos &&
       newstr.find(string("*"))==string::npos ) return Kabbala(newstr,value);

  newz.sk = m_sthelp.String2Tree(newstr);
  m_sthelp.DeleteMinus(newz.sk);
  newz.sk = m_stree.Copy(newz.sk,0);
  newz.value = Number((*p_zxl).size(),value);
  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetZnumber(int* arg,Complex* coupl,Complex value) 
{
  int numb = GetNumber(1,value);
  if (numb!=(int)(*p_zxl).size()) return (*p_zxl)[numb].value;
  
  //new Zfunc  

  ZXlist newz;

  newz.zlist  = 1;
  newz.narg   = 12;
  newz.value  = Number(numb,value);
  newz.arg    = new int[12];
  for (short int i=0;i<8;i++)  newz.arg[i] = arg[i];
  for (short int i=8;i<12;i++) newz.arg[i] = GetCnumber(coupl[i-8]);

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetXnumber(int* arg,Complex* coupl,Complex value) 
{
  int numb = GetNumber(0,value);
  if (numb!=(int)(*p_zxl).size()) return (*p_zxl)[numb].value;

  //new Zfunc (X)
  ZXlist newz;

  newz.zlist = 0;
  newz.narg  = 7;
  newz.arg   = new int[7];

  newz.value = Number(numb,value);
  for (short int i=0;i<5;i++) newz.arg[i] = arg[i];
  for (short int i=5;i<7;i++) newz.arg[i] = GetCnumber(coupl[i-5]);

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetYnumber(int* arg,Complex* coupl,Complex value) 
{
  int numb = GetNumber(4,value);
  if (numb!=(int)(*p_zxl).size()) return (*p_zxl)[numb].value;

  //new Zfunc (Y)
  ZXlist newz;

  newz.zlist  = 4;
  newz.narg   = 6;
  newz.arg    = new int[6];
  newz.value = Number(numb,value);
  for (short int i=0;i<4;i++) newz.arg[i] = arg[i];
  for (short int i=4;i<6;i++) newz.arg[i] = GetCnumber(coupl[i-4]);

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetEnumber(Complex value) 
{
  if (value==Complex(0.,0.)) return (*p_zxl)[0].value;

  int numb = GetNumber(2,value);
  if (numb!=(int)(*p_zxl).size()) return (*p_zxl)[numb].value;

  //new Zfunc E
  ZXlist newz;

  newz.zlist  = 2;
  newz.value  = Number(numb,value);

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetPnumber(Pfunc* pl,int numb) 
{
  for (size_t i=0;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==5) {
      if( ATOOLS::IsEqual((*p_zxl)[i].value.Value(),pl->value) &&
	  ((*p_flavours)[int((*p_zxl)[i].arg[0])]==(long int)(pl->fl).Kfcode()) )
	return (*p_zxl)[i].value;
    }
  }
  //new Zfunc P -> propagator

  ZXlist newz;

  newz.zlist  = 5;
  newz.value  = Number((*p_zxl).size(),pl->value);
  newz.narg   = 2;
  newz.arg    = new int[2];
  newz.arg[0] = GetFnumber((pl->fl).Kfcode());
  newz.arg[1] = numb;
  
  (*p_zxl).push_back(newz);
  
  return newz.value;
}

Kabbala String_Generator::GetMassnumber(int numb,ATOOLS::Flavour fl,Complex value) 
{
  for (size_t i=0;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==7) {
      if( ATOOLS::IsEqual((*p_zxl)[i].value.Value(),value) &&
	  ( (*p_flavours)[int((*p_zxl)[i].arg[0])]==(long int)(fl) )
	  ) 
	return (*p_zxl)[i].value;
    }
  }
  //new Zfunc P -> propagator
  ZXlist newz;

  newz.zlist = 7;
  newz.value = Number((*p_zxl).size(),value);
  newz.narg  = 2;
  newz.arg   = new int[2];

  if (fl.IsAnti()) newz.arg[0] = GetFnumber(-fl.Kfcode());
              else newz.arg[0] = GetFnumber(fl.Kfcode());

  newz.arg[1] = numb;

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetMnumber(ATOOLS::Flavour fl,Complex value) 
{
  for (size_t i=0;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==8) {
      if (ATOOLS::IsEqual((*p_zxl)[i].value.Value(),value)) return (*p_zxl)[i].value;
    }
  }
  //new Zfunc P -> propagator
  ZXlist newz;

  newz.zlist = 8;
  newz.value = Number((*p_zxl).size(),value);
  newz.narg  = 1;
  newz.arg   = new int[1];

  newz.arg[0] = fl.Kfcode();

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetSnumber(const int a1,const int a2,Complex value) 
{
  for (size_t i=0;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==3) {
      if (((*p_zxl)[i].arg[0]==a1) && ((*p_zxl)[i].arg[1]==a2)) return (*p_zxl)[i].value;
      if (((*p_zxl)[i].arg[1]==a1) && ((*p_zxl)[i].arg[0]==a2)) return (*p_zxl)[i].value;
    }
  }
  //new Zfunc S
  ZXlist newz;

  newz.zlist  = 3;
  newz.value  = Number((*p_zxl).size(),value);
  newz.narg   = 2;
  newz.arg    = new int[2];
  newz.arg[0] = a1;
  newz.arg[1] = a2;

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetScplxnumber(const int a1,const int a2,Complex value) 
{
  if (value==Complex(0.0,0.0)) return (*p_zxl)[0].value;
  for (size_t i=0;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==9) {
      if (((*p_zxl)[i].arg[0]==a1) && ((*p_zxl)[i].arg[1]==a2)) return (*p_zxl)[i].value;
      if (((*p_zxl)[i].arg[1]==a1) && ((*p_zxl)[i].arg[0]==a2)) return (*p_zxl)[i].value;
    }
  }
  ZXlist newz;

  newz.zlist  = 9;
  newz.value  = Number((*p_zxl).size(),value);
  newz.narg   = 2;
  newz.arg    = new int[2];
  newz.arg[0] = a1;
  newz.arg[1] = a2;

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetEpsnumber(int* arg,int nc,Complex value) 
{
  int numb = GetNumber(10,value);
  if (numb!=(int)(*p_zxl).size()) return (*p_zxl)[numb].value;
  numb = GetNumber(10,-value);
  if (numb!=(int)(*p_zxl).size()) return -(*p_zxl)[numb].value;
  
  //new Epsfunc  

  ZXlist newz;

  newz.zlist  = 10;
  newz.narg   = 5;
  newz.value  = Number(numb,value);
  newz.arg    = new int[5];
  for (int i=0;i<4;i++)  newz.arg[i] = arg[i];
  newz.arg[4] = nc;

  (*p_zxl).push_back(newz);

  return newz.value;
}

Kabbala String_Generator::GetSFnumber(Complex value,int zn) 
{
  for (size_t i=0;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==zn) {
      return (*p_zxl)[i].value;
    }
  }
  //new special function
  ZXlist newz;

  newz.zlist  = zn;
  newz.value  = Number((*p_zxl).size(),value);
  newz.narg   = 0;
  newz.arg    = 0;
  (*p_zxl).push_back(newz);

  return newz.value;
}

void String_Generator::Calculate(Values* val) 
{  
  if (val!=0) {
    val->Calculate();
    return;
  }

  for (size_t i=1;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].on) {
      int* arg = (*p_zxl)[i].arg;
      if ((*p_zxl)[i].zlist>=0&&(*p_zxl)[i].zlist<13) m_zuse[(*p_zxl)[i].zlist]=1;
      switch ((*p_zxl)[i].zlist) {
      case 0: (*p_zxl)[i].value = 
		Kabbala((*p_zxl)[i].value.String(),
			Xcalc(arg[0],arg[1],arg[2],arg[3],arg[4],
			      (*p_couplings)[arg[5]],(*p_couplings)[arg[6]]));
      break;
      case 1: (*p_zxl)[i].value = 
		Kabbala((*p_zxl)[i].value.String(),
			Zcalc(arg[0],arg[1],arg[2],arg[3],arg[4],
			      arg[5],arg[6],arg[7],
			      (*p_couplings)[arg[8]],(*p_couplings)[arg[9]],
			      (*p_couplings)[arg[10]],(*p_couplings)[arg[11]]));
      break;
      case 2: break; //do nothing 				     
      case 3: (*p_zxl)[i].value = 
		Kabbala((*p_zxl)[i].value.String(),Vcalc(arg[0],arg[1]));
      break;
      case 4: (*p_zxl)[i].value = 
		Kabbala((*p_zxl)[i].value.String(),
			Ycalc(arg[0],arg[1],arg[2],arg[3],
			      (*p_couplings)[arg[4]],(*p_couplings)[arg[5]]));
      break;
      case 5: (*p_zxl)[i].value = 
		Kabbala((*p_zxl)[i].value.String(), 
			Pcalc((*p_flavours)[arg[0]],arg[1]));
      break;      
      case 6: (*p_zxl)[i].value = 
	  Kabbala((*p_zxl)[i].value.String(), 
		  m_stree.Evaluate((*p_zxl)[i].sk));
      break;	
      case 7: 
	(*p_zxl)[i].value = 
		Kabbala((*p_zxl)[i].value.String(), 
			MassTermCalc(arg[1],(*p_flavours)[arg[0]]));
      break;      
      case 8: 
	(*p_zxl)[i].value = 
		Kabbala((*p_zxl)[i].value.String(), 
			1./Complex(sqr(Flavour(arg[0]).Mass()),
				   -Flavour(arg[0]).Mass()*Flavour(arg[0]).Width()));
      break;      
      case 9: (*p_zxl)[i].value = 
 		Kabbala((*p_zxl)[i].value.String(),Vcplxcalc(arg[0],arg[1]));
      break;
      case 10: (*p_zxl)[i].value = 
 		Kabbala((*p_zxl)[i].value.String(),EpsCalc(arg[0],arg[1],arg[2],arg[3],arg[4]));
      break;
      case 11: (*p_zxl)[i].value = 
 		Kabbala((*p_zxl)[i].value.String(),Ucalc(3));
      break;
      case 12: (*p_zxl)[i].value = 
 		Kabbala((*p_zxl)[i].value.String(),Ucalc(4));
      break;
      default:msg_Error()<<"Unknown Z-Type: "<<(*p_zxl)[i].zlist<<endl;
      }
    }
  }
}

int String_Generator::Massless(int i)
{  
  int* arg = (*p_zxl)[i].arg;
  switch ((*p_zxl)[i].zlist) {
  case 1: return Zmassless(arg[0],arg[1],arg[2],arg[3],arg[4],
			   arg[5],arg[6],arg[7],
			   (*p_couplings)[arg[8]],(*p_couplings)[arg[9]],
			   (*p_couplings)[arg[10]],(*p_couplings)[arg[11]]);
  break;
  }
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
    
  if ((*p_zxl)[number].value.String()==str) {
    (*p_zxl)[number].on = 1;
    return;
  }

  msg_Error()<<"Error in String_Generator::SetOn()!"<<endl;
}

Kabbala* String_Generator::GetKabbala(const string& str)
{
  if (str==string("0")) return &(*p_zxl)[0].value;

  if (str!=string("1") && str!=string("1.") && str!=string("0.5") && str!=string("2") && str!=string("2.")) {
    string tststring = str.substr(2);
    tststring = tststring.substr(0,tststring.length()-1);
    MyStrStream msstr;
    msstr<<tststring;
    int number;
    msstr>>number;
    
    if ((*p_zxl)[number].value.String()==str) {
      (*p_zxl)[number].on = 1;
      return &(*p_zxl)[number].value; 
    }
    else {
      msg_Error()<<"Error in String_Generator::GetKabbala() : "<<str<<endl;
      abort();
    }
  }

  for (size_t i=0;i<(*p_zxl).size();i++) {
    if ((*p_zxl)[i].zlist==2) {
      if (str==string("1") && (*p_zxl)[i].value.Value()==Complex(1.,0.)) {
	(*p_zxl)[i].on = 1;
	return &(*p_zxl)[i].value; 
      }
      if (str==string("1.") && (*p_zxl)[i].value.Value()==Complex(1.,0.)) {
        (*p_zxl)[i].on = 1;
        return &(*p_zxl)[i].value;
      }
      if (str==string("0.5") && (*p_zxl)[i].value.Value()==Complex(1./2.,0.)) {
	(*p_zxl)[i].on = 1;
	return &(*p_zxl)[i].value;
      } 
      if (str==string("2") && (*p_zxl)[i].value.Value()==Complex(2.,0.)) {
	(*p_zxl)[i].on = 1;
	return &(*p_zxl)[i].value;
      } 
      if (str==string("2.") && (*p_zxl)[i].value.Value()==Complex(2.,0.)) {
	(*p_zxl)[i].on = 1;
	return &(*p_zxl)[i].value;
      } 
    }
  }
  
  if (str==string("1") || str==string("1.")) {
    //new Zfunc E
    ZXlist newz;

    newz.zlist  = 2;
    newz.value  = Number((*p_zxl).size(),Complex(1.,0.));
    newz.on     = 1;

    (*p_zxl).push_back(newz);

    return &(*p_zxl)[(*p_zxl).size()-1].value;
  }

  if (str==string("2") || str==string("2.")) {
    //new Zfunc E
    ZXlist newz;

    newz.zlist  = 2;
    newz.value  = Number((*p_zxl).size(),Complex(2.,0.));
    newz.on     = 1;

    (*p_zxl).push_back(newz);

    return &(*p_zxl)[(*p_zxl).size()-1].value;
  }

  if (str==string("0.5")) {
    //new Zfunc E
    ZXlist newz;

    newz.zlist  = 2;
    newz.value  = Number((*p_zxl).size(),Complex(1./2.,0.));
    newz.on     = 1;

    (*p_zxl).push_back(newz);

    return &(*p_zxl)[(*p_zxl).size()-1].value;
  }
  msg_Error()<<"Error: No Zvalue for String: "<<str<<endl;
  return 0;
}


void String_Generator::WriteCouplings(ofstream& os)
{
  os.precision(16);
  for (int i=0;i<NumberOfCouplings();i++) {
    int cnt=0;
    for (map<string,Complex>::iterator mit=p_couplmap->begin();mit!=p_couplmap->end();mit++) {
      if (mit->second==GetCoupling(i)) {
	cnt++;
	os<<"cpl["<<i<<"]="<<mit->first<<endl;
      }
    }
    if (cnt==0) {
      os<<"cpl["<<i<<"]=Complex"<<GetCoupling(i)<<endl;
    }
  }
}
 
int String_Generator::ReadCouplings(ifstream& is)
{
  is.precision(16);
  if (!is) return 0;
  p_couplings->clear();
  string str;
  for (;is;) {
    getline(is,str); 
    if (str.find(string("cpl"))==0) {
      int a=str.find("[");
      int b=str.find("]");
      int idx = ToType<int>(str.substr(a+1,a-b-1));
      a = str.find("=");
      Complex val;
      if (str.substr(a+1,7)==string("Complex")) {
	val = ToType<Complex>(str.substr(a+8));
      }
      else {
	str = str.substr(a+1);
	if (p_couplmap->find(str)==p_couplmap->end()) 
	  THROW(critical_error,"Coupling constant not found, check MODEL settings!");
	val = (*p_couplmap)[str];
      }
      if (idx>=NumberOfCouplings()) p_couplings->resize(idx+1);
      if (GetCoupling(idx)==Complex(0.,0.)) (*p_couplings)[idx] = val;
      else if (GetCoupling(idx)!=val) {
	THROW(critical_error,"Coupling constant inconsistency, lib incompatible with MODEL settings!");
      }
    }
  }
  return NumberOfCouplings();
}

void String_Generator::UpdateCouplings(map<string,Complex> & cmap)
{
  size_t cnt=0;
  for (int i=0;i<NumberOfCouplings();i++) {
    string help = ToString(GetCoupling(i));
    if (cmap.find(help)!=cmap.end()) {
      (*p_couplings)[i] = cmap[help];
      cnt++;
    }
  }
  if (cnt!=cmap.size()) {
    THROW(critical_error,"String_Generator::UpdateCouplings() failed: Coupling constant inconsistency.");    
  }
}
