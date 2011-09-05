#include "AMEGIC++/String/String_Handler.H"
#include "AMEGIC++/String/String_Output.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;


String_Handler::String_Handler(const int &_gen_str,Basic_Sfuncs* _BS,std::map<std::string,Complex>* cplmap) 
  : gen_str(_gen_str)
{
  m_mode=1;

  own_sgen= 1;
  own_sk  = 1;
  working = 0;
  sk      = 0;
  stringsk= 0;
  val     = 0;
  if (gen_str==0) sgen = new No_String_Generator;
             else sgen = new String_Generator(_BS);
  sgen->SetCouplings(cplmap);
}


String_Handler::String_Handler(Virtual_String_Generator * _sgen) 
{
  m_mode=1;

  own_sgen= 0;
  own_sk  = 1;
  working = 0;
  sk      = 0;
  stringsk= 0;
  val     = 0;

  sgen    = _sgen;
}

String_Handler::String_Handler(Virtual_String_Generator * _sgen, sknot *** _sk) 
{
  m_mode=1;

  own_sgen= 0;
  own_sk  = 0;
  working = 1;
  sk      = _sk;
  stringsk= 0;
  val     = 0;

  sgen    = _sgen;
}

bool String_Handler::SearchValues(const int _gen_str,string & pID,Basic_Sfuncs* _BS) 
{
  string vpID = string("V")+pID;
  if (_gen_str==2) {
    if (val) delete val;
    val = Set_Values(vpID,_BS);
  }
  if (val!=0) {
    msg_Info()<<"."<<std::flush;
    val->SetCouplFlav(*sgen->GetCouplList());
    if (sgen->NumberOfCouplings()!=val->NumberOfCouplings()) {
      msg_Error()<<" Number of Coupling constants does not fit with Process Library "<<pID<<"!"<<endl;
      return 0;
    }
    working = 1;
    return 1;
  }
  else {
    return 0;
  }
}


void String_Handler::Initialize(const int& _maxgraph,const int& _maxhel)
{
  if (gen_str==0) return;
  maxgraph = _maxgraph;
  maxhel   = _maxhel;
  
  if (val==0) {
    stringsk = new string*[maxgraph];
    sk       = new sknot**[maxgraph];
    for (short int i=0;i<maxgraph;i++) {
      sk[i]       = new sknot*[maxhel];
      stringsk[i] = new string[maxhel];
      for (short int j=0;j<maxhel;j++) {
	sk[i][j]       = 0;
	stringsk[i][j] = string("");
      }
    }
  }
}

String_Handler::~String_Handler()
{
  if (sk!=0 && own_sk) {
    for (short int i=0;i<maxgraph;i++) delete[] sk[i];
    delete[] sk;
  }
  if (own_sgen)  // delete only if constructed
    delete sgen;
  KillStrings();
  if (val) {
    delete val;
  }
}

void String_Handler::KillStrings()
{
  if (stringsk) {
    for (short int i=0;i<maxgraph;i++) delete[] stringsk[i];
    delete[] stringsk;
  }
  stringsk = NULL;
}

void String_Handler::Set_String(const int &igraph,const int &ihel,
				const string& str)
{
  if (gen_str==0 || val!=0 || working==1) return;
  sthelp.Reset();
  sknot* shelp = sthelp.String2Tree(str);
  sthelp.Delete(shelp,string("Z[0]"));
  sthelp.DeleteMinus(shelp);

  stringsk[igraph][ihel] = sthelp.Tree2String(shelp,0);
}

sknot* String_Handler::MakeSknotAFill(string & str)
{
  sthelp.Reset();
  sknot* shelp = sthelp.String2Tree(str);
  sthelp.DeleteMinus(shelp);
  shelp = stree.Copy(shelp,1);

  list<sknot*> endpoint;
  stree.GetEnd(shelp,endpoint);
  for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it) { 
    (*it)->value = sgen->GetKabbala((*it)->Str());
  }

  return shelp;
}

void String_Handler::Complete(Helicity* hel)
{
  if (gen_str==0) {
    sgen->Reset();
    return;
  }
  working = 1;

  if (val!=0) return;
  msg_Info()<<"In String_Handler::Complete : this may take some time...."<<endl;

  //connect sgenZ to treeZ
  list<sknot*> endpoint;

  for (long int i=1;i<sgen->ZXMaxNumber();i++) sgen->SetOff(i);

  for (short int j=0;j<maxhel;j++) {
    for (short int i=0;i<maxgraph;i++) {
      if (stringsk[i][j].length()>0 && hel->On(j)) {
	if (hel->On(j)) sk[i][j] = MakeSknotAFill(stringsk[i][j]);
	           else sk[i][j] = 0;
	//delete stringsk
	stringsk[i][j] = string("");
      }
    }
  }

  for (long int i=sgen->ZXMaxNumber()-1;i>0;i--) {
    if (sgen->GetZXl(i)->zlist==6 && sgen->GetZXl(i)->on) {
      if (sgen->GetZXl(i)->sk!=0) {
	endpoint.clear();
	stree.GetEnd(sgen->GetZXl(i)->sk,endpoint);
	for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it) 
	  (*it)->value = sgen->GetKabbala((*it)->Str());
       }
    }
  }

  Z_Kill();  
}

Complex String_Handler::Zvalue(int igraph,int ihel) 
{
  if (val!=0) return val->Evaluate(igraph,ihel);
  if (sk[igraph][ihel]==0) return Complex(0.,0.);
  return stree.Evaluate(sk[igraph][ihel]); 
}

void String_Handler::Calculate() {
  //  PROFILE_HERE;
  sgen->Calculate(val);
}

void String_Handler::Output(Helicity* hel)
{
  if (gen_str<2 || val!=0) return;
  String_Output so(path,maxgraph,maxhel,m_mode);
  so.Output(sk,&stree,sgen,hel);
  KillStrings();
}

void String_Handler::Output(Helicity* hel, string path)
{
  if (gen_str<2 || val!=0) return;
  String_Output so(path,maxgraph,maxhel,m_mode);
  so.Output(sk,&stree,sgen,hel);
  KillStrings();
}

void String_Handler::Z_Kill()
{
  int count = 0;
  string str;
  for (long int i=1;i<sgen->ZXMaxNumber();i++) {
    if (sgen->GetZXl(i)->on==0) {
      //search in 6
      /*str = sgen->GetZXl(i)->value.String();
      int hit = 0;
      for (long int j=i+1;j<sgen->ZXMaxNumber();j++) {
	if (sgen->GetZXl(j)->zlist==6 && sgen->GetZXl(j)->on) {
	  if (sgen->GetZXl(j)->sk!=0) {
	    stree.Find(sgen->GetZXl(j)->sk,str,hit);
	    if (hit==1) break;
	  }
	}
      }
      if (hit==1) sgen->SetOn(i);
      else*/ 
      count++;
    }
  }
}

typedef Values* (*Lib_Getter_Function)(Basic_Sfuncs*);

Values* String_Handler::Set_Values(std::string& pID,Basic_Sfuncs* BS)
{
  s_loader->AddPath(rpa->gen.Variable("SHERPA_LIB_PATH"));
  Lib_Getter_Function gf = (Lib_Getter_Function)s_loader->GetLibraryFunction
    ("Proc_"+pID.substr(1),"Getter_"+pID);
  if (gf==NULL) return NULL;
  return gf(BS);
}
