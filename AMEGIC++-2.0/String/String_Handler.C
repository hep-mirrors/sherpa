#include "String_Handler.H"
#include "String_Output.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "prof.hh"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

String_Handler::String_Handler(const int &_gen_str,const string& _path,Basic_Sfuncs* _BS) 
  : gen_str(_gen_str), path(_path)
{
  working = 0;
  sk      = 0;
  val     = 0;
  own_sgen= 1;
  if (gen_str==0) sgen = new No_String_Generator;
  else {
    sgen = new String_Generator(_BS);

    string pID = path;
    for (short int i=path.length()-1;i>=0;i--) {
      if (path[i]=='/') {
	pID    = string("V")+path.substr(i+1);
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
    if (gen_str==2) val = Set_Values(pID,_BS);
    if (val!=0) {
      msg.Debugging()<<pID<<" loaded."<<endl; 
      val->SetCouplFlav();
      working = 1;
    }
  }
  if (val==0) msg.Tracking()<<"No Value-Library available !"<<endl;
}



String_Handler::String_Handler(const int &_gen_str,Basic_Sfuncs* _BS) 
  : gen_str(_gen_str)
{
  own_sgen= 1;
  working = 0;
  sk      = 0;
  val     = 0;
  if (gen_str==0) sgen = new No_String_Generator;
             else sgen = new String_Generator(_BS);
}


String_Handler::String_Handler(Virtual_String_Generator * _sgen) 
{
  own_sgen= 0;
  working = 0;
  sk      = 0;
  val     = 0;

  sgen    = _sgen;
}

bool String_Handler::SearchValues(const int _gen_str,string & pID,Basic_Sfuncs* _BS) 
{
  string vpID = string("V")+pID;
  if (_gen_str==2) val = Set_Values(vpID,_BS);
  if (val!=0) {
    msg.Tracking()<<vpID<<" loaded."<<endl;
    val->SetCouplFlav();
    working = 1;
    return 1;
  }
  else {
    msg.Tracking()<<vpID<<" not loaded."<<endl; 
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
  if (sk!=0) {
    for (short int i=0;i<maxgraph;i++) delete[] sk[i];
    delete[] sk;
  }
  if (own_sgen)  // delete only if constructed
    delete sgen;
}

void String_Handler::Set_String(const int &igraph,const int &ihel,
				const string& str)
{
  if (gen_str==0 || val!=0 || working==1) return;
  //sk[igraph][ihel] = stree.String2Tree(str);
  //stree.Delete(sk[igraph][ihel],string("Z[0]"));
  String_Tree st;
  sknot* shelp = st.String2Tree(str);
  st.Delete(shelp,string("Z[0]"));

  stringsk[igraph][ihel] = st.Tree2String(shelp,0);
  //sk[igraph][ihel] = stree.String2Tree(stringsk[igraph][ihel]);
}

sknot* String_Handler::MakeSknotAFill(string & str)
{
  String_Tree st;
  
  sknot* shelp = st.String2Tree(str);
  
  st.DeleteMinus(shelp);

  String_Tree st2;
  shelp = st2.String2Tree(st.Tree2String(shelp,0));

  st2.Cluster(shelp,0);
  st2.DeleteMinus(shelp);
  st2.Delete(shelp,string("Z[0]"));
  
  shelp = stree.Copy(shelp,1);

  list<sknot*> endpoint;
  stree.GetEnd(shelp,endpoint);
  for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it) { 
    (*it)->value = sgen->GetKabbala((*it)->Str());
    sgen->SetOn((*it)->Str());	
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
  msg.Debugging()<<"Completing the strings, this may take some time...."<<endl;

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
  int countall = 0;
  int counton = 0;
  
  for (long int i=1;i<sgen->ZXMaxNumber();i++) {
    countall++;
    if (sgen->GetZXl(i)->on) counton++;
  }

  for (long int i=sgen->ZXMaxNumber()-1;i>0;i--) {
    if (sgen->GetZXl(i)->zlist==6 && sgen->GetZXl(i)->on) {
      if (sgen->GetZXl(i)->sk!=0) {
	endpoint.clear();
	stree.GetEnd(sgen->GetZXl(i)->sk,endpoint);
	for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it) 
	  (*it)->value = sgen->GetKabbala((*it)->Str());
	/*
	string spre = stree.Tree2String(sgen->GetZXl(i)->sk,0);
	stree.Simplify(sgen->GetZXl(i)->sk);
	string safter = stree.Tree2String(sgen->GetZXl(i)->sk,0);
	if (spre!=safter) {
	  //changed endpoints
	  endpoint.clear();
	  stree.GetEnd(sgen->GetZXl(i)->sk,endpoint);
	}
	*/
	for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it) 
	  sgen->SetOn((*it)->Str());
      }
    }
  }

  countall = 0;
  counton = 0;
  
  for (long int i=1;i<sgen->ZXMaxNumber();i++) {
    countall++;
    if (sgen->GetZXl(i)->on) counton++;
  }

  Z_Kill();  
}

Complex String_Handler::Evaluate(int igraph,int ihel)
{
  if (val!=0) return val->Evaluate(igraph,ihel);
  if (sk[igraph][ihel]==0) return Complex(0.,0.);
  return stree.Evaluate(sk[igraph][ihel]); 
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
  String_Output so(path,maxgraph,maxhel);
  so.Output(sk,&stree,sgen,hel);
}

void String_Handler::Output(Helicity* hel, string path)
{
  if (gen_str<2 || val!=0) return;
  String_Output so(path,maxgraph,maxhel);
  so.Output(sk,&stree,sgen,hel);
}

void String_Handler::Z_Kill()
{
  int count = 0;
  msg.Debugging()<<"Number of Z functions: "<<sgen->ZXMaxNumber()<<endl;
  string str;
  for (long int i=1;i<sgen->ZXMaxNumber();i++) {
    if (sgen->GetZXl(i)->on==0) {
      //search in 6
      str = sgen->GetZXl(i)->value.String();
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
             else count++;
    }
  }
  
  msg.Debugging()<<count<<"/"<<sgen->ZXMaxNumber()<<" Z functions have been deleted."<<endl;
}

