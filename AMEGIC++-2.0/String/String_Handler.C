#include "String_Handler.H"
#include "String_Output.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "prof.hh"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

String_Handler::String_Handler(const int &_gen_str,const string& _path,Basic_Sfuncs* _BS) 
  : gen_str(_gen_str), path(_path)
{
  
  //testing
  /*
  String_Tree st;
  sknot* shelp = st.String2Tree("A*(B+C)-A*A*(D+E)");

  st.DeleteMinus(shelp);
  cout<<"Vor Cluster: "<<st.Tree2String(shelp,0)<<endl;
  st.Cluster(shelp,0);
  st.DeleteMinus(shelp);
  st.Delete(shelp,string("Z[0]"));
  cout<<"Nachher: "<<st.Tree2String(shelp,0)<<endl;
  abort();
  */


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
      msg.Tracking()<<pID<<" loaded."<<endl; 
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
    sk = new sknot**[maxgraph];
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
    (*it)->value = sgen->Get_Kabbala((*it)->Str());
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
  msg.Tracking()<<"Completing the strings, this may take some time...."<<endl;

  //connect sgenZ to treeZ
  list<sknot*> endpoint;

  for (long int i=1;i<sgen->ZX_Max_Number();i++) sgen->SetOff(i);

  for (short int j=0;j<maxhel;j++) {
    //    cout<<"Helicity: "<<j<<"; Knotsize = "<<stree.SknotListSize()<<endl;
    for (short int i=0;i<maxgraph;i++) {
      if (stringsk[i][j].length()>0 && hel->On(j)) {
	if (hel->On(j)) sk[i][j] = MakeSknotAFill(stringsk[i][j]);
	           else sk[i][j] = 0;
	//delete stringsk
	stringsk[i][j] = string("");
	/*
	//Simplify
	string spre = stree.Tree2String(sk[i][j],0);
	stree.Simplify(sk[i][j]);
	string safter = stree.Tree2String(sk[i][j],0);
	if (spre!=safter) {
	  //changed endpoints
	  endpoint.clear();
	  stree.GetEnd(sk[i][j],endpoint);
	}
	*/
      }
    }
  }
  int countall = 0;
  int counton = 0;
  
  for (long int i=1;i<sgen->ZX_Max_Number();i++) {
    countall++;
    if (sgen->Get_ZXl(i)->on) counton++;
  }

  msg.Tracking()<<counton<<"/"<<countall<<" direct Z-Functions."<<endl;


  for (long int i=sgen->ZX_Max_Number()-1;i>0;i--) {
    if (sgen->Get_ZXl(i)->zlist==6 && sgen->Get_ZXl(i)->on) {
      if (sgen->Get_ZXl(i)->sk!=0) {
	endpoint.clear();
	stree.GetEnd(sgen->Get_ZXl(i)->sk,endpoint);
	for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it) 
	  (*it)->value = sgen->Get_Kabbala((*it)->Str());
	/*
	string spre = stree.Tree2String(sgen->Get_ZXl(i)->sk,0);
	stree.Simplify(sgen->Get_ZXl(i)->sk);
	string safter = stree.Tree2String(sgen->Get_ZXl(i)->sk,0);
	if (spre!=safter) {
	  //changed endpoints
	  endpoint.clear();
	  stree.GetEnd(sgen->Get_ZXl(i)->sk,endpoint);
	}
	*/
	for (list<sknot*>::iterator it=endpoint.begin();it!=endpoint.end();++it) 
	  sgen->SetOn((*it)->Str());
      }
    }
  }

  countall = 0;
  counton = 0;
  
  for (long int i=1;i<sgen->ZX_Max_Number();i++) {
    countall++;
    if (sgen->Get_ZXl(i)->on) counton++;
  }

  msg.Tracking()<<counton<<"/"<<countall<<" direct+indirect Z-Functions."<<endl;

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
  //msg.Out()<<igraph<<";"<<ihel<<" : "<<stree.Tree2String(sk[igraph][ihel],0)<<endl;
  return stree.Evaluate(sk[igraph][ihel]); 
}

void String_Handler::Calculate() {
  //  PROFILE_HERE;
  sgen->Calculate(val);
}

void String_Handler::Output(Helicity* hel)
{
#ifdef Kabbala_on  
  if (gen_str<2 || val!=0) return;
  msg.Debugging()<<"String_Handler::Output() : "<<path<<endl;
  String_Output so(path,maxgraph,maxhel);
  so.Output(sk,&stree,sgen,hel);
#endif
}

void String_Handler::Output(Helicity* hel, string path)
{
#ifdef Kabbala_on  
  if (gen_str<2 || val!=0) return;
  msg.Debugging()<<"String_Handler::Output() : "<<path<<endl;
  String_Output so(path,maxgraph,maxhel);
  so.Output(sk,&stree,sgen,hel);
#endif
}

void String_Handler::Z_Kill()
{
#ifdef Kabbala_on
  int count = 0;
  msg.Tracking()<<"Number of Z functions: "<<sgen->ZX_Max_Number()<<endl;
  string str;
  for (long int i=1;i<sgen->ZX_Max_Number();i++) {
    if (sgen->Get_ZXl(i)->on==0) {
      //search in 6
      str = sgen->Get_ZXl(i)->value.String();
      int hit = 0;
      for (long int j=i+1;j<sgen->ZX_Max_Number();j++) {
	if (sgen->Get_ZXl(j)->zlist==6 && sgen->Get_ZXl(j)->on) {
	  if (sgen->Get_ZXl(j)->sk!=0) {
	    stree.Find(sgen->Get_ZXl(j)->sk,str,hit);
	    if (hit==1) break;
	  }
	}
      }
      if (hit==1) sgen->SetOn(i);
             else count++;
    }
  }
  
  msg.Tracking()<<count<<"/"<<sgen->ZX_Max_Number()<<" Z functions have been deleted."<<endl;
#endif
}

