#include "Super_Amplitude.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

class Pair {
public:
  int pold,pnew;
  
  Pair(int _pold,int _pnew) : pold(_pold), pnew(_pnew) {}
};



void Super_Amplitude::Init(string _str)
{

  str = _str;

  for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) {  
    int old = zlist.size();

    list<Zfunc*>* gzlist = (*g)->GetZlist(); 
    for (list<Zfunc*>::iterator zit=gzlist->begin();zit!=gzlist->end();++zit) {
      int hit = 0;
      for (list<Zfunc*>::iterator zit2=zlist.begin();zit2!=zlist.end();++zit2) {
	if ((*zit)->Equal==(*zit2)->Equal) {
	  hit = 1;
	  break;
	}
      }
      if (hit==0) zlist.push_back(new Zfunc(*(*zit)));
    }


    list<Pfunc*>* gplist = (*g)->GetPlist(); 
    
    vector<Pair> pairlist;

    for (list<Pfunc*>::iterator pit=gplist->begin();pit!=gplist->end();++pit) {
      Pfunc* p = *pit;
      //search for number in old list
      
      Pfunc* pequalnum  = 0;
      Pfunc* pequalprop = 0;
      
      for (list<Pfunc*>::iterator pit2=plist.begin();pit2!=plist.end();++pit2) {
	Pfunc* p2 = *pit2;
	if (p->momnum==p2->momnum && p->fl==p2->fl) pequalprop = p2;
	if (p->arg[0]==p2->arg[0])                  pequalnum  = p2;
      }
      
      if (pequalnum==pequalprop) {
	if (pequalnum==0) plist.push_back(new Pfunc(*p));  //totally new propagator
      }
      else {
	if (pequalprop==0) { //new propagator with already used number
	  //found equal prop with different number
	  int newnumb = FindNewNumber(p->arg[0]);
	  pairlist.push_back(Pair(p->arg[0],newnumb));
	  plist.push_back(new Pfunc(*p));
	  plist.back()->arg[0] = newnumb;
	}
	else { //propagator already in the list with non-equal number
	  //found equal prop
	  pairlist.push_back(Pair(p->arg[0],pequalprop->arg[0]));
	}
      }
    }
    //Print all Pairs
    /*
    for (int i=0;i<pairlist.size();i++) 
      cout<<"Change "<<pairlist[i].pold<<" --> "<<pairlist[i].pnew<<endl;
    cout<<"================================"<<endl;
    */
    //cout<<"Old is: "<<old<<"/"<<

    int i=0;
    for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit,++i) {
      if (i>=old) {
	Zfunc* z = *zit;
	for (int j=0;j<z->narg;j++) {
	  for (int k=0;k<pairlist.size();k++) {
	    if (pairlist[k].pold==z->arg[j]) {
	      z->arg[j] = pairlist[k].pnew;
	      break;
	    }
	  }
	}
	//cout<<"Pn: "<<z->pn<<endl;
	for (int j=0;j<z->pn;j++) {
	  //cout<<"Prop: "<<z->psnew[j].numb<<endl;
	  for (int k=0;k<pairlist.size();k++) {
	    if (pairlist[k].pold==z->psnew[j].numb) {
	      //cout<<"From "<<z->psnew[j].numb<<" to "<<pairlist[k].pnew<<endl;
	      z->psnew[j].numb = pairlist[k].pnew;
	      break;
	    }
	  }
	}
      }
    }
  }

  ReduceZfuncs(str);
  /*  
  cout<<"New Superamplitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  Amplitude_Group::PrintGraph();
  Single_Amplitude_Base::PrintGraph();
  abort();
  */

  sign = graphs.front()->GetSign();

  int hit = 0;

  for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) {  
    if (rpa.gen.Tracking()) {
      list<Zfunc*>* gzlist = (*g)->GetZlist(); 
      cout<<(*g)->GetSign()<<flush;
      for (list<Zfunc*>::iterator gzit=gzlist->begin();gzit!=gzlist->end();++gzit) {
	cout<<" "<<(*gzit)->str<<flush;
      }
      cout<<endl;
    }
    if (sign!=(*g)->GetSign()) {
      cout<<" Different Signs in the sub amplitudes!"<<endl;
      hit = 1;
      // *AS*      break;
    }
  }
  if (hit) SetZfuncSign(); 
}

int Super_Amplitude::NewSigns(vector<vector<int> > & zsignlists) {
  int hit =0;

  //  cout<<" in Super_Amplitude::NewSigns "<<endl;

  // don't change the first sign of each super zfunc
  for (int i=zsignlists.size()-1;i>=0;--i) {
    for (int j=zsignlists[i].size()-1;j>0;--j) {    
      if (zsignlists[i][j]==1) {
	zsignlists[i][j]=-1; 
	hit = 1;
	break;
      }
      else {
	zsignlists[i][j]=1; 
      }
    }
    if (hit) break;
  }

  /*
  cout <<" new Signlist: "<<endl;
  for (int i=0;i<zsignlists.size();++i) { 
    for (int j=0;j<zsignlists[i].size();++j) {    
      cout<<" "<<zsignlists[i][j];
    }
    cout<<endl;
  }
  */

  return hit;
}


void Super_Amplitude::SetZfuncSign()
{
  //  cout<<"In SetZfuncSign"<<endl;
  //search for suitable factor

  vector<vector<int> > zsignlists;

  for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
    zsignlists.push_back(vector<int>((*zit)->GetSize(),1));
  }

  int global_sign=GetSign();

  int ok;
  for (;;) {
    ok=1;
    // probieren
    for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) { 
      int sign=1;
      list<Zfunc*>* gzlist = (*g)->GetZlist(); 
      for (list<Zfunc*>::iterator gzit=gzlist->begin();gzit!=gzlist->end();++gzit) {
	// looking for zfunc in superampl
	int i =0;
	for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit,++i) {
	  for (int j=0;j<(*zit)->GetSize();j++) {
	    Zfunc* z = (*(*zit))[j];
	    if ((*gzit)->str==z->str) {
	      //gotcha
	      sign*=zsignlists[i][j];
	    }
	  }
	}
      }
      if (sign != global_sign * ((*g)->GetSign())) {
	ok = 0;
	break;
      }
      //      cout<<"  Global Sign = "<<global_sign<<endl;
      /*
      int i =0;
      for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit,++i) {
	for (int j=0;j<(*zit)->GetSize();++j) {
	  cout<<"   "<<(*(*zit))[j]->str<<" has sign "<<zsignlists[i][j]<<endl;
	}
      }
      */

    }
    if (ok) {
      // found permutation
      cout<<"Found a suitable permutation!!"<<endl;
      //      cout<<"  Global Sign = "<<global_sign<<endl;
      int i =0;
      for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit,++i) {
	for (int j=0;j<(*zit)->GetSize();++j) {
	  (*zit)->SetSign(j,zsignlists[i][j]);
	  //	  cout<<"   "<<(*(*zit))[j]->str<<" has sign "<<zsignlists[i][j]<<endl;
	}
      }
      break;
    }
    if (!NewSigns(zsignlists)) break;
  }
  // did not find a permutation
  if (ok==0) {
    cerr<<"Found no suitable factor in Super_Amplitude::SetZfuncSign()!"<<endl;
    abort();
  }
}

/*
void Super_Amplitude::SetZfuncSign()
{
  cout<<"In SetZfuncSign"<<endl;
  //search for suitable factor
  int ok = 0;
  for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
    if ((*zit)->GetSize()>1) {
      ok = 1;
      vector<int> zsignlist;
      for (int i=0;i<(*zit)->GetSize();i++) {
	Zfunc* z = (*(*zit))[i];
	int zsign = 0;
	for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) { 
	  list<Zfunc*>* gzlist = (*g)->GetZlist(); 
	  for (list<Zfunc*>::iterator gzit=gzlist->begin();gzit!=gzlist->end();++gzit) {
	    if ((*gzit)->str==z->str) {
	      //gotcha
	      if (zsign==0) zsign = (*g)->GetSign();
	      if (zsign!=(*g)->GetSign()) ok = 0;
	      break;
	    }
	  }
	  if (ok==0) break;
	}
	if (ok==0) break;
	      else zsignlist.push_back(zsign);
      }
      if (ok==1) {
	cout<<"Found a suitable factor!!"<<endl;
	//found a factor 
	for (int i=0;i<(*zit)->GetSize();i++) (*zit)->SetSign(i,sign*zsignlist[i]);
	break;
      }
    }
  }
  if (ok==0) {
    cerr<<"Found no suitable factor in Super_Amplitude::SetZfuncSign()!"<<endl;
    abort();
  }
}
*/

void Super_Amplitude::ReduceZfuncs(string str)
{
  String_Tree st;
  sknot* shelp = st.String2Tree(str);
  
  //cout<<"String: "<<st.Tree2String(shelp,0)<<endl;

  list<sknot*> factorlist;
  st.Factors(shelp,factorlist);

  //cout<<"Factorlist-Size: "<<factorlist.size()<<endl;

  for (list<sknot*>::iterator fit=factorlist.begin();fit!=factorlist.end();++fit) {
    //cout<<"   New factor: "<<st.Tree2String(*fit,0)<<endl;

    list<sknot*> zfunclist;
    st.Addends(*fit,zfunclist);
    //new Superfunc
    Zfunc_Group* superfunc;
    
    int first = 1;
   
    for (list<sknot*>::iterator sit=zfunclist.begin();sit!=zfunclist.end();++sit) {
      //cout<<"   New Zfunc: "<<st.Tree2String(*sit,0)<<endl; 
      int hit = 0;
      for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
	if ((*zit)->str==st.Tree2String(*sit,0)) {
	  hit = 1;
	  //cout<<"New Hit for "<<(*zit)->str<<endl;
	  if (first) {
	    first = 0;
	    superfunc = new Zfunc_Group(*(*zit));
	    superfunc->str = st.Tree2String(*fit,0);
	  }
	  superfunc->zlist.push_back(*zit);
	  superfunc->zsign.push_back(1);
	  zlist.erase(zit);
	  break;
	}
      }
      if (hit==0) {
	cerr<<"No Zfunc found in Super_Amplitude::ReduceZfuncs()!"<<endl;
	abort();
      }
    }
    zlist.push_back(superfunc);
  }  
}


int Super_Amplitude::FindNewNumber(int number)
{
  int hit;
  do {
    hit = 0;
    for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
      Pfunc* p = *pit;
      if (p->arg[0]==number) {
	hit = 1;
	break;
      }
    }
    if (hit) number++;
  }
  while (hit>0);

  return number;
}
 
void Super_Amplitude::PrintGraph() 
{
  msg.Out()<<"--------"<<amplnumber+1<<". Amplitude----------"<<endl;
  Single_Amplitude_Base::PrintGraph();

  msg.Out()<<"Overall sign "<<sign<<endl;
}

Kabbala Super_Amplitude::Single_Zvalue(Argument* args,Zfunc* z)
{ 
  Kabbala value;
  
  for (int i=0;i<z->GetSize();i++) {
    if (z->GetSign(i)==-1) value -= Single_Amplitude_Base::Single_Zvalue(arglist[i],(*z)[i]);
                      else value += Single_Amplitude_Base::Single_Zvalue(arglist[i],(*z)[i]);
  }
  
  if (buildstring && z->GetSize()>1) {
    value = (shand->Get_Generator())->Get_CZnumber(value.Value(),value.String());
  }

  return value;
}

void Super_Amplitude::Fill_Args(Zfunc* z,Argument* args,int* signlist,
				vector<int>* iz,vector<int>* ii,
				vector<vector<int> >* iargs) 
{
  for (int i=0;i<arglist.size();i++) delete[] arglist[i];
  arglist.clear();

  for (int i=0;i<z->GetSize();i++) {
    Single_Amplitude_Base::Fill_Args((*z)[i],args,signlist,iz,ii,iargs);
    Argument* saveargs = new Argument[2*z->narg];
    for (int j=0;j<2*z->narg;j++) {
      saveargs[j] = args[j]; 
      args[j].spinortype = Spinor::None;
    }
    arglist.push_back(saveargs);
  }
}


Complex Super_Amplitude::Zvalue(int ihel,int * signlist)   
  {return Single_Amplitude_Base::Zvalue(ihel,signlist);}
//{return Amplitude_Group::Zvalue(ihel,signlist);}

Complex Super_Amplitude::Zvalue(String_Handler * sh,int ihel)   
  {return Single_Amplitude_Base::Zvalue(sh,ihel);}
