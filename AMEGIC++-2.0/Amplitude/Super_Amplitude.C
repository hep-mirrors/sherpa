#include "Super_Amplitude.H"
#include "Zfunc.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;


void Super_Amplitude::Init(string _str)
{

  str = _str;

  for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) {  
    int old = zlist->size();

    Zfunc_List* gzlist = (*g)->GetZlist();
    for (Zfunc_Iterator zit=gzlist->begin();zit!=gzlist->end();++zit) {
      int hit = 0;
      for (Zfunc_Iterator zit2=zlist->begin();zit2!=zlist->end();++zit2) {
	if ((*zit)->p_equal==(*zit2)->p_equal) {
	  hit = 1;
	  break;
	}
      }
      if (hit==0) {
	zlist->push_back((*zit));
	//zlist->push_back(new Zfunc(*(*zit)));
      }
    } 

    Pfunc_List* gplist = (*g)->GetPlist(); 
    
    vector<Pair> pairlist;

    for (Pfunc_Iterator pit=gplist->begin();pit!=gplist->end();++pit) {
      Pfunc* p = *pit;
      //search for number in old list
      
      Pfunc* pequalnum  = 0;
      Pfunc* pequalprop = 0;
      
      for (Pfunc_Iterator pit2=plist.begin();pit2!=plist.end();++pit2) {
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

    int i=0;
    for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit,++i) {
      if (i>=old) {
	(*zit)->ReplaceProp(&pairlist);	
      }
    }
  }

  ReduceZfuncs(str);
 
  sign = graphs.front()->GetSign();

  int hit = 0;

  for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) {  
    if (sign!=(*g)->GetSign()) {
      hit = 1;
    }
  }
  if (hit) SetZfuncSign(); 
}

int Super_Amplitude::NewSigns(vector<vector<int> > & zsignlists) {
  int hit =0;

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

  return hit;
}


void Super_Amplitude::SetZfuncSign()
{
  //search for suitable factor

  vector<vector<int> > zsignlists;

  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {

    // if super zf (+)
      if ((*zit)->GetOp()=='+')
	  zsignlists.push_back(vector<int>((*zit)->GetSize(),1));
      else
	  zsignlists.push_back(vector<int>(1,1));
      // if pregroup zf (*)
      // if single zf
  }

  int global_sign=GetSign();

  int ok;
  for (;;) {
    ok=1;
    // try
    for (vector<Amplitude_Base*>::iterator g=graphs.begin();g!=graphs.end();++g) { 
      int sign=1;
      Zfunc_List* gzlist = (*g)->GetZlist(); 
      for (Zfunc_Iterator gzit=gzlist->begin();gzit!=gzlist->end();++gzit) {
	// looking for zfunc in superampl
	int i =0;
	for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit,++i) {
	    int zsize=(*zit)->GetSize();
	    if (!((*zit)->GetOp()=='+'))
		zsize=1;
	  for (int j=0;j<zsize;j++) {
	    Zfunc* z = (*(*zit))[j]; 
	    if (zsize==1) z=(*zit);  
	    if ((*gzit)->m_str==z->m_str) {
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
    }
    if (ok) {
      int i =0;
      for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit,++i) {
	for (int j=0;j<(*zit)->GetSize();++j) {
	    if ((*zit)->GetOp()=='+') { 
	      (*zit)->SetSign(j,zsignlists[i][j]);
	  }
	}
      }
      break;
    }
    if (!NewSigns(zsignlists)) break;
  }
  // did not find a permutation
  if (ok==0) {
    msg.Error()<<"Found no suitable factor in Super_Amplitude::SetZfuncSign()!"<<endl;
    abort();
  }
}

void Super_Amplitude::ReduceZfuncs(string str)
{
  String_Tree st;
  sknot* shelp = st.String2Tree(str);
  
  list<sknot*> factorlist;
  st.Factors(shelp,factorlist);

  for (list<sknot*>::iterator fit=factorlist.begin();fit!=factorlist.end();++fit) {
    list<sknot*> zfunclist;
    st.Addends(*fit,zfunclist);
    Zfunc_Group* superfunc;
    
    int first = 1;
   
    for (list<sknot*>::iterator sit=zfunclist.begin();sit!=zfunclist.end();++sit) {
      int hit = 0;
      for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
	if ((*zit)->m_str==st.Tree2String(*sit,0)) {
	  hit = 1;
	  if (first) {
	    first = 0;
	    superfunc = new Zfunc_Group(*(*zit));
	    superfunc->m_str = st.Tree2String(*fit,0);
	  }
	  superfunc->m_zlist.push_back(*zit);
	  superfunc->m_zsigns.push_back(1);
	  zlist->erase(zit);
	  break;
	}
      }
      if (hit==0) {
	cerr<<"No Zfunc found in Super_Amplitude::ReduceZfuncs()!"<<endl;
	abort();
      }
    }
    if(superfunc->GetSize()==1){
      zlist->push_back((*superfunc)[0]);
      delete superfunc;
    }
    else zlist->push_back(superfunc);
  }  
}


int Super_Amplitude::FindNewNumber(int number)
{
  int hit;
  do {
    hit = 0;
    for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
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
  if (!rpa.gen.Tracking()) return;
  msg.Out()<<"--------"<<amplnumber+1<<". Amplitude----------"<<endl;
  Single_Amplitude_Base::PrintGraph();

  msg.Out()<<"Overall sign "<<sign<<endl;
}

Complex Super_Amplitude::Zvalue(int ihel,int * signlist)   
  {return Single_Amplitude_Base::Zvalue(ihel,signlist);}

Complex Super_Amplitude::Zvalue(String_Handler * sh,int ihel)   
  {return Single_Amplitude_Base::Zvalue(sh,ihel);}
