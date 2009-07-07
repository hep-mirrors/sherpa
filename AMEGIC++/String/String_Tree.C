#include "AMEGIC++/String/String_Tree.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include <cstring>

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

string sknot::emptystring = string("");

sknot String_Tree::zero;
const int String_Tree::block_size = 256;

String_Tree::String_Tree() 
{
  CSC.Init();
  zero.op    = 0;
  zero.SetString(string("0"));
  zero.left  = 0;
  zero.right = 0;
  skpos       =-1;
}

String_Tree::~String_Tree() 
{
  for (vector<sknot*>::iterator it=sblocks.begin();it!=sblocks.end();++it) 
    delete[] (*it);
}

sknot* String_Tree::newsk() 
{ 
  ++skpos;
  if (skpos%block_size==0) {
    if (skpos/block_size==(int)sblocks.size()) {
      sblocks.push_back(new sknot[block_size]);
    }    
  }

  return &sblocks[skpos/block_size][skpos%block_size];
}

void String_Tree::popsk()
{
  sblocks[skpos/block_size][skpos%block_size].KillString();
  --skpos;
}

void String_Tree::Reset() 
{
  for (int i=0;i<=skpos;i++) sblocks[i/block_size][i%block_size].KillString();
  skpos=-1;
}

void String_Tree::CleanValues()
{
  for (int i=0;i<=skpos;i++) sblocks[i/block_size][i%block_size].KillValue();
}

sknot* String_Tree::String2Tree(string term,int fixed) 
{
  if (term==string("")) return 0;

  sknot* m = newsk();

  long int kl;
  long int i;
  int sw1;

  for (;;) {
    sw1 = 1;
    //plus search
    kl = 0;
    for (i=0;i<(int)term.length();i++) {
      switch (term[i]) {
      case '(':kl++;break;
      case ')':kl--;break;
      case '+':if (kl==0) {
	if (i!=0) {
	  m->op = '+';
	  m->left  = String2Tree(term.substr(0,i),fixed);
	  m->right = String2Tree(term.substr(i+1),fixed);
	  return m;
	}
      }
      break;			
      case '-':if (kl==0) {
	if (i!=0) {
	  m->op = '+';
	  m->left  = String2Tree(term.substr(0,i),fixed);
	  m->right = String2Tree(term.substr(i),fixed);
	  return m;
	}
      }
      break;
      }
    }
    //times search
    kl = 0;
    for (i=0;i<(int)term.length();i++) {
      switch (term[i]) {
      case '(':kl++;break;
      case ')':kl--;break;
      case '*':
	if (kl==0) {
	  m->op = '*';
	  m->right = String2Tree(term.substr(0,i),fixed);
	  m->left  = String2Tree(term.substr(i+1),fixed);
	  //change *1
	  if (m->right->op==0) {
	    if (m->right->Str()==string("1") || m->right->Str()==string("1.")) {
	      m = m->left;
	      return m;
	    }
	  }
	  if (m->left->op==0) {
	    if (m->left->Str()==string("1") || m->left->Str()==string("1.")) {
	      m = m->right;
	      return m;
	    }
	  } 
	  
	  return m;
	}	
	break;
      }
    }
    //divided search !!!!!!!!!!!!! (a/b/c goes wrong) 
    kl = 0;
    for (i=0;i<(int)term.length();i++) {
      switch (term[i]) {
      case '(':kl++;break;
      case ')':kl--;break;
      case '/':if (kl==0) {
	m->op = '/';
	m->left = String2Tree(term.substr(0,i),fixed);
	m->right  = String2Tree(term.substr(i+1),fixed);
	return m;
      }
      break;
      }
    }
 
    //strip ()
    if ((term[0]=='(') && (term[term.length()-1]==')')) {
      sw1  = 0;
      term = term.substr(1,term.length()-2);
    }
    else {
      if (term[0]=='-') {
	if (term.substr(1)==string("") || term.substr(1)==string("0")) {
	  term = string("0");
	  return Leaf(term,m,fixed);
	}
	m->op = '-';
	m->left  = String2Tree(string("0"),fixed);
	m->right = String2Tree(term.substr(1),fixed);
	return m;
      }
      if (term[0]=='+') {
	sw1 = 0;
	term = term.substr(1);
      }
      if (term.substr(0,3)==string("exp")) {
	m->op = 'e';
	m->left  = String2Tree(string("0"),fixed);
	m->right = String2Tree(term.substr(4,term.length()-5),fixed);
	return m;
      }	
      if (term.substr(0,3)==string("sin")) {
	m->op = 's';
	m->left  = String2Tree(string("0"),fixed);
	m->right = String2Tree(term.substr(4,term.length()-5),fixed);
	return m;
      }	
    }
    
    if (sw1) break;
  }
   
  return Leaf(term,m,fixed);
}

sknot* String_Tree::Leaf(string& term,sknot* m,int fixed)
{
  if (fixed) {
    // search for Leaf
    for (list<sknot*>::iterator it=leaflist.begin();it!=leaflist.end();++it) {
      if ((*it)->op==0) {
	if ((*it)->Str()==term) {
	  //kill the last in the list (since it is equal m)
	  popsk();
	  return (*it); 
	}
      }
    }
  }
  if (term==string("1.")) term = string("1");
  m->op    = 0;
  m->SetString(term);
  m->left  = 0;
  m->right = 0;

  //add to leaflist
  if (fixed) leaflist.push_back(m);

  return m;
}


Complex String_Tree::eval(sknot* m)
{
  switch (m->op) {
  case '+': return eval(m->left)+eval(m->right);
  case '-': return eval(m->left)-eval(m->right);
  case '*': return eval(m->left)*eval(m->right);
  default: 
    if (m->Str()==string("CF"))  return Complex(4./3.,0.);
    if (m->Str()==string("i"))   return Complex(0.,1.);
    if (m->Str()==string("2"))   return Complex(2.,0.);
    if (m->Str()==string("1"))   return Complex(1.,0.);
    if (m->Str()==string("3"))   return Complex(3.,0.);
    if (m->Str()==string("8"))   return Complex(8.,0.);
    if (m->Str()==string("1/2")) return Complex(1./2.,0.);
    if (m->Str()==string("1/3")) return Complex(1./3.,0.);
    if (m->Str()==string("0.5")) return Complex(1./2.,0.);
    if (m->Str()==string("0.33")) return Complex(1./3.,0.);
    return 0.;
  }
}


Complex String_Tree::evalcolor(sknot* m)
{
  switch (m->op) {
  case '+': return evalcolor(m->left)+evalcolor(m->right);
  case '-': return evalcolor(m->left)-evalcolor(m->right);
  case '*': return evalcolor(m->left)*evalcolor(m->right);
  default: 
    if (m->Str()==string("CF"))  return Complex(CSC.CF,0.);
    if (m->Str()==string("i"))   return Complex(0.,1.);
    if (m->Str()==string("2"))   return Complex(2.,0.);
    if (m->Str()==string("1"))   return Complex(1.,0.);
    if (m->Str()==string("1/2")) return Complex(1./2.,0.);
    if (m->Str()==string("0.5")) return Complex(1./2.,0.);
    if (m->Str()==string("Nc"))  return Complex(CSC.Nc,0.);
    if (m->Str()==string("iNc")) return Complex(1./CSC.Nc,0.);
    if (m->Str()==string("Ng"))  return Complex(CSC.Nc*CSC.Nc-1.,0.);
    if (m->Str()!="0") cout<<METHOD<<" : possible evaluation error! "<<m->Str()<<endl;
    return 0.;
  }
}


void String_Tree::Find(sknot* m,const string&str,int& hit)
{
  if (hit==1) return;
  if (m->op==0) {
    if (str==m->Str()) hit = 1;
    return;
  }
  Find(m->left,str,hit);
  Find(m->right,str,hit);
}

void String_Tree::Addends(sknot* m,list<sknot*> &addend_list)
{
  if (m==0) return;

  if (m->op=='+' || m->op=='-') {
    Addends(m->left,addend_list);
    Addends(m->right,addend_list);
    return;
  }

  addend_list.push_back(m);
}

void String_Tree::Factors(sknot* m,list<sknot*> &factor_list)
{
  if (m==0) return;

  if (m->op=='*') {
    Factors(m->left,factor_list);
    Factors(m->right,factor_list);
    return;
  }

  factor_list.push_back(m);
}

void String_Tree::GetEnd(sknot* m,list<sknot*> &endpoint)
{
  if (m==0) return;
  if (m->op==0) {
    endpoint.push_back(m);
    return;
  }
  
  GetEnd(m->left,endpoint);
  GetEnd(m->right,endpoint);
}

Complex String_Tree::Evaluate(sknot* m)
{
  //if (m==0) return 0;
  switch (m->op) {
  case '+': return Evaluate(m->left)+Evaluate(m->right);
  case '-': return Evaluate(m->left)-Evaluate(m->right);
  case '*': return Evaluate(m->left)*Evaluate(m->right);
#ifndef __GNUC__
  case 'e': return std::exp(Evaluate(m->right));
  case 's': return std::sin(Evaluate(m->right));
#else
  case 'e': return std::exp(Evaluate(m->right));
  case 's': return std::sin(Evaluate(m->right));
#endif
  default:
    if (m->Str()==string("0")) return Complex(0.,0.);
    return m->value->Value();
  }
}


string String_Tree::Tree2String(sknot* m,sknot* g)
{
  if (m->op!=0) {
    char cc[2];
    string sop;
    cc[1] = 0;cc[0]=m->op;
    sop = string(cc);

    if (g==0 || m->op=='*') {   
      if (m->op=='*') {
	if (m->left->Str()==string("1")) return Tree2String(m->right,m);
	if (m->right->Str()==string("1")) return Tree2String(m->left,m);
	if (m->left->Str()==string("0") || m->right->Str()==string("0")) 
	  return string("0");
      }
      return Tree2String(m->left,m)+sop+Tree2String(m->right,m);
    }

    if (m->right->Str()==string("0")) {
      if (m->left->op==0) return Tree2String(m->left,m);
      sop = string("");
    }

    if (g->op=='+') {
      if (m->op=='-' && m->left->op==0) {
	if (m->left->Str()==string("0")) 
	  return string("(")+sop
	    +Tree2String(m->right,m)+string(")");
      }
      else 
	return Tree2String(m->left,m)+sop+Tree2String(m->right,m);
    }

    if (g->op=='-') {
      if (m==g->left) return Tree2String(m->left,m)+sop+Tree2String(m->right,m);
      if (m==g->right) return string("(")+Tree2String(m->left,m)+sop
			                 +Tree2String(m->right,m)+string(")");
   
    }
    if (m->op=='+' || m->op=='-') {
      return string("(")+Tree2String(m->left,m)+sop
	                +Tree2String(m->right,m)+string(")");
    }
    if (m->op=='e') 
      return string("exp(") + Tree2String(m->right,m) + string(")");
    if (m->op=='s') 
      return string("sin(") + Tree2String(m->right,m) + string(")");

  }

  if (m->left!=0 || m->right!=0) {
    msg_Error()<<"Error in Tree2String: "<<m->Str()<<endl;
    abort();
  }

  if (m->Str()==string("0")) return string("");

  if (m->Str()==string("2")) return string("2."); 
 
  if (m->Str()==string("1")) return string("1."); 
 
  return m->Str();
}

string String_Tree::Tree2Tex(sknot* m,sknot* g)
{
  if (m->op!=0) {
    char cc[2];
    string sop;
    cc[1] = 0;cc[0]=m->op;
    sop = string(cc);
    if (m->op=='*') sop = string(" ");

    if (m->op=='/')       
      return string("\\frac{")+Tree2Tex(m->left,m)+
	string("}{")+Tree2Tex(m->right,m)+string("}");
    if (g==0 || m->op=='*')       
      return Tree2Tex(m->left,m)+sop+Tree2Tex(m->right,m);
    if (g->op=='+') {
      if (m->op=='-' && m->left->op==0) {
	if (m->left->Str()==string("0")) 
	  return string("(")+Tree2Tex(m->left,m)+sop
	    +Tree2Tex(m->right,m)+string(")");
      }
      else 
	return Tree2Tex(m->left,m)+sop+Tree2Tex(m->right,m);
    }
    if (m->op=='+' || m->op=='-') 
      return string("(")+Tree2Tex(m->left,m)+sop
	+Tree2Tex(m->right,m)+string(")");
  }
  if (m->Str()==string("0")) return string("");
  
  return m->Str();
}

void String_Tree::CollectLeafs(sknot* leaf,vector<sknot*>& sklist,int full) 
{
  sknot* s = leaf; 

  while (s->op=='*') {
    sklist.push_back(s->right);
    s = s->left;
  }

  sklist.push_back(s);

  //setting Tree2String's
  if (full) {
    for (size_t i=0;i<sklist.size();i++) {
      sklist[i]->SetString(Tree2String(sklist[i],0));
    }
  }

  for (size_t i=0;i<sklist.size();i++) 
    for (size_t j=i+1;j<sklist.size();j++) {
      if (full) {
	if (sklist[i]->Str()>sklist[j]->Str()) {
	  //change
	  s         = sklist[i];
	  sklist[i] = sklist[j];
	  sklist[j] = s;
	}
      }
      else {
	if (sklist[i]->op!=0 && sklist[j]->op==0) {
	  //change
	  s         = sklist[i];
	  sklist[i] = sklist[j];
	  sklist[j] = s;
	}
	else {
	  if (sklist[i]->op==0 && sklist[j]->op==0) {
	    if (sklist[i]->Str()>sklist[j]->Str()) {
	      //change
	      s         = sklist[i];
	      sklist[i] = sklist[j];
	      sklist[j] = s;
	    }
	  }
	}
      }
    }
}


int String_Tree::CountFactorNumber(sknot* leaf1,vector<sknot*>*& list1,
				   sknot* leaf2,vector<sknot*>*& list2,
				   int full)
{
  list1 = new vector<sknot*>;
  list2 = new vector<sknot*>;

  CollectLeafs(leaf1,*list1,full);
  CollectLeafs(leaf2,*list2,full);
  
  int count = 0;
  int i,j;
  i = j = 0;

  do {
    if (full==0) {
      if ((*list1)[i]->op!=0) break;
      if ((*list2)[j]->op!=0) break;
    }
    if ((*list1)[i]->Str()==(*list2)[j]->Str()) {
      //put them to the beginning
      sknot* s1 = (*list1)[i];
      for (int k=i;k>count;k--) (*list1)[k] = (*list1)[k-1];
      (*list1)[count] = s1;
      sknot* s2 = (*list2)[j];
      for (int k=j;k>count;k--) (*list2)[k] = (*list2)[k-1];
      (*list2)[count] = s2;

      count++;i++;j++;
    }
    else {
      if ((*list1)[i]->Str()<(*list2)[j]->Str()) i++;
                                            else j++;
    }
  }
  while (i<(int)(list1->size()) && j<(int)(list2->size()));

  //equal lists will not be clustered at this stage
  if (count==(int)(list1->size()) || count==(int)(list2->size())) count = 0;

  // *AS* *TG*
  // beschleunigung clustern bei superamplitude
  if (full) {
    if (count!=(int)(list1->size())-1 || count!=(int)(list2->size())-1 )
      count =0;
      }


  return count;
}

void String_Tree::Cluster(sknot* m,sknot* g,int full)
{
  if (m==0) return;
  if (m->op==0) return;

  int sw = 0;
  if (g==0) sw = 1;
  else {
    if (g->op!='+' && g->op!='-') sw = 1;
  }

  int hit;
  do {
    hit = 0;
    if ((m->op=='+' || m->op=='-') && sw) {
      //head of PM-chain
      sknot* s1 = m;
      
      int winner = 0;
      sknot* winknot1=0,*winknot2=0;

      vector<sknot*>* winnerlist1 = 0; 
      vector<sknot*>* winnerlist2 = 0; 

      int minus1 = 1;
      int winnerminus = 1;
     
      do {
	sknot* s2  = s1;
	int minus2 = minus1;
	do {
	  if (s2->op=='-') minus2 *= -1;
	  s2 = s2->right;
	  	  
	  sknot* leaf1 = s1->left;
	  sknot* leaf2;	

	  if (s2->op=='+' || s2->op=='-') leaf2 = s2->left;
	  else leaf2 = s2;
	  vector<sknot*>* list1; 
	  vector<sknot*>* list2; 
	  int n = CountFactorNumber(leaf1,list1,leaf2,list2,full);
	  if (n>winner) {
	    winner = n;
	    winknot1 = s1;
	    winknot2 = s2;
	    if (winnerlist1!=0) {
	      delete winnerlist1;
	      delete winnerlist2;
	    }

	    winnerlist1 = list1;
	    winnerlist2 = list2;  
	    if (minus1!=minus2) winnerminus = -1;
	                   else winnerminus = 1;
	  }
	  else {
	    delete list1;
	    delete list2;
	  }
	  if (full && winner) break;
	  
	}
	while (s2->op=='+' || s2->op=='-');
	
	if (s1->op=='-') minus1 *= -1;
	s1 = s1->right;
	
      if (full && winner) break;

      }
      while (s1->op=='+' || s1->op=='-');

      hit = winner;
      if (winner>0) {
	sknot* leaf1,*leaf2;
	leaf1 = winknot1->left;
	if (winknot2->op=='+' || winknot2->op=='-') leaf2 = winknot2->left;
	                                       else leaf2 = winknot2;

	if (winknot1->right!=leaf2) {
	}
	else {
	  leaf2 = winknot1;
	}
	//remove first leaf from list
	if (winknot1->op=='-') {
	  /*sknot *sk = newsk();
	  sk->op    = 0;
	  sk->SetString(string("0"));
	  sk->left  = 0;
	  sk->right = 0;*/
	  winknot1->left = &zero;
	  //test
	  if (leaf2==winknot1) leaf2 = winknot1->right;
	}
	else {
	  winknot1->op    = winknot1->right->op;
	  winknot1->left  = winknot1->right->left;	  
	  winknot1->right = winknot1->right->right;
	}

	sknot* si = leaf2;
	sknot* sj = leaf1;
	sknot* sprev = 0;
	if (winner<(int)(winnerlist1->size())) {
	  for (int i=0;i<winner;i++) {
	    si->right = (*winnerlist1)[i];    
	    sprev = si;
	    si = si->left;
	    sj = sj->left;
	  }
	  sknot* spm  = leaf1;
	  if (winnerminus==-1) spm->op  = '-';
	                  else spm->op  = '+';
	  sprev->left = spm;
	  //left side = winner2
	  if (winnerlist2->size()-winner>1) {
	    spm->left = si;
	    for (size_t i=winner;i<winnerlist2->size()-1;i++) {
	      si->right = (*winnerlist2)[i];   
	      sprev = si;
	      si    = si->left;
	    }
	    sprev->left = (*winnerlist2)[winnerlist2->size()-1];	      
	  }
	  else spm->left = (*winnerlist2)[winnerlist2->size()-1];
	  
	  if (winnerlist1->size()-winner>1) {
	    spm->right = sj;
	    for (int i=winner;i<(int)(winnerlist1->size())-1;i++) {
	      sj->right = (*winnerlist1)[i];   
	      sprev     = sj;
	      sj        = sj->left;
	    }
	    sprev->left = (*winnerlist1)[winnerlist1->size()-1];	      
	  }
	  else spm->right = (*winnerlist1)[winnerlist1->size()-1];	      
	}
	else {
	  abort();
	}
	
	delete winnerlist1;
	delete winnerlist2;
      }    
    }
  }
  while (hit>0);

  Cluster(m->left,m);
  Cluster(m->right,m);
}

void String_Tree::DeleteMinus(sknot* &m)
{
  //  PROFILE_HERE;
  if (m->op=='*'){
    int nminus = CountMinus(m);
    if (nminus%2) {
      sknot* totalminus = newsk();
      totalminus->op    = '-';
      totalminus->left  = &zero;
      totalminus->right = m;
      
      m = totalminus;      
    }
  }
  int hit;
  do {
    hit = 0;
    SingleDeleteMinus(m,hit);
  }
  while (hit>0); 
  LinearOrderPM(m);
}

int String_Tree::CountMinus(sknot* &m) 
{
  if (m==0) return 0;
  if (m->op==0) return 0;
  if (m->op=='+' || m->op=='-') return 0;

  int count = 0;

  if (m->op=='*' && m->left->op=='-') {
    //Endpoint
    if (m->left->left->op==0 &&  m->left->left->Str()==string("0")) {
      //manipulation
      m->left = m->left->right;
      count++;
    }
  }
  if (m->op=='*' && m->right->op=='-') {
    //Endpoint
    if (m->right->left->op==0 &&  m->right->left->Str()==string("0")) {
      //manipulation
      m->right = m->right->right;
      count++;
    }
  }

  return count+CountMinus(m->left)+CountMinus(m->right);
}

void String_Tree::SingleDeleteMinus(sknot* &m,int& hit)
{
  if (m==0) return;
  if (m->op==0) return;

  if (m->op=='+' || m->op=='-') {
    if (m->left->op=='*') {
      int nminus = CountMinus(m->left);
      //manipulation....
      if (nminus%2) {
	if (m->op=='+') {
	  m->op = '-';
	  sknot* shelp = m->left;
	  m->left      = m->right;
	  m->right     = shelp;
	}
	else {
	  sknot* newplus = newsk();
	  newplus->op    = '+';
	  newplus->left  = m->left;
	  newplus->right = m->right;

	  m->left  = &zero;
	  m->right = newplus;
	}
      }

      if (nminus>0) {
	hit = 1;
	return;
      }
    }
    
    if (m->right->op=='*') {
      int nminus = CountMinus(m->right);
      //manipulation....
      if (nminus%2) {
	if (m->op=='+') m->op = '-';
	           else m->op = '+';
      }
      if (nminus>0) {
	hit = 1;
	return;
      }
    }
    
    if (m->left->op=='-') { 
      if (m->left->left->op==0 &&  m->left->left->Str()==string("0")) {
	if (m->op=='+') {
	  m->op = '-';
	  sknot* shelp = m->left->right;
	  m->left      = m->right;
	  m->right     = shelp;
	}
	else {
	  sknot* shelp   = m->right;
	  m->right       = m->left;
	  m->right->op   = '+';
	  m->left        = m->right->left;
	  m->right->left = shelp;
	}
	hit = 1;
	return;
      }
    }

    if (m->right->op=='-') { 
      if (m->right->left->op==0 &&  m->right->left->Str()==string("0")) {
	m->right = m->right->right;
	if (m->op=='-') m->op = '+';
	           else m->op = '-';
	hit = 1;
	return;
      }
    }
  }

  SingleDeleteMinus(m->left,hit);
  if (hit==0) SingleDeleteMinus(m->right,hit);
}

void String_Tree::Delete(sknot* &m,const string& delstr) 
{
  int hit;
  do {
    Single_Delete(m,0,delstr);
    if (m->op==0) {
      if (m->Str()==delstr) m->SetString(string("0"));
    }
    hit = Tree2String(m,0).find(delstr);
  }
  while (hit>=0);
}

void String_Tree::Simplify(sknot*& m)
{
  if (m==0) return;
  if (m->op==0) return;

  if (m->op=='+' || m->op=='-') {
    Complex vleft  = Evaluate(m->left);
    Complex vright = Evaluate(m->right);
    
    if (ATOOLS::IsZero(vleft/(vleft+vright))) {
      //kill left part
      if (m->op=='-') {
	if (m->left->op!=0) m->left  = String2Tree(string("0"));
      }
      else m = m->right;
    }

    if (ATOOLS::IsZero(abs(vright/(vleft+vright))))
      //kill right part
      m = m->left;
  }
  
  Simplify(m->left);
  Simplify(m->right);
}


void String_Tree::Single_Delete(sknot* &m,sknot* g,const string& delstr) 
{
  if (m==0) return;
  if (m->op==0) return;

  int hit = 0;

  if (m->left->op==0) {
    if (m->left->Str()==delstr) {
      hit = 1;
      switch (m->op) {
      case '*': m = m->left;break;
      case '+': m = m->right;break;
      case '-':  
	if (g==0) m->left->SetString(string("0"));
	else {
	  switch (g->op) {
	  case '*': m->left->SetString(string("0"));break;
	  case '+': 
	    if (g->left==m) m->left->SetString(string("0"));
	    else {
	      g->op = '-';
	      m = m->right;
	    }
	    break;
	  case '-': 
	    if (g->left==m) m->left->SetString(string("0"));
	    else {
	      g->op = '+';
	      m = m->right;
	    }
	    break;
	  }
	}
	break;
      }
    }
  }
  if (hit==0) {
    if (m->right->op==0) {
      if (m->right->Str()==delstr) {
	if (m->op=='*') m = m->right;
	else {
	  if (m->op=='e') {
	    m = m->left;
	    m->SetString(string("1"));
	  }
	  else {
	    if (m->op=='s') {
	      m = m->right;
	    }
	    else {
	      if (m->left->Str() == string("0")) m = m->right;
	                                    else m = m->left;
	    }
	  }
	}
      }
    }
  }

  Single_Delete(m->left,m,delstr);
  Single_Delete(m->right,m,delstr);
}

void String_Tree::DetermineDepth(sknot* m,char oldop,int& currdepth)
{
  if (m==0) return;
  
  if (m->op!=0) {
    if ((m->op=='+' || m->op=='-') && oldop=='*') {
      currdepth++;
      oldop = '+';
    }
    else {
      if (m->op=='*' && oldop=='+') {
	currdepth++;
	oldop = '*';
      }
    }
    int savedepth = currdepth;
    DetermineDepth(m->left,oldop,currdepth);
    DetermineDepth(m->right,oldop,savedepth); 
    if (savedepth>currdepth) currdepth = savedepth;
  }
}

void String_Tree::ExpandToDepth(sknot* m,int depth,list<sknot*>& addlist)
{
  Addends(m,addlist);

  int hit;
  do {
    hit = 0;
    for (list<sknot*>::iterator it=addlist.begin();it!=addlist.end();) {
      char oldop = '+';
      int  currdepth = 0;
      DetermineDepth(*it,oldop,currdepth);
      if (currdepth>depth) {
	int dummy=0;
	Single_Expand(*it,dummy);
	Addends(*it,addlist);
	it = addlist.erase(it);
      }
      else ++it;
    }
  }
  while (hit>0);
}



void String_Tree::Expand(sknot* m)
{
  int hit;
  do {
    hit = 0;
    Single_Expand(m,hit);
  }
  while (hit==1);
}

void String_Tree::Single_Expand(sknot* m,int& hit)
{
  if (m==0) return;
  if (hit==1) return;
  if (m->op=='*') {
    int changed = 0;
    if (m->left->op=='+' || m->left->op=='-') {
      sknot* c = m->left;
      m->left  = m->right;
      m->right = c;
      changed = 1;
    }
    if ((m->right->op=='+') || (m->right->op=='-')) {
      if (m->right->op=='-' && m->right->left->op==0) {
	if (m->right->left->Str()==string("0")) {
	  sknot* c1 = m->left;
	  m->left = m->right->left;
	  m->right->left = c1;
	  m->op          = m->right->op;
	  m->right->op = '*';
	  hit  = 1;
	}
      }
      if (hit!=1) {
	sknot* m_new = newsk();
	sknot* c1    = m->left;
	sknot* c2    = Copy(m->left);
	
	m->left        = m_new;
	if (changed) {
	  m->left->left  = m->right->left;
	  m->left->right = c1;
	  
	  m->right->left  = m->right->right;
	  m->right->right = c2;
	}
	else {
	  m->left->right = m->right->left;
	  m->left->left  = c1;
	  m->right->left = c2;
	}
	
	m->op        = m->right->op;
	m->left->op  = '*';
	m->right->op = '*';
	hit = 1;
      }
    }
  }
  if (m->op=='+' || m->op=='-') {
    if (m->right->op=='-') {
      if (m->right->left->op==0) {
	if (m->right->left->Str()==string("0")) {
	  if (m->op=='+') m->op = '-';
	             else m->op = '+';
	  m->right = m->right->right;
	  hit = 1;
	}
      }
    }
  }

  Single_Expand(m->left,hit);
  Single_Expand(m->right,hit);
} 


void String_Tree::Linear(sknot* m)
{
  if (m==0) return;
  if (m->op=='*') {
    for (;;) {
      if ((m->left->op=='*') && (m->right->op=='*')) {
	sknot* c       = m->right->left;
	m->right->left = m->left;
	m->left        = m->right;
	m->right       = c;
      }
      else {
	if ((m->right->op=='*')) {
	  sknot* c = m->right;
	  m->right = m->left;
	  m->left  = c;
	}
	else break;
      }
    }
  }
  Linear(m->left);
  Linear(m->right);
}

void String_Tree::LinearOrderPM(sknot* m)
{
  LinearPM(m);
  OrderPM(m,0);
  Linear(m);
}

void String_Tree::DetermineLeafAndSign(sknot* m,vector<sknot*>&pmleafs,vector<int>&pmsigns,int& currentsign)
{
  if (m==0) return;
  if (m->op==0) return;

  if (m->op=='+' || m->op=='-') {
    pmleafs.push_back(m->left);
    pmsigns.push_back(currentsign);   

    if (m->op=='-') currentsign*= -1;
 
    if (m->right->op!='+' && m->right->op!='-') {
      pmleafs.push_back(m->right);
      pmsigns.push_back(currentsign);  
    }
    else DetermineLeafAndSign(m->right,pmleafs,pmsigns,currentsign);
  }
}

void String_Tree::SetLeafAndSign(sknot* m,vector<sknot*>&pmleafs,vector<int>&pmsigns,int& count)
{
  if (m==0) return;
  if (m->op==0) return;

  if (m->op=='+' || m->op=='-') {
    m->left = pmleafs[count];
    if (pmsigns[count]==1) m->op = '+';
                      else m->op = '-';
    count++;
   
    if (count==(int)pmleafs.size()-1) {
      m->right = pmleafs[count];
      count++;
    }
    else SetLeafAndSign(m->right,pmleafs,pmsigns,count);
  }
}

void String_Tree::DeleteEquals(vector<sknot*>&pmleafs,vector<int>&pmsigns)
{
  int signchange = 0;

  for (int i=0;i<(int)pmleafs.size()-1;i++) {
    if (pmsigns[i+1]==-1 && pmsigns[i]==+1) {
      signchange = i+1;
      break;
    }
  }
  
  if (signchange==0) return;

  int hit;
 
  do {
    hit = 0;
    for (int i=0;i<signchange;i++) {
      for (int j=signchange;j<(int)pmleafs.size();j++) {
	if (pmleafs[i]->op==0 && pmleafs[j]->op==0) {
	  if (pmleafs[i]->Str()==pmleafs[j]->Str()) {
	    //kick them
	    for (int k=i;k<j-1;k++) {
	      pmleafs[k] = pmleafs[k+1];
	      pmsigns[k] = pmsigns[k+1];
	    }
	    for (size_t k=j-1;k<pmleafs.size()-2;k++) {
	      pmleafs[k] = pmleafs[k+2];
	      pmsigns[k] = pmsigns[k+2];
	    }    
	    pmleafs.pop_back();
	    pmleafs.pop_back();
	    pmsigns.pop_back();
	    pmsigns.pop_back();
	    signchange--;
	    if (signchange==0) return;
	    hit = 1;
	    break;
	  }
	}
      }
      if (hit) break;
    }
  }
  while (hit>0);
      

}
void String_Tree::OrderPM(sknot* m,sknot* g)
{
  if (m==0) return;
  if (m->op==0) return;

  int sw = 0;
  if (g==0) sw = 1;
  else {
    if (g->op!='+' && g->op!='-') sw = 1;
  }

  if ((m->op=='+' || m->op=='-') && sw) {
    if (m->right->op=='+' || m->right->op=='-') {
      //Head of a plus minus line......
      vector<sknot*> pmleafs;
      vector<int>    pmsigns;
      int currentsign = 1;
      DetermineLeafAndSign(m,pmleafs,pmsigns,currentsign);
      //order leaf and sign
      for (size_t i=0;i<pmleafs.size();i++)
	for (size_t j=i+1;j<pmleafs.size();j++) {
	  if (pmsigns[i]==-1 && pmsigns[j]==1) {
	    //exchange
	    sknot* shelp = pmleafs[i];
	    pmleafs[i]   = pmleafs[j];
	    pmleafs[j]   = shelp;
	    pmsigns[i]   = 1;
	    pmsigns[j]   = -1;
	  }
	}
      
      DeleteEquals(pmleafs,pmsigns);

      if (pmleafs.empty()) {
	m->op = 0;
	m->left = 0;
	m->right = 0;
	m->SetString(string("Z[0]"));
      }
      else {
	sknot* startknot = m;
	if (pmsigns[0]==-1) {
	  //complete minus list
	  m->left   = &zero;
	  m->op     = '-'; 
	  startknot = m->right;
	}

	for (int i=0;i<(int)pmleafs.size()-1;i++) {
	  if (pmsigns[i+1]==-1 && pmsigns[i]==+1) pmsigns[i] = -1;
	                                     else pmsigns[i] = +1;
	}
	 
	if (pmleafs.size()==1) {
	  startknot->op    = pmleafs[0]->op;
	  startknot->left  = pmleafs[0]->left;
	  startknot->right = pmleafs[0]->right;
	  startknot->SetString(pmleafs[0]->Str());
	}
	else {
	  int count = 0;
	  SetLeafAndSign(startknot,pmleafs,pmsigns,count);
	}
      }
    }
  }

  OrderPM(m->left,m);
  OrderPM(m->right,m);
}

void String_Tree::LinearPM(sknot* m)
{
  if (m==0) return;
  
  if ((m->op=='+') || (m->op=='-')) {
    for (;;) {
      if (m->left->op=='+') {
	sknot* c        = m->right;
	m->right        = m->left;
	m->left         = m->right->right;
	m->right->right = c;
	
	if (m->op=='-') {
	  m->right->op = '-';
	  m->op        = '+';
	}
      }
      else {
	if (m->left->op=='-') {
	  sknot* c        = m->right;
	  m->right        = m->left;
	  m->left         = m->right->left;
	  m->right->left  = m->right->right;
	  m->right->right = c;
	  
	  if (m->op=='+') m->op = '-';
	             else m->right->op = '+';
	}
	else break;
      }
    }
  }

  if (m->op!=0) {
    LinearPM(m->left);
    LinearPM(m->right);
  }
}

void String_Tree::Sort(sknot* m)
{
  if (m==0) return;
  if (m->op=='*') {
    int hit;
    sknot* akt = m;
    do {
      hit = 0;
      do {
	if (m->left->op=='*') {
	  if (strcmp((m->right->Str()).c_str(),(m->left->right->Str()).c_str())<0) {
	    sknot* help = m->right;
	    m->right = m->left->right;
	    m->left->right = help;
	    hit  = 1;
	  }
	}
	else {
	  if (m->left->op==0) {
	    if (strcmp((m->right->Str()).c_str(),(m->left->Str()).c_str())<0) {
	      sknot* help     = m->right;
	      m->right = m->left;
	      m->left  = help;
	      hit      = 1;	  
	    }
	  }
	}
	m = m->left;
      }
      while (m->op=='*');
      m = akt;
    }
    while (hit);
  }

  Sort(m->left);
  Sort(m->right);
}

sknot* String_Tree::Copy(sknot* m,int fixed) 
{
  if (m==0) return m;
  
  if (m->op!=0) {
    sknot* s = newsk(); 
    s->op  = m->op;
    s->left  = Copy(m->left);
    s->right = Copy(m->right);
    return s;
  }

  if (fixed) {
    // search for Leaf
    for (list<sknot*>::iterator it=leaflist.begin();it!=leaflist.end();++it) {
      if ((*it)->op==0) {
	if ((*it)->Str()==m->Str()) return (*it); 
      }
    }
  }

  if(m->Str()==string("0")){
    if (fixed) leaflist.push_back(&zero);
    return &zero;
  }

  //totally new sknot
  sknot* s = newsk(); 
  s->op    = 0;
  s->SetString(m->Str());
  s->left  = 0;
  s->right = 0;
  //add to leaflist
  if (fixed) leaflist.push_back(s);
  
  return s;
} 
