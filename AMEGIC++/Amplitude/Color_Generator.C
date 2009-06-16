#include "AMEGIC++/Amplitude/Color_Generator.H"

using namespace AMEGIC;
using namespace MODEL;
using namespace std;

void Color_Generator::CFConvert(int N, int& dummy, Point* p)
{
  if ((p->left==0) && (p->right==0)) return;
    
  if (p->Color->Type()!=cf::None) {
    Color_Function* CFh  = new Color_Function;
    *CFh = *(p->Color); 
    
    //pointer on the head of the CFh list
    Color_Function* clhead = CFh;
    
    //if (p->Color->Next()==0) CFh->Next = 0;
    
    //convert the Colorlist CFh
    while (CFh) {
      if (CFh->Type()!=cf::None) {
	int partarg[3] = {-1,-1,-1}; 
	for (short int i=0;i<3;i++) {
	  if ((CFh->Type()==cf::D || CFh->Type()==cf::G) && i==2) break;
	  switch (CFh->ParticleArg(i)) {
	  case 0: partarg[i]  = p->number;break;
	  case 1: partarg[i]  = p->left->number;break;
	  case 2: partarg[i]  = p->right->number;break;
	  case 3: partarg[i]  = p->middle->number;break;
	  default: partarg[i] = N+1+dummy/2;dummy++;
	  }
	}
	CFh->SetParticleArg(partarg[0],partarg[1],partarg[2]);
      }
      CFh = CFh->Next();
    }   
    
    //connect Colorlist (clhead) and CFlist 
    if (CFlist) {
      CFlist->Append(clhead);
    }
    else {
      CFlist = clhead;
    }
  }
  CFConvert(N,dummy,p->right);
  if (p->middle) CFConvert(N,dummy,p->middle);
  CFConvert(N,dummy,p->left);
}

string Color_Generator::CF2String(Color_Function* cflist)
{
  Color_Function* CFh = cflist;

  string stringchain("");
  
  while (CFh) {
    if (stringchain.length()>1) stringchain += string("*");
    stringchain += CFh->String();
    CFh = CFh->Next();
  }
  return stringchain;
}

void Color_Generator::FillString(int N, Color_Function* cflist,int& prop) 
{
  char ca = 'A';
  char ci = 'i';

  Color_Function* CFh = cflist;
  
  while (CFh) {
    //if (CFh->Type()==cf::None) break;    
    for (short int i=0;i<3;i++) {
      if ((CFh->Type()==cf::D || CFh->Type()==cf::G) && i==2) break;
      if ((CFh->StringArg(i)>=48 && CFh->StringArg(i)<=52)) {
	char chelp(0);
	int arg = CFh->ParticleArg(i);
	switch (CFh->Type()) {
	case cf::F: {
	  if (arg<99) {
	    if (arg>N) chelp = ca+(prop++)+N;
	      else chelp = ca+arg;
	  }
	  else chelp = ca+(prop++)+N;
	}
	  break;
	case cf::T: 
	  if (i==0) {
	    if (arg<99) chelp = ca+arg;
	                   else chelp = ca+(prop++)+N;
	  }
	  else {
	    if (arg<99) chelp = ci+arg;
	                       else chelp = ci+(prop++)+N;
	  }
	  break; 
	case cf::D: 
	  if (arg<99) chelp = ci+arg;
	                     else chelp = ci+(prop++)+N;
	  
	  break;
	case cf::G: 
	  if (arg<99) chelp = ca+arg;
	                     else chelp = ca+(prop++)+N;
	  
	  break;

	default: break;
	}
	
	Color_Function* CFh2 = CFh;
	while (CFh2) {
	  CFh2->Replace(CFh->ParticleArg(i),chelp);
	  CFh2 = CFh2->Next();
	}
      }
    }
    CFh = CFh->Next();
  }
}


void Color_Generator::CFBuildString(int N)
{
  int prop = 0;
  
  FillString(N,CFlist,prop);

  //Conjugated Color list

  Color_Function * CFh  = CFlist;
  Color_Function * CFh2 = 0;

  if (CFh) CFh2 = new Color_Function(*CFh);  // copy complete list
  CCFlist = CFh2;

  while (CFh2) {
    //clear strarg's
    CFh2->SetStringArg('0','0','0');

    //exchange indizes -> conjugated matrix
    if (CFh2->Type()==cf::T) {
      CFh2->Conjugate();
    }
    CFh2 = CFh2->Next();
  }
  prop+=N-3;
  FillString(N,CCFlist,prop);
}

void Color_Generator::CFKill() 
{
  // replace all type=10 TP(P,i,j) with delta(i,j)
  Color_Function* c;
  Color_Function* c2;
  c = CFlist;
  int replace,with;
  replace = with  = 0;
  while (c) {
    if (c->Type() == cf::D || c->Type() == cf::G) {
      replace = -1;
      if (c->ParticleArg(0)>99) {
	replace = c->ParticleArg(0);
	with    = c->ParticleArg(1);
      }
      if (c->ParticleArg(1)>99) {
	replace = c->ParticleArg(1);
	with    = c->ParticleArg(0);
      }
      if (replace!=-1) {
	c2 = CFlist;
	while (c2) {
	  if (c2!=c) c2->Replace(replace,with);
	  c2 = c2->Next();	    
	}
      }
    }
    c = c->Next();
  }
  // kill all TP with propagator
    
  Color_Function* last;
  last = CFlist;
  c = CFlist;
  while (c) {
    if ((c->Type()==cf::D || c->Type()==cf::G) && ((c->ParticleArg(0)>99) || (c->ParticleArg(1)>99))) { 
      if (c==CFlist) {
	c = CFlist = last->Erase(0);
	last = CFlist;
      }
      else {
	c=c->Erase(last);
      }
    }
    else {
      last = c;
      c = c->Next();
    }
  }
}
