#include "Color_Generator.H"

using namespace AMEGIC;
using namespace std;

void Color_Generator::CFConvert(Point* p)
{
  if ((p->left==0) && (p->right==0)) return;

  if (p->Color->type!=cf::None) {
    Color_Function* CFh  = new Color_Function;
    *CFh = *(p->Color); 
    
    //pointer on the head of the CFh list
    Color_Function* clhead = CFh;
    
    if (p->Color->Next==0) CFh->Next = 0;
    
    //convert the Colorlist CFh
    while (CFh) {
      
      if (CFh->type!=cf::None) {
	for (short int i=0;i<3;i++) {
	  if ((CFh->type==cf::D || CFh->type==cf::G) && i==2) break;
	  switch (CFh->partarg[i]) {
	  case 0: CFh->partarg[i] = p->number;break;
	  case 1: CFh->partarg[i] = p->left->number;break;
	  case 2: CFh->partarg[i] = p->right->number;break;
	  case 3: CFh->partarg[i] = p->middle->number;break;
	  }
	}
      }
      CFh = CFh->Next;
    }   
    
    //connect Colorlist (clhead) and CFlist 
    if (CFlist) {
      Color_Function* help = CFlist;
      while (help->Next) help = help->Next;
      help->Next = clhead;  
      if (clhead->Next!=0) {
	Color_Function* cforig = clhead->Next;
	Color_Function* cfcopy = help->Next;
	while (cforig) {
	  cfcopy->Next = new Color_Function(*cforig);
	  cfcopy = cfcopy->Next;
	  cforig = cforig->Next;
	}
      }
    }
    else {
      CFlist = clhead;
      if (clhead->Next!=0) {
	Color_Function* cforig = clhead->Next;
	Color_Function* cfcopy = CFlist;
	while (cforig) {
	  cfcopy->Next = new Color_Function(*cforig);
	  cfcopy = cfcopy->Next;
	  cforig = cforig->Next;
	}
      }
    }
  }
  CFConvert(p->right);
  if (p->middle) CFConvert(p->middle);
  CFConvert(p->left);
}

string Color_Generator::CF2String(Color_Function* cflist)
{
  Color_Function* CFh = cflist;

  string stringchain("");
  
  while (CFh) {
    if (stringchain.length()>1) stringchain += string("*");
    stringchain += CFh->String();
    CFh = CFh->Next;
  }

  return stringchain;
}

void Color_Generator::FillString(int N, Color_Function* cflist,int& prop) 
{
  char ca = 'A';
  char ci = 'i';

  Color_Function* CFh = cflist;
  
  while (CFh) {
    //if (CFh->type==cf::None) break;    
    for (short int i=0;i<3;i++) {
      if ((CFh->type==cf::D || CFh->type==cf::G) && i==2) break;
      if ((CFh->strarg[i]>=48 && CFh->strarg[i]<=52)) {
	char chelp;
	switch (CFh->type) {
	case cf::F: {
	  if (CFh->partarg[i]<99) chelp = ca+CFh->partarg[i];
  	                     else chelp = ca+(prop++)+N;
	}
	  break;
	case cf::T: 
	  if (i==0) {
	  if (CFh->partarg[i]<99) chelp = ca+CFh->partarg[i];
  	                     else chelp = ca+(prop++)+N;
	  }
	  else {
	    if (CFh->partarg[i]<99) chelp = ci+CFh->partarg[i];
	                       else chelp = ci+(prop++)+N;
	  }
	  break; 
	case cf::D: 
	  if (CFh->partarg[i]<99) chelp = ci+CFh->partarg[i];
	                     else chelp = ci+(prop++)+N;
	  
	  break;
	case cf::G: 
	  if (CFh->partarg[i]<99) chelp = ca+CFh->partarg[i];
	                     else chelp = ca+(prop++)+N;
	  
	  break;
	default :
	  break;
	}
	
	Color_Function* CFh2 = CFh;
	while (CFh2) {
	  //if (CFh2->type==cf::None) break;
	  for (short int j=0;j<3;j++) {
	    if ((CFh2->type==cf::D || CFh2->type==cf::G) && j==2) break;
	    if (CFh2->partarg[j]==CFh->partarg[i]) CFh2->strarg[j] = chelp;
	  }
	  CFh2 = CFh2->Next;
	}
      }
    }
    CFh = CFh->Next;
  }
}


void Color_Generator::CFBuildString(int N)
{
  int prop = 0;
  
  FillString(N,CFlist,prop);

  //Conjugated Color list

  Color_Function* CFh = CFlist;

  while (CFh) {

    Color_Function* CFh2 = new Color_Function;
    *CFh2 = *CFh;
    CFh2->Next = 0;

    //clear strarg's
    for (short int i=0;i<3;i++) CFh2->strarg[i] = '0';

    //exchange indizes -> conjugated matrix
    if (CFh2->type==cf::T) {
      int help = CFh2->partarg[1];
      CFh2->partarg[1] = CFh2->partarg[2];
      CFh2->partarg[2] = help;
      
      //not necessary anymore
      /*
      char help2 = CFh2->strarg[1];
      CFh2->strarg[1] = CFh2->strarg[2];
      CFh2->strarg[2] = help2;  
      */    
    }

    if (CCFlist) {
      Color_Function* help;
      help = CCFlist;
      while (help->Next) help = help->Next; 
      help->Next = CFh2;  
    }
    else CCFlist = CFh2;

    CFh = CFh->Next;
  }
  
  prop+=2; //Savety
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
  short int i;
  while (c) {
    if (c->type == cf::D || c->type == cf::G) {
      replace = -1;
      if (c->partarg[0]>99) {
	replace = c->partarg[0];
	with    = c->partarg[1];
      }
      if (c->partarg[1]>99) {
	replace = c->partarg[1];
	with    = c->partarg[0];
      }
      if (replace!=-1) {
	c2 = CFlist;
	while (c2) {
	  if (c2!=c) {
	    for (i=0;i<3;i++) {
	      if ((c2->type==cf::D || c2->type==cf::G) && i==2) break;
	      if (c2->partarg[i]==replace) c2->partarg[i]=with;
	    }
	  }
	  c2 = c2->Next;	    
	}
      }
    }
    c = c->Next;
  }
  // kill all TP with propagator
    
  Color_Function* last;
  last = CFlist;
  c = CFlist;
  while (c) {
    if ((c->type==cf::D || c->type==cf::G) && ((c->partarg[0]>99) || (c->partarg[1]>99))) {
      if (c==CFlist) {
	CFlist = c->Next;
	c = c->Next;
	delete last;
	last = CFlist;
      }
      else {
	last->Next = c->Next;
	c2 = c;
	c = c->Next;
	delete c2;
      }
    }
    else {
      last = c;
      c = c->Next;
    }
  }
}
















