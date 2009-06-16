#include "AMEGIC++/Main/Topology.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Topology::Topology() 
{
  ntop = 0;
  top = 0;
}

Topology::Topology(int Nmax) 
{
  //Topology
  ntop = Nmax+1;
  Build_All(ntop);
}

Topology::~Topology() 
{
  int i,j;
  if (top) {
    for (i=0;i<ntop;i++) {
      for (j=0;j<top[i].number;j++) 
	delete[] top[i].p[j];
      delete[] top[i].p;
    }
    delete[] top;
  }
}

void Topology::Build_All(int N)
{
  top = new Single_Topology[N];
  top[0].number  = 1;
  top[0].depth   = 1;
  top[0].p    = new Point*[1];
  top[0].p[0] = new Point[1]; 
  top[0].p[0][0].left  = 0;
  top[0].p[0][0].right = 0;

  for(int i=1;i<N;i++) Build_Single(i+1,top);
}

void Topology::Print(Point* p)
{
  if (!msg_LevelIsDebugging()) return;
  if (p==0) {msg_Out()<<"End."<<endl;return;}

  msg_Out()<<"Left - ";
  Print(p->left);
  msg_Out()<<"Right - ";  
  Print(p->right);
  if (p->middle!=0) { 
    msg_Out()<<"Middle - ";  
    Print(p->middle);
  }
}

Point* Topology::Copy(Point* po,Point* pc,int& ll)
{
  pc[ll] = *po;

  if (po->left==0) {
    pc[ll].left   = 0;
    pc[ll].right  = 0;
    pc[ll].middle = 0;
    ll++;
    return &pc[ll-1];
  }
  int lsave = ll;
  ll++;
  pc[lsave].left   = Copy(po->left,pc,ll);
  pc[lsave].right  = Copy(po->right,pc,ll); 
  if (po->middle!=0) pc[lsave].middle = Copy(po->middle,pc,ll); 
  return &pc[lsave];
}

void Topology::Build_Single(int nlegs,Single_Topology* t) 
{

  int newnumber = 0;
  // t[0] = 1jet
  // t[1] = 2jet
  // t[2] = 3jet etc.

  int i,j,k;
  
  for (i=0;i<nlegs-1;i++) {
     newnumber += t[i].number * t[nlegs-i-2].number;
  }
  t[nlegs-1].number = newnumber;
  t[nlegs-1].depth  = nlegs;
  t[nlegs-1].p = new Point*[newnumber];
  for (i=0;i<newnumber;i++) t[nlegs-1].p[i] = new Point[2*nlegs-1]; 

  newnumber = 0;
  int ll; 

  for (i=0;i<nlegs-1;i++) {
    // i -> left
    for (j=0;j<t[i].number;j++) {
       //nlegs-i-2 -> right
      for (k=0;k<t[nlegs-i-2].number;k++) {
	ll = 1;
	t[nlegs-1].p[newnumber][0].left  = Copy(t[i].p[j],t[nlegs-1].p[newnumber],ll);
	t[nlegs-1].p[newnumber][0].right = Copy(t[nlegs-i-2].p[k],t[nlegs-1].p[newnumber],ll); 
	t[nlegs-1].depth = ll;
	newnumber++;
      }
    }	   
  }
} 









