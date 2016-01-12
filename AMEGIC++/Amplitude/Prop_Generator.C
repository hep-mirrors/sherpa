#include "AMEGIC++/Amplitude/Prop_Generator.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include "AMEGIC++/Amplitude/Zfunc.H"
#include "AMEGIC++/Main/Point.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include <cassert>

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

void Prop_Generator::Convert(Point* p)
{
  if ((p->left==0) && (p->right==0)) return;  

  if (p->number>99) {
    Pfunc* Ph = new Pfunc;
    Ph->on = 1;
    Ph->haspol = p->m; 
    Ph->zerowidth = p->zwf; 

    Ph->fl = p->fl;
    if (p->middle) Ph->argnum = 4;
              else Ph->argnum = 3; 

    Ph->arg = new int[Ph->argnum];
    Ph->arg[0] = p->number;
    Ph->arg[1] = p->left->number;
    Ph->arg[2] = p->right->number;

    if (p->middle) Ph->arg[3] = p->middle->number;

    plist.push_back(Ph);

    if (p->nextra>0) {
      for (short int i=0;i<p->nextra;i++) {
	Pfunc* Ph2  = new Pfunc;
	Ph2->haspol = 0;    
	Ph2->zerowidth = p->zwf; 
	Ph2->on     = 1;	
	Ph2->fl = p->extrafl[i];
	Ph2->argnum = Ph->argnum;
	Ph2->arg = new int[Ph2->argnum];
	Ph2->arg[0] = p->number*(p->extrafl[i]).Kfcode();
	Ph2->arg[1] = p->left->number;
	Ph2->arg[2] = p->right->number;
	if (p->middle) Ph2->arg[3] = p->middle->number;

	plist.push_back(Ph2);
      }
    }
  }
  Convert(p->right);
  Convert(p->left);
  if (p->middle) Convert(p->middle);
}

void Prop_Generator::Fill()
{
  int sw1;
  for (;;) {
    sw1 = 1;
    for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
      Pfunc* p = *pit;
      for (int i=1;i<p->argnum;i++) {
	if (p->arg[i]>99) {
	  sw1 = 0;
	  int* harg;
	  harg = new int[p->argnum];
	  for (int j=0;j<p->argnum;j++) harg[j] = p->arg[j];
	  Pfunc* ph = NULL;
	  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
	    ph = *pit;
	    if (ph->arg[0]==harg[i]) break;
	  }
	  assert(ph !=NULL);
	  p->argnum += ph->argnum-2;
	  delete[] p->arg;
	  p->arg = new int[p->argnum];
	  for (int j=0;j<i;j++) p->arg[j] = harg[j];
	  for (int j=i;j<i+ph->argnum-1;j++) p->arg[j] = ph->arg[j-i+1];
	  for (int j=i+ph->argnum-1;j<p->argnum;j++) p->arg[j] = harg[j-ph->argnum+2];
	  delete[] harg;
	  break;
	}
      }
    }
    if (sw1) break;
  }
}

void Prop_Generator::Kill(Zfunc_List& zlist) 
{
  //neglect props without sum

  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc* p = *pit;
    int sw1 = 1;
    //if (p->arg[0]>199) p->on = 0;
    for (Zfunc_Iterator zit=zlist.begin();zit!=zlist.end();++zit) {
      Zfunc* z = (*zit);
      for (int i=0;i<z->m_narg;i++) {
	if (z->p_arguments[i]==p->arg[0]) {
	  sw1 = 0;
	  break;
	}
      }
      if (sw1==0) break;
    }

    if (sw1 && !p->fl.IsScalar()) p->on = 0;
  }
}

void Prop_Generator::Get(Pfunc_List& _plist) {
  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) 
    _plist.push_back(*pit);
}
