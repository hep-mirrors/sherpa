#include "Zfunc.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

#define new_calc


void Zfunc::ReplaceProp(vector<Pair>* pairlist)
{
  for(int i=0;i<narg;i++){
    for (int k=0;k<pairlist->size();k++) {
      if ((*pairlist)[k].pold==arg[i]) {
	arg[i] = (*pairlist)[k].pnew;
	break;
      }
    }
  }
  for (int i=0;i<pn;i++) {
    for (int k=0;k<pairlist->size();k++) {
      if ((*pairlist)[k].pold==psnew[i].numb) {
	psnew[i].numb = (*pairlist)[k].pnew;
	break;
      }
    }
  }
}

void Zfunc::ClearCalcList()
{
  if (!value.empty()) {
    value.clear();
    for (int i=0;i<calclist.size();i++) {
      delete[] calclist[i];
    }
    calclist.clear();
  }
}

void Zfunc::Print() 
{
  msg.Out()<<"Z(["<<type<<"],";
  msg.Out()<<"[";
  for (int i=0;i<narg-1;i++) msg.Out()<<arg[i]<<";";
  
  if (narg>0) msg.Out()<<arg[narg-1];
  msg.Out()<<"][";
  msg.Out().precision(2);
  for (int i=0;i<ncoupl-1;i++) {
    if ( !AMATOOLS::IsZero(real(coupl[i])) &&
	    AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<real(coupl[i])<<";";
    if (  AMATOOLS::IsZero(real(coupl[i])) &&
	  !AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<imag(coupl[i])<<" I;";
    if ( !AMATOOLS::IsZero(real(coupl[i])) &&
	 !AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<real(coupl[i])<<"+"<<imag(coupl[i])<<" I;";
    if (  AMATOOLS::IsZero(real(coupl[i])) &&
	  AMATOOLS::IsZero(imag(coupl[i])) )
      msg.Out()<<"0;";
  }
  if ( !AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
       AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
      msg.Out()<<real(coupl[ncoupl-1])<<"])";
  if (  AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
	!AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
    msg.Out()<<imag(coupl[ncoupl-1])<<" I])";
  if ( !AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
       !AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
    msg.Out()<<real(coupl[ncoupl-1])<<"+"<<imag(coupl[ncoupl-1])<<" I])";
  if (  AMATOOLS::IsZero(real(coupl[ncoupl-1])) &&
	AMATOOLS::IsZero(imag(coupl[ncoupl-1])) )
	msg.Out()<<"0])";
  msg.Out()<<endl;
  msg.Out().precision(6);
}

Zfunc_Group::Zfunc_Group(Zfunc& z1,Zfunc& z2,int si,list<Pfunc*>* pl)
{
  int n99=0;                    //Marker 99
  int i1=0;
  for(int i=0;i<z1.narg;i++){
    if (z1.arg[i]==si) i1++;
    if (z1.arg[i]==99) n99++;   //polarisation dummys not needed in Zfunc_Group
  }
  int i2=0;
  for(int i=0;i<z2.narg;i++){
    if (z2.arg[i]==si) i2++;
    if (z2.arg[i]==99) n99++;
  }
  if((i1==0&&i2>0) || (i2==0&&i1>0)){
    msg.Error()<<"Error in Zfunc_Group(Z*Z-Constructor): sum index"<<endl;
    abort();
  }

  type   = zl::Unknown;
  narg   = z1.narg + z2.narg - i1 - i2 - n99;
  ncoupl = z1.ncoupl+z2.ncoupl;
  pn=0;

  int ii=0;
  if(narg>0){
    arg = new int[narg];
    while (ii<narg) {
      for (int i=0;i<z1.narg;i++) 
	if( z1.arg[i]!=si && z1.arg[i]!=99 ){
	  arg[ii] = z1.arg[i];
	  ii++;
	}
      for (int i=0;i<z2.narg;i++) 
	if( z2.arg[i]!=si && z2.arg[i]!=99 ){
	  arg[ii] = z2.arg[i];
	  ii++;
	}
    }
  }

  if(z1.GetOp()=='*')pn+=z1.pn;
  if(z1.GetOp()==0){
    for(int i=0;i<z1.pn;i++){
      for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc*  p = *pit;
	if (p->arg[0]==z1.psnew[i].numb && p->on==0){
	  pn++;
	  break;
	}
      } 
    }
  }
  int k=pn;

  if(z2.GetOp()=='*')pn+=z2.pn;
  if(z2.GetOp()==0){
    for(int i=0;i<z2.pn;i++){
      for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc*  p = *pit;
	if (p->arg[0]==z2.psnew[i].numb && p->on==0){
	  pn++;
	  break;
	}
      } 
    }
  }

  pn++;

  if (pn>0) {
    psnew = new Argument[pn];
    if((z1.GetOp()=='*')&&(z1.pn>0)){
      for (int i=0;i<k;i++) psnew[i] = z1.psnew[i];
      k=z1.pn;
      delete[] z1.psnew;
      z1.pn=0;
    }
    if(z1.GetOp()==0){
      int j=0;
      for(int i=0;i<z1.pn;i++){
	for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
	  Pfunc*  p = *pit;
	  if (p->arg[0]==z1.psnew[i].numb && p->on==0){
	    psnew[j] = z1.psnew[i];
	    j++;
	    break;
	  }
	} 
      }
    }
    if((z2.GetOp()=='*')&&(z2.pn>0)){
      for (int i=0;i<z2.pn;i++) psnew[i+k] = z2.psnew[i];
      delete[] z2.psnew;
      z2.pn=0;      
    }
    if(z2.GetOp()==0){
      int j=k;
      for(int i=0;i<z2.pn;i++){
	for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
	  Pfunc*  p = *pit;
	  if (p->arg[0]==z2.psnew[i].numb && p->on==0){
	    psnew[j] = z2.psnew[i];
	    j++;
	    break;
	  }
	} 
      }
    }
    psnew[pn-1].numb=si;
  }
	
  if(ncoupl>0){
    coupl = new Complex[ncoupl];
    for (int i=0;i<z1.ncoupl;i++) coupl[i] = z1.coupl[i];
    for (int i=z1.ncoupl;i<z1.ncoupl+z2.ncoupl;i++) coupl[i] = z2.coupl[i-z1.ncoupl];
  }

  Equal      = this;
  sign=1;
  sumindex = si;
  op='*';
  zlist.push_back(&z1);
  zsign.push_back(1);
  zlist.push_back(&z2);
  zsign.push_back(1);
 //z1.Print();
  //z2.Print();
  //Print();
}


void Zfunc_Group::ReplaceProp(vector<Pair>* pairlist)
{
  for (int k=0;k<pairlist->size();k++) {
    if ((*pairlist)[k].pold==sumindex) {
      sumindex = (*pairlist)[k].pnew;
      break;
    }
  }
  Zfunc::ReplaceProp(pairlist);
  for(int i=0;i<zlist.size();i++) zlist[i]->ReplaceProp(pairlist);
}

void Zfunc_Group::ClearCalcList()
{
  Zfunc::ClearCalcList();
  for(int i=0;i<zlist.size();i++) zlist[i]->ClearCalcList();
}

void Zfunc_Group::Print() 
{
  msg.Out()<<"SZ(["<<type<<"],";
  msg.Out()<<"[";
  for (int i=0;i<narg-1;i++) msg.Out()<<arg[i]<<";";
  
  if (narg>0) msg.Out()<<arg[narg-1];
  msg.Out()<<"])";
  msg.Out()<<endl;

  if(op=='+'){  
    for (int i=0;i<zlist.size();i++) {
      if (zsign[i]==-1) {
	cout<<"   - "<<zlist[i]->psnew[0].numb<<" * ";zlist[i]->Print();
      }
      else {
	if (zlist[i]->psnew!=NULL) {
	  cout<<"   + "<<zlist[i]->psnew[0].numb<<" * ";
	  zlist[i]->Print();
	}
	else {
	  cout<<" ??? "<<" * ";zlist[i]->Print();
	}
      }
    }
  }
  if(op=='*'){
    for (int i=0;i<zlist.size();i++) {
      if(i>0)cout<<"  *";else cout<<" ->";
      zlist[i]->Print();
    }
    cout<<"Sum over "<<sumindex<<endl;
  }  
}
