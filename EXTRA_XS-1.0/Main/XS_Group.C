#include "XS_Group.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "MathTools.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace EXTRAXS;
using namespace PHASIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace std;

XS_Group::XS_Group(int _nin,int _nout,string _name) 
{
  Init(_nin,_nout,0);
  name   = _name; 
  taumin = 0.;
  taumax = 1.;
}

XS_Group::~XS_Group()
{
  for(int i=xs.size();i>0;i--) {
    if (xs[i-1]) delete xs[i-1];
  }
}

void XS_Group::Add(XS_Base * _xs) 
{
  if (xs.size()==0) {
    nin  = _xs->Nin();
    nout = _xs->Nout();
  }
  else {
    if ( (nin != _xs->Nin()) || (nout != _xs->Nout())) {
      AORGTOOLS::msg.Error()<<"Error : Cannot add xs "<<_xs->Name()
			    <<" to group "<<name<<" ! "<<endl
			    <<"   Inconsistent number of external legs."<<endl 
			    <<"  Before : ("<<nin<<" -> "<<nout<<" )"<<endl
			    <<"  Now    : ("<<_xs->Nin()<<" -> "<<_xs->Nout()<<" )"<<endl;
      return;
    }
  }  
  AORGTOOLS::msg.Tracking()<<"Add xs "<<_xs->Name()<<" to group "<<name<<" ! "<<endl; 
  xs.push_back(_xs);
};


void XS_Group::AddPoint(double value) 
{
  n++;
  totalsum    += value;
  totalsumsqr += value*value;
  if (value>max) max = value;
  for (int i=0;i<xs.size();i++) {
    if (dabs(last)>0.) xs[i]->AddPoint(value*xs[i]->Last()/last);
    else               xs[i]->AddPoint(0.);
  }  
}

void XS_Group::SetTotalXS()  { 
  totalxs  = totalsum/n; 
  totalerr = sqrt( (totalsumsqr/n - 
		    (AMATOOLS::sqr(totalsum)-totalsumsqr)/n/(n-1) )  / n); 
  if (sel) sel->Output();
  if (isr->On()) msg.Events()<<"  ISR range : "
			     <<isr->SprimeMin()<<" ... "<<isr->SprimeMax()<<endl;
  max = 0.;
  for (int i=0;i<xs.size();i++) {
    xs[i]->SetTotalXS();
    max += xs[i]->Max();
  }
  msg.Events()<<"------------------------------------------------------------------------------"<<endl;
  msg.Events()<<"Total XS for "<<name<<"("<<xs.size()<<") : "<<totalxs*rpa.Picobarn()<<" pb";
  msg.Events()<<" +/- "<<totalerr/totalxs*100.<<"%"<<endl;
}


void XS_Group::SelectOne()
{
  if (totalxs==0) selected = xs[int(Ran.get()*xs.size())];
  else {
    double disc = max * Ran.get();
    for (int i=0;i<xs.size();i++) {
      disc -= xs[i]->Max();
      if (disc<0.) {
	selected = xs[i];
	selected->SelectOne();
	return;
      }
    }
    if (disc>0.) { 
      msg.Error()<<"Error in XS_Group::SelectOne() : ";
      msg.Error()<<"Total xsec, max = "<<totalxs<<", "<<max<<endl;
      return;
    }
  }
}

void XS_Group::DeSelect() {
  selected = 0;
  for (int i=0;i<xs.size();i++) xs[i]->DeSelect();
}

XS_Base * XS_Group::Selected() { 
  if (selected==0)    return 0; 
  if (selected==this) return this;
  return selected->Selected(); 
}    


double XS_Group::Differential(vec4d * p) {
  s = (p[0]+p[1]).abs2();
  t = (p[0]-p[2]).abs2();
  u = (p[0]-p[3]).abs2();
  return Differential(s,t,u);
};

double XS_Group::Differential(double s,double t,double u) { 
  last = 0.;
  for (int i=0;i<xs.size();i++) last += xs[i]->DSigma(s,t,u);
  if ((!(last<=0)) && (!(last>0))) {
    msg.Error()<<"---- XS_Group::Differential -------------------"<<endl;
  }
  return last;
}

double XS_Group::Differential2() { 
  double tmp = 0.;
  for (int i=0;i<xs.size();i++) tmp += xs[i]->DSigma2();

  if ((!(tmp<=0)) && (!(tmp>0))) {
    msg.Error()<<"---- XS_Group::Differential -------------------"<<endl;
  }
  last += tmp;
  return tmp;
}

double XS_Group::DSigma(double s,double t,double u)
{
  last = 0;
  for (int i=0;i<xs.size();i++) last += xs[i]->DSigma(s,t,u);
  return last;
}

double XS_Group::DSigma2()
{
  double tmp = 0.;
  for (int i=0;i<xs.size();i++) tmp += xs[i]->DSigma2();
  last += tmp;
  return tmp;
}

void XS_Group::SetISR(ISR::ISR_Handler * _isr) { 
  msg.Debugging()<<"XS_Group::SetISR("<<_isr->Type()<<") for "<<name<<"  : "<<_isr<<endl;
  isr = _isr; 

//   isr->SetSprimeMin(taumin*s);
//   isr->SetSprimeMax(taumax*s);
  for (int i=0;i<xs.size();i++) xs[i]->SetISR(isr);
}
