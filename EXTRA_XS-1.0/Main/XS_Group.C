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
using namespace ATOOLS;
using namespace std;

XS_Group::XS_Group(int _nin,int _nout,Flavour * _fl,
		   PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
		   ATOOLS::Selector_Data * _seldata,
		   int _scalescheme,int _kfactorscheme,double _scalefactor) :
  XS_Base(_nin,_nout,_fl,_isr,_beam,_seldata,_scalescheme,_kfactorscheme,_scalefactor),
  p_xsselector(NULL), m_atoms(0), m_fsrchannels(false)
{
  p_selected = NULL;
}

XS_Group::XS_Group(int _nin,int _nout,Flavour * _fl) :
  XS_Base(_nin,_nout,_fl), p_xsselector(NULL), m_atoms(0)
{
  p_selected = NULL;
}

XS_Group::XS_Group(int _nin,int _nout,std::string _name) :
  XS_Base(_nin,_nout,NULL), p_xsselector(NULL), m_atoms(0)
{
  m_name     = _name;
  p_selected = NULL;
}

XS_Group::XS_Group() : p_xsselector(NULL), m_atoms(0)
{
  p_selected = NULL;
}

XS_Group::~XS_Group()
{
  for(int i=m_xsecs.size();i>0;i--) {
    if (m_xsecs[i-1]) delete m_xsecs[i-1];
  }
  if (p_xsselector) { delete p_xsselector; p_xsselector = NULL; }
}

void XS_Group::Add(XS_Base * _xsec) 
{
  if (m_xsecs.size()==0) {
    m_nin  = _xsec->Nin();
    m_nout = _xsec->Nout();
    
    p_fl  = new Flavour[m_nin+m_nout];
    for (short int i=0;i<m_nin+m_nout;i++) {
      p_fl[i]  = _xsec->Flavs()[i];
    }
  }
  else {
    if ( (m_nin!=_xsec->Nin()) || (m_nout!=_xsec->Nout())) {
      ATOOLS::msg.Error()<<"Error : Cannot add Process "<<_xsec->Name()
			    <<" to group "<<m_name<<" ! "<<endl
			    <<"   Inconsistent number of external legs."<<endl 
			    <<"  Before : ("<<m_nin<<" -> "<<m_nout<<" )"<<endl
			    <<"  Now    : ("<<_xsec->Nin()<<" -> "<<_xsec->Nout()<<" )"<<endl;
      return;
    }
  }  
  ATOOLS::msg.Debugging()<<"Add xs "<<_xsec->Name()<<" to group "<<m_name<<" ! "<<endl; 
  m_xsecs.push_back(_xsec);
}

bool XS_Group::Delete(XS_Base *_xsec) 
{
  for (std::vector<XS_Base*>::iterator xsit=m_xsecs.begin();xsit!=m_xsecs.end();++xsit) {
    if (*xsit==_xsec) {
      delete *xsit;
      m_xsecs.erase(xsit);
      return true;
    }
  }
  return false;
}

void XS_Group::Clear() 
{
  while (m_xsecs.size()>0) {
    delete *m_xsecs.begin();
    m_xsecs.erase(m_xsecs.begin());
  }
}

void XS_Group::SelectOne()
{
  DeSelect();
  if (m_totalxs==0) p_selected = m_xsecs[int(ran.Get()*m_xsecs.size())];
  else {
    double disc;
    if (m_atoms) {
      // select according to total xsecs.
      disc = m_totalxs * ran.Get();
      for (int i=0;i<m_xsecs.size();i++) {
	disc -= m_xsecs[i]->Total();
	if (disc<0.) {
	  p_selected = m_xsecs[i];
	  p_selected->SelectOne();
	  return;
	}
      }
      if (disc>0.) { 
	msg.Error()<<"Error in Process_Group::SelectOne() : "
		   <<"Total xsec, max = "<<m_totalxs<<", "<<m_max<<endl;
	return;
      }
    }
    else {
      disc = m_max * ran.Get();
      for (int i=0;i<m_xsecs.size();i++) {
	disc -= m_xsecs[i]->Max();
	if (disc<0.) {
	  p_selected = m_xsecs[i];
	  p_selected->SelectOne();
	  return;
	}
      }
      if (disc>0.) { 
	msg.Error()<<"Error in Process_Group::SelectOne() : "
		   <<"Total xsec, max = "<<m_totalxs<<", "<<m_max<<endl;
	return;
      }
    }
  }
}



void XS_Group::DeSelect() {
  p_selected = 0;
  for (int i=0;i<m_xsecs.size();i++) m_xsecs[i]->DeSelect();
}



XS_Base * XS_Group::Selected() { 
  if (p_selected==this) return this;
  return p_selected->Selected(); 
}    

void XS_Group::SetISR(PDF::ISR_Handler * _isr) {
  p_isr = _isr;
  for (int i=0;i<m_xsecs.size();i++) m_xsecs[i]->SetISR(_isr);
}

bool XS_Group::CalculateTotalXSec()
{
  m_n=0;
  m_last=m_lastlumi=m_lastdxs=0.0;
  m_totalxs=m_totalsum=m_totalsumsqr=m_totalerr=0.0;
  if (p_isr) {
    if (m_nin==2) {
      if ( (p_fl[0].Mass() != p_isr->Flav(0).Mass()) ||
	   (p_fl[1].Mass() != p_isr->Flav(1).Mass()) ) p_isr->SetPartonMasses(p_fl);
    }
    for (int i=0;i<m_xsecs.size();i++) m_xsecs[i]->SetISR(p_isr);
  }

  if (!m_fsrchannels) {
    CreateFSRChannels();
    p_ps->CreateIntegrators();
  }

  m_totalxs = p_ps->Integrate()/ATOOLS::rpa.Picobarn(); 
  if (!(ATOOLS::IsZero((m_n*m_totalxs-m_totalsum)/(m_n*m_totalxs+m_totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<endl
	       <<"  "<<m_name<<" : "<<m_totalxs<<" vs. "<<m_totalsum/m_n<<endl;
  }
  SetTotalXS();
  if (m_totalxs>0.) return 1;
  return 0;
}


void XS_Group::SetTotalXS()  { 
  m_totalxs  = m_totalsum/m_n; 
  m_totalerr = sqrt( (m_n*m_totalsumsqr-ATOOLS::sqr(m_totalsum))/(m_n-1))/m_n;
  if (p_sel) p_sel->Output();

  m_max = 0.;
  for (int i=0;i<m_xsecs.size();i++) {
    m_xsecs[i]->SetTotalXS();
    m_max += m_xsecs[i]->Max();
  }
  msg.Events()<<"-----------------------------------------------------------------------"<<endl
	      <<"Total XS for "<<m_name<<"("<<m_xsecs.size()<<") : "<<m_totalxs*rpa.Picobarn()<<" pb"
	      <<" +/- "<<m_totalerr/m_totalxs*100.<<"%,"<<endl
	      <<"      max = "<<m_max<<endl;
}


bool XS_Group::OneEvent() {
  if (m_atoms) {
    SelectOne();
    return p_selected->OneEvent();
  }
  return p_ps->OneEvent();
}



void XS_Group::AddPoint(const double value) 
{
  m_n++;
  m_totalsum    += value;
  m_totalsumsqr += value*value;
  if (value>m_max) m_max = value;

  for (int i=0;i<m_xsecs.size();i++) {
    if (dabs(m_last)>0.) {
      m_xsecs[i]->AddPoint(value*m_xsecs[i]->Last()/m_last);
    }
    else {
      m_xsecs[i]->AddPoint(0.);
    }  
  }
}

double XS_Group::Differential(double s,double t,double u)
{
  m_last = 0;
  for (int i=0;i<m_xsecs.size();i++) m_last += m_xsecs[i]->Differential(s,t,u);
  if ((!(m_last<=0)) && (!(m_last>0))) {
    msg.Error()<<"---- X_Group::Differential -------------------"<<endl;
  }
  return m_last;
}



double XS_Group::Differential2()
{
  if (p_isr) {
    if (p_isr->On()==0) return 0.;
    double tmp = 0.;
    for (int i=0;i<m_xsecs.size();i++) tmp += m_xsecs[i]->Differential2();

    if ((!(tmp<=0)) && (!(tmp>0))) {
      msg.Error()<<"---- X_Group::Differential -------------------"<<endl;
    }
    m_last += tmp;
    return tmp;
  }
  return 0.;
}

