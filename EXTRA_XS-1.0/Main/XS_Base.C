#include "XS_Base.H"

#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Message.H"
#include <stdio.h>

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


XS_Base::XS_Base(int _nin,int _nout,Flavour * _fl,
		 PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
		 ATOOLS::Selector_Data * _seldata,
		 int _scalescheme,int _kfactorscheme,double _scalefactor) :
  m_nin(_nin), m_nout(_nout),m_name(std::string("")),
  m_scalescheme(_scalescheme), m_kfactorscheme(_kfactorscheme), m_scalefactor(_scalefactor),
  m_n(0), m_last(0.), m_lastlumi(0.), m_lastdxs(0.),
  m_totalxs(0.),m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.),
  p_selected(NULL), p_beam(_beam), p_isr(_isr), p_sel(NULL), p_ps(NULL), 
  p_fl(NULL), p_colours(NULL), p_moms(NULL)
{
  Init(_fl);
  ResetSelector(_seldata);
  p_ps   = new PHASIC::Phase_Space_Handler(this,p_isr,p_beam);
  p_moms = new ATOOLS::Vec4D[m_nin+m_nout];
}

XS_Base::XS_Base(int _nin,int _nout,Flavour * _fl) :
  m_nin(_nin), m_nout(_nout),m_name(std::string("")),
  m_scalescheme(0), m_kfactorscheme(0), m_scalefactor(1.),
  m_n(0), m_last(0.), m_lastlumi(0.), m_lastdxs(0.),
  m_totalxs(0.),m_totalsum (0.), m_totalsumsqr(0.), m_totalerr(0.),
  p_selected(NULL), p_sel(NULL), p_beam(NULL), p_isr(NULL), p_ps(NULL), 
  p_fl(NULL), p_colours(NULL), p_moms(NULL)
{
  Init(_fl);
  p_sel = new No_Selector();
  p_moms = new ATOOLS::Vec4D[m_nin+m_nout];
}

XS_Base::~XS_Base() {
  if (p_fl)       { delete [] p_fl;   p_fl   = 0; }
  if (p_moms)     { delete [] p_moms; p_moms = 0; }
  if (p_sel)     { delete    p_sel;  p_sel  = 0; }
  if (p_colours) { 
    for (int i=0;i<m_nin+m_nout;i++) delete p_colours[i];
    delete [] p_colours; p_colours = 0;
  }
}

void XS_Base::Init(Flavour * _fl)
{
  p_fl = new Flavour[m_nin+m_nout];
  if (_fl) {
    for (short int i=0;i<m_nin+m_nout;i++) p_fl[i]  = _fl[i];
    GenerateName();
  }
  p_colours = new int*[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) { 
    p_colours[i]    = new int[2]; 
    p_colours[i][0] = p_colours[i][1] = 0; 
  }  

  double massin = 0., massout =0.;
  for (int i=0;i<m_nin;i++)      massin  += p_fl[i].Mass();
  for (int i=m_nin;i<m_nout;i++) massout += p_fl[i].Mass();
  if (massin>massout) m_thres = ATOOLS::sqr(massin);
                 else m_thres = ATOOLS::sqr(massout);
}

void XS_Base::GenerateName() {
  char help[20];
  sprintf(help,"%i",m_nin);
  m_name       = string(help);
  m_name      += string("_");
  sprintf(help,"%i",m_nout);
  m_name      += string(help);
  m_name      += string("_");

  for (int i=0;i<m_nin;i++) {
    m_name    += string(p_fl[i].Name());
    if ((p_fl[i].Kfcode()==kf::e)   ||
	(p_fl[i].Kfcode()==kf::mu)  ||
	(p_fl[i].Kfcode()==kf::tau) ||
	(p_fl[i].Kfcode()==kf::Hmin)) {
      //kill last
      m_name.erase(m_name.length()-1,1);
      if (p_fl[i].IsAnti()) m_name += string("+");
                       else m_name += string("-");      
    }
    else {
      if (p_fl[i].IsAnti()) m_name += string("b"); 
    }
    m_name += string("_");
  }
  m_name.erase(m_name.length()-1,1);

  m_name      += string(" -> ");
  for (int i=m_nin;i<m_nin+m_nout;i++) {
    m_name    += string(p_fl[i].Name());
    if ((p_fl[i].Kfcode()==kf::e)   ||
	(p_fl[i].Kfcode()==kf::mu)  ||
	(p_fl[i].Kfcode()==kf::tau) ||
	(p_fl[i].Kfcode()==kf::Hmin)) {
      //kill last
      m_name.erase(m_name.length()-1,1);
      if (p_fl[i].IsAnti()) m_name += string("+");
                       else m_name += string("-");      
    }
    else {
      if (p_fl[i].IsAnti()) m_name += string("b"); 
    }
    m_name += string("_");
  }
  m_name.erase(m_name.length()-1,1);
}

double XS_Base::Differential(ATOOLS::Vec4D * p) {
  for (int i=0;i<m_nin+m_nout;i++) p_moms[i] = p[i];
  m_s = (p[0]+p[1]).Abs2();
  m_t = (p[0]-p[2]).Abs2();
  m_u = (p[0]-p[3]).Abs2();
  return Differential(m_s,m_t,m_u);
}

bool XS_Base::SetColours(ATOOLS::Vec4D * p) {
  for (int i=0;i<m_nin+m_nout;i++) p_moms[i] = p[i];
  m_s = (p[0]+p[1]).Abs2();
  m_t = (p[0]-p[2]).Abs2();
  m_u = (p[0]-p[3]).Abs2();
  // scale can be overwritten in SetColors(s,t,u)
  m_scale = sqr(p[2][1])+sqr(p[2][2]);  
  return SetColours(m_s,m_t,m_u);
}

double XS_Base::Scale(ATOOLS::Vec4D * p) {
  for (int i=0;i<m_nin+m_nout;i++) p_moms[i] = p[i];
  if (m_nin==1) return p[0].Abs2();
  m_s = (p[0]+p[1]).Abs2();
  m_t = (p[0]-p[2]).Abs2();
  m_u = (p[0]-p[3]).Abs2();

  double pt2;
  switch (m_scalescheme) {
  case 1  :
    if (m_nin+m_nout==4) {
      pt2 = ATOOLS::sqr(p[2][1])+ATOOLS::sqr(p[2][2]);
      //pt2 = 2.*m_s*m_t*m_u/(m_s*m_s+m_t*m_t+m_u*m_u);
    }
    return m_scale = pt2;
  case 2  :
    return m_scale;
  default :
    return m_scale = m_s;
  }
}

double XS_Base::KFactor(double _scale) {
  switch (m_kfactorscheme) {
  case 1  :
    return pow(as->AlphaS(_scale * m_scalefactor)/
	       as->AlphaS(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())),m_nin+m_nout-2);
  default :
    return 1.;
  }
}



void XS_Base::SwapInOrder() {
  Flavour help = p_fl[0];
  p_fl[0] = p_fl[1];
  p_fl[1] = help;
  Vec4D mom = p_moms[0];
  p_moms[0] = p_moms[1];
  p_moms[1] = mom;
  m_swaped = 1;
}

void XS_Base::RestoreInOrder() {
  if (m_swaped) {
    Flavour help = p_fl[0];
    p_fl[0] = p_fl[1];
    p_fl[1] = help;
    m_swaped = 0;
  }
}

void XS_Base::ResetSelector(ATOOLS::Selector_Data *_seldata)
{
  if (p_sel!=NULL) delete p_sel;
  if (_seldata) p_sel = new Combined_Selector(m_nin,m_nout,p_fl,_seldata);
  else {
    msg.Error()<<"Potential Error in Single_Process "<<m_name<<endl
	       <<"   No selection cuts specified. Init No_Selector !"<<endl;
    p_sel = new No_Selector();
  }
}


