#include "Single_Broker.H"
#include "XS_Base.H"
#include "Single_XS.H"
#include "Process_Base.H"
#include "Run_Parameter.H"
#include "Standard_Selector.H"
#include "Message.H"
#include "Combined_Selector.H"

#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Polarisation.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace BEAM;
using namespace ISR;
using namespace std;

Single_Broker::Single_Broker(EXTRAXS::XS_Base * _xsec) {
  xsec     = _xsec;
  nin      = _xsec->Nin();
  nout     = _xsec->Nout();
  name     = _xsec->Name();

  {
    Flavour * _fl = _xsec->Flavs();
    flin   = new Flavour[nin];
    flout  = new Flavour[nout];
    for (short int i=0;i<nin;i++) {
      flin[i]  = _fl[i];
    }
    for (short int i=nin;i<nin+nout;i++) {
      flout[i-nin]  = _fl[i];
    }
    Polarisation pol;
    int                   is_massless_pol  = pol.Massless_Vectors(nin,flin);
    if (!is_massless_pol) is_massless_pol  = pol.Massless_Vectors(nout,flout);
    int                   nmassive_pols    = pol.Massive_Vectors(nin,flin);
    nmassive_pols                         += pol.Massive_Vectors(nout,flout);
    nvec = nin+nout+is_massless_pol+nmassive_pols;
    delete [] flin;
    delete [] flout;
  }

  fl  = new Flavour[nin+nout];
  for (short int i=0;i<nin+nout;i++) {
    fl[i]  = _xsec->Flavs()[i];
  }

  selected = this;

  kfactorscheme = rpa.me.KFactorScheme();
  scalescheme   = rpa.me.ScaleScheme();

  CreateSelector();

  n        = 0;
  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = lastlumi    = 0.;

  taumin   = 0.;
  taumax   = 1.;
}

Single_Broker::Single_Broker(int _nin,int _nout, APHYTOOLS::Flavour * _fl,
			     ISR::ISR_Handler * _isr,BEAM::Beam_Handler * _beam,
			     APHYTOOLS::Selector_Data * _sel, int _kfactorscheme, int _scalescheme) {
  nin  = _nin;
  nout = _nout;
  kfactorscheme = _kfactorscheme;
  scalescheme   = _scalescheme;

  fl  = new Flavour[nin+nout];
  for (short int i=0;i<nin+nout;i++) {
    fl[i]  = _fl[i];
  }

  isr  = _isr;
  beam = _beam;

  xsec     = new EXTRAXS::Single_XS(nin, nout, fl);

  Initialize(_sel);

  taumin   = 0.;
  taumax   = 1.;
}
  
Single_Broker::~Single_Broker() {
  if (xsec)     { delete xsec;                   }

  if (fl)       { delete [] fl;    fl       = 0; }
  if (flin)     { delete [] flin;  flin     = 0; }
  if (flout)    { delete [] flout; flout    = 0; }
  if (pl)       { delete [] pl;    pl       = 0; }
  if (plin)     { delete [] plin;  plin     = 0; }
  if (plout)    { delete [] plout; plout    = 0; }
  if (b)        { delete [] b;     b        = 0; }
  if (moms)     { delete [] moms;  moms     = 0; }
 
  if (sel)      { delete sel;      sel      = 0; }
  if (cuts)     { delete cuts;     cuts     = 0; }
  if (ps)       { delete ps;       ps       = 0; }
}

void Single_Broker::Initialize(APHYTOOLS::Selector_Data * _seldata) {
  cuts = 0;
  InitCuts();
  ps   = new Phase_Space_Handler(this,isr,beam);

  if (_seldata) sel = new Combined_Selector(nin,nout,fl,_seldata);
  else {
    msg.Error()<<"Potential Error in Single_Broker "<<name<<endl
	       <<" No selection cuts specified. Init No_Selector !"<<endl;
    sel = new No_Selector();
  }

  moms = 0;

  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = 0.;
  lastlumi = 1.;

  selected = this;
  atoms    = 1;
  analyse  = 0;
}

bool Single_Broker::SetUpIntegrator(ISR_Handler * _isr,Beam_Handler * _beam) {
  isr  = _isr;
  beam = _beam;

//   isr->SetSprimeMin(taumin*s);
//   isr->SetSprimeMax(taumax*s);

  ps  = new Phase_Space_Handler(this,isr,beam);
  CreateSelector();
  ps->CreateIntegrators();
  return 1;
}

bool Single_Broker::CalculateTotalXSec() {
  msg.Out()<<"In Single_Broker::CalculateTotalXSec() for "<<name<<endl; 
  totalxs = ps->Integrate()/AORGTOOLS::rpa.Picobarn();
  if (!(AMATOOLS::IsZero((n*totalxs-totalsum)/(n*totalxs+totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<endl;
    msg.Error()<<"  "<<name<<" : "<<totalxs<<" vs. "<<totalsum/n<<endl;
  }
  SetTotalXS(0);
  if (totalxs>0.) return 1;
  return 0;
}

void Single_Broker::AddPoint(const double value)
{
  n++;
  totalsum    += value;
  totalsumsqr += value*value;
  if (value>max) max = value;
} 

void Single_Broker::SetTotalXS(int hint)  { 
  totalxs  = totalsum/n; 
  totalerr = sqrt( (totalsumsqr/n - 
		    (AMATOOLS::sqr(totalsum)-totalsumsqr)/n/(n-1) )  / n); 
  AORGTOOLS::msg.Events()<<"      xs for "<<name<<" : "
			 <<totalxs*AORGTOOLS::rpa.Picobarn()<<" pb"
			 <<" +/- "<<totalerr/totalxs*100.<<"%"<<endl;
}

bool Single_Broker::OneEvent() {
  return ps->OneEvent();  
}

double Single_Broker::Differential(AMATOOLS::vec4d * p) { 
  return DSigma(p, 0); 
}

double Single_Broker::Differential2() { 
  return DSigma2(); 
}

double Single_Broker::DSigma(AMATOOLS::vec4d * p, bool lookup) { 
  s = (p[0]+p[1]).abs2();
  t = (p[0]-p[2]).abs2();
  u = (p[0]-p[3]).abs2();
  lastdxs = xsec->operator()(s,t,u);
  lastlumi = isr->Weight(fl);
  return last = lastdxs * lastlumi;
}

double Single_Broker::DSigma2() { 
  if (fl[0]==fl[1]) return 0.;
  double tmp = lastdxs * isr->Weight2(fl); 
  last += tmp;
  return tmp;
}

bool Single_Broker::SameEvent() { 
  AORGTOOLS::msg.Error()<<"Virtual Method : Single_Broker::SameEvent()."<<std::endl; 
  return 0.; 
}

double Single_Broker::WeightedEvent() { 
  AORGTOOLS::msg.Error()<<"Virtual Method : Single_Broker::WeightedEvent()."<<std::endl; 
  return 0.; 
}

double Single_Broker::KFactor(double _scale) { 
  return xsec->KFactor(_scale); 
}


