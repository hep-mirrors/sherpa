#include "XS_Base.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace BEAM;
using namespace ISR;
using namespace std;

XS_Base::~XS_Base() {
  if (fl)      { delete [] fl; fl = 0; }
  if (ps)      { delete ps; ps = 0; }
  if (sel)     { delete sel; sel = 0; }
  if (colours) { 
    for (int i=0;i<nin+nout;i++) delete colours[i];
    delete [] colours; colours = 0;
  }
};

void XS_Base::Init(int _nin,int _nout,Flavour * _fl) {
  msg.Debugging()<<"In XS_Base::Init("<<_nin<<","<<_nout<<")"<<endl;
  nin     = _nin; nout = _nout;
  fl      = new Flavour[nin+nout];
  if (_fl) {
    for (short int i=0;i<nin+nout;i++) {
      fl[i]  = _fl[i];
      msg.Debugging()<<fl[i]<<" ";
    }
  }
  msg.Debugging()<<endl;

  colours = new int*[nin+nout];
  for (int i=0;i<nin+nout;i++) colours[i] = new int[2];

  selected = this;

  n        = 0;
  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = lastlumi    = 0.;

  taumin = 0.; taumax = 1.;
}

void XS_Base::ISRInfo(int i,int & type,double & mass,double & width) {
  msg.Debugging()<<"In XS_Base::ISRInfo("<<i<<")"<<endl;
  if ((i<0) || (i>=isr_types.size())) {
    msg.Error()<<"ERROR : XS_Base::ISRInfo("<<i<<") out of bounds !"<<endl;
    abort();
  }
  type  = isr_types[i];
  mass  = isr_masses[i];
  width = isr_widths[i];
  msg.Debugging()<<"Out XS_Base::ISRInfo("<<i<<") :"<<type<<","<<mass<<endl;
};


bool XS_Base::SetUpIntegrator(ISR_Handler * _isr,Beam_Handler * _beam) {
  isr  = _isr;
  beam = _beam;

//   isr->SetSprimeMin(taumin*s);
//   isr->SetSprimeMax(taumax*s);

  ps  = new Phase_Space_Handler(this,isr,beam);
  CreateSelector();
  ps->CreateIntegrators();
  return 1;
}

bool XS_Base::CalculateTotalXSec() {
  msg.Out()<<"In XS_Base::CalculateTotalXSec() for "<<name<<endl; 
  totalxs = ps->Integrate()/AORGTOOLS::rpa.Picobarn();
  if (!(AMATOOLS::IsZero((n*totalxs-totalsum)/(n*totalxs+totalsum)))) {
    msg.Error()<<"Result of PS-Integrator and internal summation to not coincide!"<<endl;
    msg.Error()<<"  "<<name<<" : "<<totalxs<<" vs. "<<totalsum/n<<endl;
  }
  SetTotalXS();
  if (totalxs>0.) return 1;
  return 0;
}

void XS_Base::AddPoint(double value)
{
  n++;
  totalsum    += value;
  totalsumsqr += value*value;
  if (value>max) max = value;
} 

void XS_Base::SetTotalXS()  { 
  totalxs  = totalsum/n; 
  totalerr = sqrt( (totalsumsqr/n - 
		    (AMATOOLS::sqr(totalsum)-totalsumsqr)/n/(n-1) )  / n); 
  AORGTOOLS::msg.Events()<<"      xs for "<<name<<" : "
			 <<totalxs*AORGTOOLS::rpa.Picobarn()<<" pb"
			 <<" +/- "<<totalerr/totalxs*100.<<"%"<<endl;
}

bool XS_Base::OneEvent() {
  return ps->OneEvent();  
}


double XS_Base::Differential(vec4d * p) {
  s = (p[0]+p[1]).abs2();
  t = (p[0]-p[2]).abs2();
  u = (p[0]-p[3]).abs2();
  return Differential(s,t,u);
};

bool XS_Base::SetColours(AMATOOLS::vec4d * p) {
  s = (p[0]+p[1]).abs2();
  t = (p[0]-p[2]).abs2();
  u = (p[0]-p[3]).abs2();
  return SetColours(s,t,u);
}


double XS_Base::Differential(double s,double t,double u) { 
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::Differential()."<<std::endl; 
  return 0.; 
}

double XS_Base::Differential2() { 
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::Differential2()."<<std::endl; 
  return 0.; 
}

double XS_Base::DSigma(double s,double t,double u) { 
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::DSigma()."<<std::endl; 
  return 0.; 
}

double XS_Base::DSigma2() { 
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::DSigma2()."<<std::endl; 
  return 0.; 
}

double XS_Base::operator()(double,double,double) {
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::operator()."<<std::endl; 
  return 0.;
}

bool XS_Base::SetColours(double,double,double) {
  AORGTOOLS::msg.Error()<<"Virtual Method : XS_Base::SetColours()."<<std::endl; 
  return 0.; 
}





double XS_Base::Scale(vec4d * p)     { return scale = (p[0]+p[1]).abs2(); }

double XS_Base::KFactor(double _scale) { return 1.; }

void XS_Base::CreateSelector() { sel = new No_Selector(); }

void XS_Base::SetBeam(BEAM::Beam_Handler * _beam) { 
  msg.Debugging()<<"Error : Virtual method XS_Base::SetISR()"<<endl;
}

void XS_Base::SetISR(ISR::ISR_Handler * _isr) { 
  msg.Debugging()<<"Error : Virtual method XS_Base::SetISR()"<<endl;
}

XS_Base * XS_Base::Selected()  { 
  msg.Debugging()<<"Error : Virtual method XS_Base::Selected()"<<endl;
}

