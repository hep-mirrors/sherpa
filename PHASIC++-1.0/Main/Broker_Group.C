#include "Broker_Group.H"
#include "Single_Broker.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Combined_Selector.H"
#include "Selector.H"
#include "Flavour.H"

using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace PHASIC;

Broker_Group::Broker_Group(EXTRAXS::XS_Base * _xsec) 
{
  xsecs = _xsec;
  name  = _xsec->Name(); 
  nin   = _xsec->Nin();
  nout  = _xsec->Nout();

  {
    flin   = new Flavour[nin];
    flout  = new Flavour[nout];
    for (short int i=0;i<nin;i++) {
      flin[i]  = (_xsec->Flavs())[i];
    }
    for (short int i=nin;i<nin+nout;i++) {
      flout[i-nin]  = (_xsec->Flavs())[i];
    }
    Polarisation pol;
    int                   is_massless_pol  = pol.Massless_Vectors(nin,flin);
    if (!is_massless_pol) is_massless_pol  = pol.Massless_Vectors(nout,flout);
    int                   nmassive_pols    = pol.Massive_Vectors(nin,flin);
    nmassive_pols                         += pol.Massive_Vectors(nout,flout);
    nvec = nin+nout+is_massless_pol+nmassive_pols;
  }

  fl  = new Flavour[nin+nout];
  for (short int i=0;i<nin+nout;i++) {
    fl[i]  = (_xsec->Flavs())[i]; 
  }

  selected = this;

  kfactorscheme = rpa.me.KFactorScheme();
  scalescheme   = rpa.me.ScaleScheme();

  CreateSelector();

  n        = 0;
  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = lastlumi    = 0.;

  taumin = 0.; taumax = 1.;

  atoms  = 0;
}

Broker_Group::~Broker_Group()
{
  if (flin)  delete [] flin;
  if (flout) delete [] flout;
  if (fl)    delete [] flin;
  if (ps)    delete ps;
  if (sel)   delete sel;

  for(int i=brokers.size();i>0;i--) {
    if (brokers[i-1]) delete brokers[i-1];
  }
}

void Broker_Group::Initialize(APHYTOOLS::Selector_Data * _seldata) {
  cuts = 0;
  InitCuts();
  ps = new Phase_Space_Handler(this,isr,beam);

  if (_seldata) sel = new APHYTOOLS::Combined_Selector(nin,nout,fl,_seldata);
  else {
    msg.Error()<<"Potential Error in Broker_Group "<<name<<std::endl
	       <<" No selection cuts specified. Init No_Selector !"<<std::endl;
    sel = new APHYTOOLS::No_Selector();
  }

  moms = 0;

  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = 0.;
  lastlumi = 1.;

  selected = this;
  analyse  = 0;
}

void Broker_Group::Add(EXTRAXS::XS_Base * _xsec) 
{
  if (brokers.size()==0) {
    nin  = _xsec->Nin();
    nout = _xsec->Nout();
  }
  else {
    if ( (nin != _xsec->Nin()) || (nout != _xsec->Nout())) {
      AORGTOOLS::msg.Error()<<"Error : Cannot add Broker for "<<_xsec->Name()
			    <<" to group "<<name<<" ! "<<std::endl
			    <<" Inconsistent number of external legs."<<std::endl 
			    <<" Before : ("<<nin<<" -> "<<nout<<" )"<<std::endl
			    <<" Now    : ("<<_xsec->Nin()<<" -> "<<_xsec->Nout()
			    <<" )"<<std::endl;
      return;
     }
  }  
  AORGTOOLS::msg.Tracking()<<"Add Broker for "<<_xsec->Name()<<" to group "<<name<<" ! "<<std::endl; 
  Single_Broker * newbroker = new Single_Broker(_xsec);
  brokers.push_back(newbroker);
};


void Broker_Group::AddPoint(double value) 
{
  n++;
  totalsum    += value;
  totalsumsqr += value*value;
  if (value>max) max = value;
  for (int i=0;i<brokers.size();i++) {
    if (dabs(last)>0.) brokers[i]->AddPoint(value*brokers[i]->Last()/last);
    else               brokers[i]->AddPoint(0.);
  }  
}

void Broker_Group::SetTotalXS(int hint)  { 
  msg.Events()<<std::endl;
  msg.Events()<<"--------------------------------------------------------------------------"<<std::endl;
  totalxs  = totalsum/n; 
  totalerr = sqrt( (totalsumsqr/n - 
		    (AMATOOLS::sqr(totalsum)-totalsumsqr)/n/(n-1) ) / n); 
  if (sel) sel->Output();
  if (isr->On()) msg.Events()<<"  ISR range : "
			     <<isr->SprimeMin()<<" ... "<<isr->SprimeMax()<<std::endl;
  max = 0.;
  for (int i=0;i<brokers.size();i++) {
    brokers[i]->SetTotalXS(0);
    max += brokers[i]->Max();
  }
  msg.Events()<<"=========================================================================="<<std::endl;
  msg.Events()<<"Total XS for "<<name<<"("<<brokers.size()<<") : "<<totalxs*rpa.Picobarn()<<" pb";
  msg.Events()<<" +/- "<<totalerr/totalxs*100.<<"%"<<std::endl;
  msg.Events()<<"=========================================================================="<<std::endl;
  msg.Events()<<"      max = "<<max<<std::endl;

}

void Broker_Group::SelectOne()
{
  if (totalxs==0) selected = brokers[int(ran.Get()*brokers.size())];
  else {
    double disc = max * ran.Get();
    for (int i=0;i<brokers.size();i++) {
      disc -= brokers[i]->Max();
      if (disc<0.) {
	selected = brokers[i];
	selected->SelectOne();
	return;
      }
    }
    if (disc>0.) { 
      msg.Error()<<"Error in Broker_Group::SelectOne() : ";
      msg.Error()<<"Total xsec, max = "<<totalxs<<", "<<max<<std::endl;
      return;
    }
  }
}

void Broker_Group::DeSelect() {
  selected = 0;
  for (int i=0;i<brokers.size();i++) brokers[i]->DeSelect();
}

AMEGIC::Process_Base * Broker_Group::Selected() { 
  if (selected==0)    return 0; 
  if (selected==this) return this;
  return selected->Selected(); 
}    

double Broker_Group::Differential(AMATOOLS::Vec4D * p) 
{
  last = 0.;
  for (int i=0;i<brokers.size();i++) last += brokers[i]->Differential(p);
  if ((!(last<=0)) && (!(last>0))) {
    msg.Error()<<"--- Broker_Group::Differential ---"<<std::endl;
  }
  return last;
};

double Broker_Group::Differential2() { 
  double tmp = 0.;
  for (int i=0;i<brokers.size();i++) tmp += brokers[i]->DSigma2();
  if ((!(tmp<=0)) && (!(tmp>0))) {
    msg.Error()<<"--- Broker_Group::Differential ---"<<std::endl;
  }
  last += tmp;
  return tmp;
}

double Broker_Group::DSigma(AMATOOLS::Vec4D * p, bool lookup) 
{ 
  last = 0;
  for (int i=0;i<brokers.size();i++) last += brokers[i]->DSigma(p, lookup);
  return last;
}

double Broker_Group::DSigma2()
{
  double tmp = 0.;
  for (int i=0;i<brokers.size();i++) tmp += brokers[i]->DSigma2();
  last += tmp;
  return tmp;
}

void Broker_Group::SetISR(ISR::ISR_Handler * _isr) { 
  msg.Debugging()<<"Broker_Group::SetISR("<<_isr->Type()<<") for "<<name<<"  : "<<_isr<<std::endl;
  isr = _isr; 

//   isr->SetSprimeMin(taumin*s);
//   isr->SetSprimeMax(taumax*s);

  for (int i=0;i<brokers.size();i++) brokers[i]->SetISR(isr);
}

void Broker_Group::SetBeam(BEAM::Beam_Handler * _beam) { 
  msg.Debugging()<<"Broker_Group::SetBeam("<<_beam->Type()<<") for "<<name<<"  : "<<_beam<<std::endl;
  beam = _beam; 

//   isr->SetSprimeMin(taumin*s);
//   isr->SetSprimeMax(taumax*s);

  for (int i=0;i<brokers.size();i++) brokers[i]->SetBeam(beam);
}

bool Broker_Group::CalculateTotalXSec() {
  msg.Tracking()<<"Broker_Group::CalculateTotalXSec()"<<std::endl;
  if (atoms) {
    bool okay = 1;
    for (int i=0;i<brokers.size();i++) {
      msg.Tracking()<<"Broker_Group::CalculateTotalXSec for "<<brokers[i]->Name()<<std::endl;
      if (!(brokers[i]->CalculateTotalXSec())) okay = 0;
    }
    return okay;
  }
  else {
    if (nin==2) {
      if ( (fl[0].Mass() != rpa.gen.Beam1().Mass()) ||
	   (fl[1].Mass() != rpa.gen.Beam2().Mass()) ) isr->SetPartonMasses(fl);
    }
    sel->BuildCuts(cuts);
    totalxs = ps->Integrate()/AORGTOOLS::rpa.Picobarn(); 
    if (!(AMATOOLS::IsZero((n*totalxs-totalsum)/(n*totalxs+totalsum)))) {
      msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<std::endl
		 <<"  "<<name<<" : "<<totalxs<<" vs. "<<totalsum/n<<std::endl;
    }
    SetTotalXS(0);
    if (totalxs>0.) return 1;
  }
}

bool Broker_Group::SetUpIntegrator(ISR::ISR_Handler * _isr,BEAM::Beam_Handler * _beam) {
  SetISR(_isr);
  SetBeam(_beam);
  if (atoms) {
    for (int i=0; i<brokers.size(); i++) brokers[i]->SetUpIntegrator(isr, beam);
    return 1;
  }
  sel->BuildCuts(cuts);
  if (nin==2) {
    if ( (fl[0].Mass() != rpa.gen.Beam1().Mass()) ||
	 (fl[1].Mass() != rpa.gen.Beam2().Mass()) ) isr->SetPartonMasses(fl);
  }
  ps  = new Phase_Space_Handler(this,_isr,_beam);
  if (!ps->CreateIntegrators()) return 0;
  return 1;
}

bool Broker_Group::OneEvent() {
  if (atoms) {
    SelectOne();
    return selected->OneEvent();
  }
  else return ps->OneEvent();
}

bool Broker_Group::SameEvent() {
  if (atoms) {
    if (selected)
      return selected->SameEvent();
    msg.Error()<<" ERROR in bool Process_Group::SameEvent() "<<std::endl;
    return 0;
  }
  else return ps->SameEvent();
}

double Broker_Group::WeightedEvent() {
  if (atoms) {
    SelectOne();
    return selected->WeightedEvent();
  }
  else return ps->WeightedEvent();
}

void Broker_Group::SetScale(double _scale)
{
  scale = _scale;
  for (int i=0;i<brokers.size();i++) brokers[i]->SetScale(scale); 
} 

double Broker_Group::KFactor(double _scale) 
{
  return brokers[0]->KFactor(_scale); 
}
