#include "ISR_Handler.H"
#include "Intact.H"
#include "Structure_Function.H"
#include "Run_Parameter.H" 
#include "Message.H"
#include <stdio.h>


using namespace ATOOLS;
using namespace PDF;
using namespace std;

ISR_Handler::ISR_Handler(ISR_Base ** _ISRBase,double * _splimits) :
  p_ISRBase(_ISRBase), m_mass12(0.), m_mass22(0.), m_x1(1.), m_x2(1.)
{
  m_mode = 0;
  for (short int i=0;i<2;i++) {
    if (p_ISRBase[i]->On()) m_mode += i+1;
  }
  m_mass12     = sqr(p_ISRBase[0]->Flav().Mass());
  m_mass22     = sqr(p_ISRBase[1]->Flav().Mass());
  Init(_splimits);
}

void ISR_Handler::Init(double * _splimits) {
  m_type = p_ISRBase[0]->Type() + std::string("*") + p_ISRBase[1]->Type();

  double s      = sqr(ATOOLS::rpa.gen.Ecms());
  m_splimits[0] = s*_splimits[0];
  m_splimits[1] = ATOOLS::Min(s*_splimits[1],s*Upper1()*Upper2());
  m_splimits[2] = s;
  m_fixed_smin = m_splimits[0];
  m_fixed_smax = m_splimits[1];
  m_ylimits[0]  = -10.;
  m_ylimits[1]  = 10.;
  m_exponent[0] = .5;
  m_exponent[1] = .98 * p_ISRBase[0]->Exponent() * p_ISRBase[1]->Exponent();

  if (m_mode>0) msg.Tracking()<<"ISR is on:  ";
           else msg.Tracking()<<"ISR is off: ";
  msg.Tracking()<<"type = "<<m_type<<" for "
		 <<p_ISRBase[0]->Flav()<<" / "<<p_ISRBase[1]->Flav()<<endl
		 <<"            Range = "<<m_splimits[0]<<" ... "<<m_splimits[1]<<" from "<<m_splimits[2]<<endl;
}


void   ISR_Handler::SetSprimeMin(double _spl)       { m_splimits[0]  = Max(m_fixed_smin,_spl); }
void   ISR_Handler::SetSprimeMax(double _spl)       { m_splimits[1]  = Min(m_fixed_smax,_spl); }
void   ISR_Handler::SetFixedSprimeMin(double _spl)  
{ 
  m_fixed_smin  = Max(m_fixed_smin,_spl);
  m_splimits[0] = Max(m_splimits[0],_spl);
}
void   ISR_Handler::SetFixedSprimeMax(double _spl)  
{
  m_fixed_smax  = Min(m_fixed_smax,_spl);
  m_splimits[1] = Min(m_splimits[1],_spl);
}


ISR_Handler::~ISR_Handler() {
  if (p_ISRBase) {
    for (int i=0;i<2;i++) {
      if (p_ISRBase[i]) delete p_ISRBase[i];  
    }
    delete[] p_ISRBase; p_ISRBase = 0;
  }
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour * _bunches,
				   ATOOLS::Flavour * _partons) {
  
    bool fit = 1;
  for (int i=0;i<2;i++) {
    if (p_ISRBase[i]->On()) {
      if (_bunches[i] != PDF(i)->Bunch()) { fit = 0; break; }
      fit = 0;
      for (unsigned int j = 0;j<(PDF(i)->Partons()).size();j++) {
	if (_partons[i] == (PDF(i)->Partons())[j]) {
	  fit = 1;
	  break; 
	}
      }
      if (fit == 0) break;
    }
    else {
      if (_partons[i]!=p_ISRBase[i]->Flav()) {
	fit = 0;
	break;
      }
    }
  }
  return fit;
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour * _partons) {
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (p_ISRBase[i]->On()) {
      fit = 0;
      for (unsigned int j = 0;j<(PDF(i)->Partons()).size();j++) {
	if (_partons[i] == (PDF(i)->Partons())[j]) {
	  fit = 1;
	  break; 
	}
      }
      if (fit == 0) break;
    }
    else {
      if (_partons[i]!=p_ISRBase[i]->Flav()) {
	fit = 0;
	break;
      }
    }
  }
  return fit;
}

void ISR_Handler::SetPartonMasses(Flavour * fl) {
  m_mass12     = sqr(fl[0].Mass());
  m_mass22     = sqr(fl[1].Mass());
  double E     = ATOOLS::rpa.gen.Ecms();
  double x     = 1./2.+(m_mass12-m_mass22)/(2.*E*E);
  double E1    = x*E;
  double E2    = E-E1;
  m_fiXVECs[0] = Vec4D(E1,0.,0., sqrt(sqr(E1)-m_mass12));
  m_fiXVECs[1] = Vec4D(E2,0.,0.,-sqrt(sqr(E1)-m_mass12));
}


bool ISR_Handler::MakeISR(Vec4D * p,double sprime,double y) 
{
  if (m_mode==0) {
    m_x1 = m_x2 = 1.;
    p[0] = m_fiXVECs[0];
    p[1] = m_fiXVECs[1];
    return 1;
  }
  else {
    if ( (sprime<m_splimits[0]) || (sprime>m_splimits[1]) ) {
      msg.Error()<<"MakeISR : sprime out of bounds !!!"<<endl
		 <<"   "<<m_splimits[0]<<"<"<<sprime<<"<"<<m_splimits[1]<<"<"<<m_splimits[2]<<endl;
      return 0;
    }

    double E      = sqrt(m_splimits[2]);
    double Eprime = sqrt(sprime);
    double x      = 1./2.+(m_mass12-m_mass22)/(2.*sprime);
    double E1     = x*Eprime;
    double E2     = Eprime-E1;
    p[0]          = Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_mass12));
    p[1]          = Vec4D(E2,(-1.)*Vec3D(p[0]));

    E1            = exp(y);  
    E2            = exp(-y);  

    m_CMSBoost    = Poincare(Vec4D(E1+E2,0.,0.,E1-E2));
    
    Vec4D p1      = p[0];
    Vec4D p2      = p[1];
    m_CMSBoost.BoostBack(p1);
    m_CMSBoost.BoostBack(p2);
    m_x1          = 2.*p1[0]/E;
    m_x2          = 2.*p2[0]/E;

    if (m_mode==1) m_x2 = 1.;
    if (m_mode==2) m_x1 = 1.;
  
    return 1;
  }
}

/* ----------------------------------------------------------------

   Weight calculation 

   ---------------------------------------------------------------- */


bool ISR_Handler::CalculateWeight(double scale) 
{
  switch (m_mode) {
  case 3 :
    if ( (p_ISRBase[0]->CalculateWeight(m_x1,scale)) && 
	 (p_ISRBase[1]->CalculateWeight(m_x2,scale)) ) return 1;
    break;
  case 2 :
    if (p_ISRBase[1]->CalculateWeight(m_x2,scale))     return 1;
    break;
  case 1 :
    if (p_ISRBase[0]->CalculateWeight(m_x1,scale))     return 1;
    break;
  }
  return 0;
};

bool ISR_Handler::CalculateWeight2(double scale) 
{
  if (m_mode != 3) { 
    msg.Error()<<"ISR_Handler::CalculateWeight2 called for one ISR only."<<endl;
    abort();
  }
  if ( (p_ISRBase[0]->CalculateWeight(m_x2,scale)) && 
       (p_ISRBase[1]->CalculateWeight(m_x1,scale)) ) { 
    return 1;
  }
  return 0;
};

double ISR_Handler::Weight(Flavour * flin)
{
  if (m_mode!=3 || (CheckRemnantKinematics(flin[0],m_x1,0) &&
      CheckRemnantKinematics(flin[1],m_x2,1))) 
    return (p_ISRBase[0]->Weight(flin[0]) * p_ISRBase[1]->Weight(flin[1]));
  return 0.;
}

double ISR_Handler::Weight2(Flavour* flin)
{
  if (CheckRemnantKinematics(flin[0],m_x1,1) &&
      CheckRemnantKinematics(flin[1],m_x2,0)) 
    return (p_ISRBase[0]->Weight(flin[1]) * p_ISRBase[1]->Weight(flin[0]));
  return 0.;
}



/* ----------------------------------------------------------------

   Boosts

   ---------------------------------------------------------------- */


void  ISR_Handler::BoostInCMS(Vec4D* p,int n) {
  for (int i=0; i<n; ++i) m_CMSBoost.Boost(p[i]);
}

void  ISR_Handler::BoostInLab(Vec4D* p,int n) {
  for (int i=0; i<n; ++i) m_CMSBoost.BoostBack(p[i]);
}


const Flavour ISR_Handler::DiQuark(const Flavour & fl1,const Flavour & fl2) 
{
  // lightes flavour with that content
  int kf1=fl1.Kfcode();
  int kf2=fl2.Kfcode();
  Flavour diquark;
  if (kf1>kf2) diquark =(kf::code)(kf1*1000 + kf2*100 + 1);
  else if (kf1<kf2)  diquark =(kf::code)(kf2*1000 + kf1*100 + 1);
  else diquark =(kf::code)(kf2*1000 + kf1*100 + 3);
  if (fl1.IsAnti()) diquark=diquark.Bar();
  return diquark;    
}


bool ISR_Handler::CheckRemnantKinematics(const ATOOLS::Flavour & fl,double x,int nbeam)
{
  if (x<.99) return true;

  double mf   = fl.PSMass();
  double msum = 0.;

  double erem = (1. -x -1.e-6)*ATOOLS::rpa.gen.Ecms();
  
  Flavour bunch=p_ISRBase[nbeam]->Flav();
  if (!bunch.IsHadron()) return true;
  
  int hadint=(bunch.Kfcode()-bunch.Kfcode()/10000)/10;

  if ((hadint<=100)||(hadint>=1000)) return true;

  Flavour constit[3];
  constit[0]=ATOOLS::Flavour(ATOOLS::kf::code(hadint)/100);
  constit[1]=ATOOLS::Flavour(ATOOLS::kf::code((hadint-(hadint/100)*100)/10));
  constit[2]=ATOOLS::Flavour(ATOOLS::kf::code(hadint-(hadint/10)*10));
  if (bunch.IsAnti()) {
    for (int i=0;i<3;i++) constit[i]=constit[i].Bar();
  }

  // valence quark
  for (int i=0;i<3;++i) if (constit[i]==fl) {
    for (int j=i+1;j<3;++j) {
      constit[j-1]=constit[j];
    }
    Flavour diquark=DiQuark(constit[0], constit[1]);
    msum+=diquark.PSMass();
    break;
  }

  if (msum==0) {
    // gluon
    if (fl.IsGluon()) {
      Flavour diquark=DiQuark(constit[0], constit[2]);
      msum+=constit[1].PSMass()+diquark.PSMass();    
    }
    else {
      // sea quark
      Flavour diquark=DiQuark(constit[0], constit[2]);
      msum+=constit[1].PSMass()+diquark.PSMass();
      msum+=fl.PSMass();
    }
  }

  if (erem < mf + msum) return false;
  return true;
}
