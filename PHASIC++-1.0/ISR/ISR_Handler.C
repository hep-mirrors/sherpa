#include "ISR_Handler.H"
#include "No_ISR.H"
#include "Structure_Function.H"
#include "Run_Parameter.H" 
#include "Message.H"
#include <stdio.h>


using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace ISR;
using namespace std;

ISR_Handler::ISR_Handler(int * isrtypes,Flavour * beams,Flavour * partons,
			 double * _splimits)
{
  ISRBase  = new ISR_Base*[2];
  mode = 0;
  for (short int i=0;i<2;i++) {
    switch (isrtypes[i]) {
    case ISR_Type::No : 
      ISRBase[i] = new No_ISR_at_all(beams[i]);break;
    case ISR_Type::Extended_Struc : 
      ISRBase[i] = new Structure_Function(beams[i]);
      mode      += i+1;
      break;
    default:
      msg.Error()<<"No ISR found for beam ("<<i+1<<"). Initialize No ISR."<<endl;
      ISRBase[i] = new No_ISR_at_all(beams[i]);break;
    }
  }
  type = ISRBase[0]->Type() + std::string("*") + ISRBase[1]->Type();

  double s    = sqr(AORGTOOLS::rpa.gen.Ecms());
  splimits[0] = smin = _splimits[0];
  splimits[1] = smax = AMATOOLS::Min(_splimits[1],s*Upper1()*Upper2());
  splimits[2] = s;

//   cout<<"splimits[0]  "<<splimits[0]<<endl;
//   cout<<"splimits[1]  "<<splimits[1]<<endl;
//   cout<<"splimits[2]  "<<splimits[2]<<endl;
//   cout<<"upper1  "<<Upper1()<<endl;
//   cout<<"upper2  "<<Upper2()<<endl;

  ylimits[0]  = -10.;
  ylimits[1]  = 10.;
  exponent[0] = .5;
  exponent[1] = .98 * ISRBase[0]->Exponent() * ISRBase[1]->Exponent();

  SetPartonMasses(partons);


  if (mode>0) msg.Debugging()<<"ISR is on;  ";
         else msg.Debugging()<<"ISR is off; ";
  msg.Debugging()<<"type = "<<type<<" for "<<beams[0]<<" / "<<beams[1]<<endl;
}

ISR_Handler::~ISR_Handler() {
  if (ISRBase) {
    for (int i=0;i<2;i++) {
      if (ISRBase[i]) delete ISRBase[i]; 
    }
    delete[] ISRBase; ISRBase = 0;
  }
}

bool ISR_Handler::CheckConsistency(APHYTOOLS::Flavour * _beams,
				   APHYTOOLS::Flavour * _partons) {
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (ISRBase[i]->On()) {
      if (_beams[i] != PDF(i)->Beam()) { fit = 0; break; }
      fit = 0;
      for (int j = 0;j<(PDF(i)->Partons()).size();j++) {
	if (_partons[i] == (PDF(i)->Partons())[j]) {
	  fit = 1;
	  break; 
	}
      }
      if (fit == 0) break;
    }
  }
  return fit;
}

void ISR_Handler::SetPartonMasses(Flavour * _fl) { 
  mass12      = sqr(_fl[0].Mass());
  mass22      = sqr(_fl[1].Mass());
  double E    = AORGTOOLS::rpa.gen.Ecms();
  double x    = 1./2.+(mass12-mass22)/(2.*E*E);
  double E1   = x*E;
  double E2   = E-E1;
  fiXVECs[0]  = Vec4D(E1,0.,0., sqrt(sqr(E1)-mass12));
  fiXVECs[1]  = Vec4D(E2,0.,0.,-sqrt(sqr(E1)-mass12));
}

bool ISR_Handler::MakeISR(Vec4D * p,double sprime,double y) 
{
  if (mode==0) {
    x1   = x2 = 1.;
    p[0] = fiXVECs[0];
    p[1] = fiXVECs[1];
    return 1;
  }
  else {
    if ( (sprime<splimits[0]) || (sprime>splimits[1]) ) {
      msg.Error()<<"MakeISR : sprime out of bounds !!!"<<endl
		 <<"   "<<splimits[0]<<"<"<<sprime<<"<"<<splimits[1]<<"<"<<splimits[2]<<endl;
      return 0;
    }

    double E      = sqrt(splimits[2]);
    double Eprime = sqrt(sprime);
    double x      = 1./2.+(mass12-mass22)/(2.*sprime);
    // Energies in the c.m. system
    double E1     = x*Eprime;
    double E2     = Eprime-E1;
    
    // initial state momenta in CMS frame
    p[0]          = Vec4D(E1,0.,0.,sqrt(sqr(E1)-mass12));
    p[1]          = Vec4D(E2,(-1.)*Vec3D(p[0]));
    // Energies in the lab system
    E1            = exp(y);  
    E2            = exp(-y);  

    // establish boost
    CMSBoost = Poincare(Vec4D(E1+E2,0.,0.,E1-E2));
    
    // calculate real x1,2
    Vec4D p1 = p[0];
    Vec4D p2 = p[1];
    CMSBoost.BoostBack(p1);
    CMSBoost.BoostBack(p2);
    x1       = 2.*p1[0]/E;
    x2       = 2.*p2[0]/E;

//     AORGTOOLS::msg.Debugging()<<"sprime,y,x1,x2 = "<<sprime<<" , "<<y<<" , "<<x1<<" , "<<x2<<std::endl;

    if (mode==1) {
      x2 = 1.;
    }
    if (mode==2) {
      x1 = 1.;
    }
    return 1;
  }
}

/* ----------------------------------------------------------------

   Weight calculation 

   ---------------------------------------------------------------- */

bool ISR_Handler::CalculateWeight(double scale) 
{
  switch (mode) {
  case 3 :
//     msg.Debugging()<<" x1 : "<<x1<<endl;
//     msg.Debugging()<<" x2 : "<<x2<<endl;
    if ( (ISRBase[0]->CalculateWeight(x1,scale)) && 
	 (ISRBase[1]->CalculateWeight(x2,scale)) ) return 1;
    break;
  case 2 :
    if (ISRBase[1]->CalculateWeight(x2,scale))     return 1;
    break;
  case 1 :
    if (ISRBase[0]->CalculateWeight(x1,scale))     return 1;
    break;
  }
  return 0;
};

bool ISR_Handler::CalculateWeight2(double scale) 
{
  if (mode != 3) { 
    msg.Error()<<"ISR_Handler::CalculateWeight2 called for one ISR only."<<endl;
    abort();
  }
  if ( (ISRBase[0]->CalculateWeight(x2,scale)) && 
       (ISRBase[1]->CalculateWeight(x1,scale)) ) { 
    return 1;
  }
  return 0;
};

double ISR_Handler::Weight(Flavour * flin)
{
//   AORGTOOLS::msg.Debugging()<<"weight  = "<<ISRBase[0]->Weight(flin[0]) * ISRBase[1]->Weight(flin[1])<<std::endl;
  return (ISRBase[0]->Weight(flin[0]) * ISRBase[1]->Weight(flin[1]));
}

double ISR_Handler::Weight2(Flavour* flin)
{
//   AORGTOOLS::msg.Debugging()<<"weight2 = "<<ISRBase[0]->Weight(flin[1]) * ISRBase[1]->Weight(flin[0])<<std::endl;
  return (ISRBase[0]->Weight(flin[1]) * ISRBase[1]->Weight(flin[0]));
}

/* ----------------------------------------------------------------

   Boosts

   ---------------------------------------------------------------- */

void  ISR_Handler::BoostInCMS(Vec4D* p,int n) {
  for (int i=0; i<n; ++i) CMSBoost.Boost(p[i]);
}

void  ISR_Handler::BoostInLab(Vec4D* p,int n) {
  for (int i=0; i<n; ++i) CMSBoost.BoostBack(p[i]);
}
