#include "Beam_Handler.H"
#include "No_Beam.H"
#include "Laser_Backscattering.H"
#include "Run_Parameter.H" 
#include "Message.H"
#include <stdio.h>


using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace BEAM;
using namespace std;

/*
std::ostream& operator<<(std::ostream &,AMEGIC::Pol_Info&);


std::ostream& operator<<(std::ostream& str,AMEGIC::Pol_Info& ) {
    str<<"      "<<part->momentum()<<" "<<part->momentum().abs2()
       <<" colours : ("<<part->flow(1)<<","<<part->flow(2)<<")"<<std::endl;
    return str;
}
*/

Beam_Handler::Beam_Handler(int * beamtypes,Flavour * bunches,double * plbunches,
			   Flavour * beams,double * plbeams, 
			   double * _splimits)
{

  mode = 0;
  BeamBase  = new Beam_Base*[2];
  for (short int i=0;i<2;i++) {
    switch (beamtypes[i]) {
    case Beam_Type::No : 
      cout<<" new No_Beam_at_all "<<endl;
      if (plbunches)  BeamBase[i] = new No_Beam_at_all(bunches[i],plbunches[i]);
      else  BeamBase[i] = new No_Beam_at_all(bunches[i],0.);
      break; 
      cout<<" done."<<endl;
    case Beam_Type::Laser_Back : 
	
      if (bunches && plbunches)      msg.Out()<<"In Beam_Handler : "<<i<<":"<<bunches[i]<<" with polarization: "<<plbunches[i]<<endl;
      BeamBase[i] = new Laser_Backscattering(bunches[i],plbunches[i]);
      mode       += i+1;
      break;
    default:
      msg.Error()<<"No Beam found for beam ("<<i+1<<"), type = "
		 <<beamtypes[i]<<". Initialize No Beam."<<endl;
      BeamBase[i] = new No_Beam_at_all(bunches[i],plbunches[i]);break;
    }
  }
  type = BeamBase[0]->Type() + std::string("*") + BeamBase[1]->Type();

  double s    = sqr(AORGTOOLS::rpa.gen.Ecms());
  splimits[0] = smin = _splimits[0];
  splimits[1] = smax = AMATOOLS::Min(_splimits[1],s*Upper1()*Upper2());
  splimits[2] = s;
  ylimits[0]  = -10.;
  ylimits[1]  = 10.;

  exponent[0] = .5;
  exponent[1] = .98 * ( BeamBase[0]->Exponent() + BeamBase[1]->Exponent());

  SetBeamMasses(beams);

  if (mode>0) msg.Debugging()<<"Beam is on;  "<<splimits[0]<<" .. "<<splimits[1]<<" .. "<<splimits[2]<<" .. "<<endl;
         else msg.Debugging()<<"Beam is off; ";
  msg.Debugging()<<"type = "<<type<<" for "<<bunches[0]<<" / "<<bunches[1]<<endl;
}

bool Beam_Handler::CheckConsistency(APHYTOOLS::Flavour * _bunches,
				    APHYTOOLS::Flavour * _beams) {
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (BeamBase[i]->Type() == string("Laser_Backscattering")) {
      if (! ( ((_bunches[i] == Flavour(kf::e)) || (_bunches[i] == Flavour(kf::e).bar())) &&
	      (_beams[i] == Flavour(kf::photon))         ) ) {
	fit = 0;
	break;
      }
    }
    if (BeamBase[i]->Type() == string("Beam_Strahlung")) {
      if (! ( ((_bunches[i] == Flavour(kf::e)) || (_bunches[i] == Flavour(kf::e).bar())) &&
	      (_beams[i] == _bunches[i])         ) ) {
	fit = 0;
	break;
      }
    }
  }
  return fit;
}


void Beam_Handler::SetBeamMasses(Flavour * beams) {
  mass12      = sqr(beams[0].mass());
  mass22      = sqr(beams[1].mass());
  double E    = AORGTOOLS::rpa.gen.Ecms();
  double x    = 1./2.+(mass12-mass22)/(2.*E*E);
  double E1   = x*E;
  double E2   = E-E1;
  fixvecs[0]  = vec4d(E1,0.,0., sqrt(sqr(E1)-mass12));
  fixvecs[1]  = vec4d(E2,0.,0.,-sqrt(sqr(E1)-mass12));
}



Beam_Handler::~Beam_Handler() { }

bool Beam_Handler::MakeBeams(vec4d * p,double sprime,double y) 
{
  if (mode==0) {
    x1   = x2 = 1.;
    p[0] = fixvecs[0];
    p[1] = fixvecs[1];
    return 1;
  }
  else {
    if ( (sprime<splimits[0]) || (sprime>splimits[1]) ) {
	//msg.Debugging()<<"Beam_Handler::MakeBeam("<<sprime<<", "<<y<<") out of bounds!"<<std::endl
	//	     <<"   "<<splimits[0]<<" < "<<sprime<<" < "<<splimits[1]<<std::endl;
      return 0;
    }

    double E      = sqrt(splimits[2]);
    double Eprime = sqrt(sprime);
    double x      = 1./2.+(mass12-mass22)/(2.*sprime);
    double E1     = x*Eprime;
    double E2     = Eprime-E1;
    
    // initial state momenta in CMS frame
    p[0]          = vec4d(E1,0.,0.,sqrt(sqr(E1)-mass12));
    p[1]          = vec4d(E2,(-1.)*vec3d(p[0]));
    E1            = std::exp(y);  
    E2            = std::exp(-y);  

    // establish boost
    CMSBoost = Poincare(vec4d(E1+E2,0.,0.,E1-E2));
    
    // calculate real x1,2
    vec4d p1 = p[0];
    vec4d p2 = p[1];
    CMSBoost.BoostBack(p1);
    CMSBoost.BoostBack(p2);
    x1       = 2.*p1[0]/E;
    x2       = 2.*p2[0]/E;

    //AORGTOOLS::msg.Out()<<"*** x1/2 : "<<x1<<" , "<<x2<<endl;
    if (mode==1) {
	//msg.Out()<<"Mode = 1 : x2 should be 1, is "<<x2<<endl;
      x2 = 1.;
    }
    if (mode==2) {
	//msg.Out()<<"Mode = 2 : x1 should be 1, is "<<x1<<endl;
      x1 = 1.;
    }

    //msg.Out()<<p1<<"/"<<p2<<endl;
    return 1;
  }
}

/* ----------------------------------------------------------------

   Weight calculation 

   ---------------------------------------------------------------- */


bool Beam_Handler::CalculateWeight(double scale) 
{
  switch (mode) {
  case 3 :
    if ( (BeamBase[0]->CalculateWeight(x1,scale)) && 
	 (BeamBase[1]->CalculateWeight(x2,scale)) ) return 1;
    break;
  case 2 :
    if (BeamBase[1]->CalculateWeight(x2,scale))     return 1;
    break;
  case 1 :
    if (BeamBase[0]->CalculateWeight(x1,scale))     return 1;
    break;
  }
  return 0;
};


double Beam_Handler::Weight(Flavour * flin)
{
  double weight = (BeamBase[0]->Weight(flin[0]) * BeamBase[1]->Weight(flin[1]));
  //cout<<"BH : "<<weight<<endl;
  return weight;
}

double Beam_Handler::Weight()
{
  double weight = (BeamBase[0]->Weight() * BeamBase[1]->Weight());
  //cout<<"BH : "<<weight<<endl;
  return weight;
}

/* ----------------------------------------------------------------

   Boosts

   ---------------------------------------------------------------- */


void  Beam_Handler::BoostInCMS(vec4d* p,int n) {
  for (int i=0; i<n; ++i) CMSBoost.Boost(p[i]);
}

void  Beam_Handler::BoostInLab(vec4d* p,int n) {
  for (int i=0; i<n; ++i) CMSBoost.BoostBack(p[i]);
}
