#include "Polarisation.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"

using namespace AMEGIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace std;

Polarisation::Polarisation()
{
  nmass     = 0;
  Mass_Norm = 1.;
  npol      = 0;
  mass_pol  = 0;
  no        = 0;
}

Polarisation::~Polarisation()
{
  if (mass_pol!=0) {
    for (short int i=0;i<no;i++) delete[] mass_pol[i];
    delete[] mass_pol;
  }
}


double Polarisation::Spin_Average(int nin,Flavour* flin)
{
  //including colours
  //unpolarized
  double Norm = 1.;
  for (short int i=0;i<nin;i++) {
    if (flin[i].isquark())   Norm *= 3.;    
    if (flin[i].issquark())  Norm *= 3.;// sQuark is boson !!    
    if (flin[i].isgluon())   Norm *= 8.;
    if (flin[i].isgluino())  Norm *= 8.;// Gluino is fermion !!

    if (flin[i].isfermion()) Norm *= 2.;
    if (flin[i].isvector()) {
      if (AMATOOLS::IsZero(flin[i].mass())) Norm *= 2.;
                                       else Norm *= 3.;
    }
  }
  return 1./Norm;
}

void Polarisation::Add_Extern_Polarisations(Basic_Sfuncs* BS,Flavour* fl,Helicity *hel)
{
#ifdef Explicit_Pols
  for(short int i=0;i<BS->GetNmomenta();i++)
    if(fl[i].isvector())BS->Build_Polarisations(i,hel->p_type[i],hel->angle[i]);
#endif
}

int Polarisation::Massless_Vectors(int N,Flavour* fl)
{
#ifndef Explicit_Pols
  for(short int i=0;i<N;i++) {
    if (fl[i].isvector() && AMATOOLS::IsZero(fl[i].mass())) {
      npol = 1;
      break;
    }
  }
  return npol;
#else
  return 0;
#endif
}

int Polarisation::Massive_Vectors(int N,Flavour* fl)
{
#ifndef Explicit_Pols
  int nmass_old = nmass;
  for(short int i=0;i<N;i++) {
    if (fl[i].isvector() && !AMATOOLS::IsZero(fl[i].mass())) {
      nmass+=2;
      AORGTOOLS::msg.Debugging()<<"Mass_Norm changed for : "<<fl[i]<<endl;
      Mass_Norm *= 3./2./sqr(fl[i].mass());//4 Pi for Phasespace
    }
  } 
  return nmass-nmass_old;
#else
  return 0;
#endif
}  

void Polarisation::Attach(int N, Flavour* fl)
{
#ifndef Explicit_Pols
  if (nmass>0) {
    mass_pol = new int*[N];
    no = N;
    for (short int i=0;i<N;i++) mass_pol[i] = new int[2];
    int count = N+npol;
    for(short int i=0;i<N;i++) {
      if (fl[i].isvector() && !AMATOOLS::IsZero(fl[i].mass())) {
	mass_pol[i][0] = count;
	mass_pol[i][1] = count+1;
	count+=2;
      }
      else {
	mass_pol[i][0] = -1; 
	mass_pol[i][1] = -1;
      } 
    } 
  }
#endif
}

void Polarisation::Reset_Gauge_Vectors(int N,vec4d* p,vec4d gauge)
{
#ifndef Explicit_Pols
  if (npol==1) p[N] = vec4d(vec3d(gauge).abs(),vec3d(gauge));
#endif
}

void Polarisation::Set_Gauge_Vectors(int N,vec4d* p,vec4d gauge)
{
#ifndef Explicit_Pols
  if (npol==1) p[N] = vec4d(vec3d(gauge).abs(),vec3d(gauge));
  if (nmass>0) {
    for (short int i=0;i<N;i++) {
      if (mass_pol[i][0]>=0) {
	vec4d r1,r2;
        vec4d pisave = p[i];
	double mass = sqrt(p[i].abs2());	
	double C = 2.*Ran.get()-1.;
	double S = sqrt(1.-C*C);
	double F = 2.*M_PI*Ran.get();
	r1 = mass/2.*vec4d(1.,S*::sin(F),S*::cos(F),C);
	r2 = vec4d(r1[0],(-1.)*vec3d(r1));
	vec4d help;
	//r2 boost
	help[0] = (p[i][0]*r2[0]+vec3d(p[i])*vec3d(r2))/mass;
	double c1 = (r2[0]+help[0])/(mass+p[i][0]);
	p[mass_pol[i][0]] = vec4d(help[0],vec3d(r2)+c1*vec3d(p[i]));  
	//r1 boost
	help[0] = (p[i][0]*r1[0]+vec3d(p[i])*vec3d(r1))/mass;
        c1 = (r1[0]+help[0])/(mass+p[i][0]);
	p[mass_pol[i][1]] = vec4d(help[0],vec3d(r1)+c1*vec3d(p[i]));
      }
    }
  }  
#endif
}

double Polarisation::Massless_Norm(int N,Flavour* fl,Basic_Sfuncs* BS)
{
#ifndef Explicit_Pols
  double norm = 1.;
  if (npol==1) {
    for (short int i=0;i<N;i++) {
      if (fl[i].isvector() && AMATOOLS::IsZero(fl[i].mass()) ) {
        for (short int j=i+1;j<N+1;j++) {
	  if ((fl[j].isvector() && AMATOOLS::IsZero(fl[j].mass()) ) || 
	      (fl[j]==Flavour(kf::pol))) {
	    norm *= BS->N(i,j);
	    break;
	  }
	}
      }
    }
  }
  return norm;
#else
  return 1.;
#endif
}

void Polarisation::Replace_Numbers(int N,Flavour* fl,Single_Amplitude* n)
{
#ifndef Explicit_Pols
  //Test!!!!!!!!!!!!!!!!!!!!!
  //n->MPolconvert(12,10);
  //n->MPolconvert(14,12);

  if (nmass>0) {
    for (short int i=0;i<N;i++) {
      if (mass_pol[i][0]>=0) {
	n->MPolconvert(i,mass_pol[i][0]+20);
	n->MPolconvert(i+20,mass_pol[i][1]+20);
	n->Prop_Replace(fl[i],i,mass_pol[i][0],mass_pol[i][1]);
      }
    }
  }
  //incoming massless bosons
  for (short int i=0;i<N;i++) {
    if (fl[i].isvector() && AMATOOLS::IsZero(fl[i].mass())) {
      if (!( (fl[i+1].isvector() && AMATOOLS::IsZero(fl[i+1].mass())) ||
	     fl[i+1]==Flavour(kf::pol) ) ) {
	//search next boson or pol
	for (short int j=i+1;j<N+1;j++) {
	  if ( (fl[j].isvector() && AMATOOLS::IsZero(fl[j].mass())) ||
	       fl[j]==Flavour(kf::pol) ) {
	    n->MPolconvert(i+10+1,j+10);
	    break;
	  }
	}	      
	
      }
    }
  }
#else
  for (short int i=0;i<N;i++) {
    if (fl[i].isvector() && AMATOOLS::IsZero(fl[i].mass()))  n->MPolconvert(i+10+1,99);
    if (fl[i].isvector() && !AMATOOLS::IsZero(fl[i].mass())) n->MPolconvert(i+20,99);
  }  
#endif
}
  

