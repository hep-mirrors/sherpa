#include "AMEGIC++/Main/Polarisation.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AMEGIC;
using namespace ATOOLS;
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
  CSC.Init();
  for (int i=0;i<nin;i++) {
    if (abs(flin[i].StrongCharge())==3)  Norm *= CSC.Nc; 
    if (flin[i].StrongCharge()==8)       Norm *= CSC.Nc*CSC.Nc-1.; 

    if (flin[i].IsFermion()) Norm *= 2.;
    if (flin[i].IsVector()) {
      if (ATOOLS::IsZero(flin[i].Mass())) Norm *= 2.;
                                       else Norm *= 3.;
    }
    if (flin[i].IsTensor()) Norm *= 5.;
  }
  return 1./Norm;
}

void Polarisation::Add_Extern_Polarisations(Basic_Sfuncs* BS,Flavour* fl,Helicity *hel)
{
#ifdef Explicit_Pols
  for(short int i=0;i<BS->GetNmomenta();i++){
    if(fl[i].IsVector())BS->BuildPolarisations(i,hel->PolTypes()[i],hel->PolAngles()[i]);
    if(fl[i].IsTensor())BS->BuildTensorPolarisations(i);
  } 
#endif
}

int Polarisation::Massless_Vectors(int N,Flavour* fl)
{
#ifndef Explicit_Pols
  for(short int i=0;i<N;i++) {
    if (fl[i].IsVector() && ATOOLS::IsZero(fl[i].Mass())) {
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
    if (fl[i].IsVector() && !ATOOLS::IsZero(fl[i].Mass())) {
      nmass+=2;
      Mass_Norm *= 3./2./sqr(fl[i].Mass());
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
      if (fl[i].IsVector() && !ATOOLS::IsZero(fl[i].Mass())) {
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

void Polarisation::Reset_Gauge_Vectors(int N,Vec4D* p,Vec4D gauge)
{
#ifndef Explicit_Pols
  if (npol==1) p[N] = Vec4D(Vec3D(gauge).Abs(),Vec3D(gauge));
#endif
}

void Polarisation::Set_Gauge_Vectors(int N,Vec4D* p,Vec4D gauge)
{
#ifndef Explicit_Pols
  if (npol==1) p[N] = Vec4D(Vec3D(gauge).Abs(),Vec3D(gauge));
  if (nmass>0) {
    for (short int i=0;i<N;i++) {
      if (mass_pol[i][0]>=0) {
	Vec4D r1,r2;
        Vec4D pisave = p[i];
	double mass = sqrt(p[i].Abs2());	
	double C = 2.*ran.Get()-1.;
	double S = sqrt(1.-C*C);
	double F = 2.*M_PI*ran.Get();
	r1 = mass/2.*Vec4D(1.,S*::sin(F),S*::cos(F),C);
	r2 = Vec4D(r1[0],(-1.)*Vec3D(r1));
	Vec4D help;
	//r2 boost
	help[0] = (p[i][0]*r2[0]+Vec3D(p[i])*Vec3D(r2))/mass;
	double c1 = (r2[0]+help[0])/(mass+p[i][0]);
	p[mass_pol[i][0]] = Vec4D(help[0],Vec3D(r2)+c1*Vec3D(p[i]));  
	//r1 boost
	help[0] = (p[i][0]*r1[0]+Vec3D(p[i])*Vec3D(r1))/mass;
        c1 = (r1[0]+help[0])/(mass+p[i][0]);
	p[mass_pol[i][1]] = Vec4D(help[0],Vec3D(r1)+c1*Vec3D(p[i]));
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
      if (fl[i].IsVector() && ATOOLS::IsZero(fl[i].Mass()) ) {
        for (short int j=i+1;j<N+1;j++) {
	  if ((fl[j].IsVector() && ATOOLS::IsZero(fl[j].Mass()) ) || 
	      (fl[j]==Flavour(kf_pol))) {
	    norm *= BS->Norm(i,j);
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
    if (fl[i].IsVector() && ATOOLS::IsZero(fl[i].Mass())) {
      if (!( (fl[i+1].IsVector() && ATOOLS::IsZero(fl[i+1].Mass())) ||
	     fl[i+1]==Flavour(kf_pol) ) ) {
	//search next boson or pol
	for (short int j=i+1;j<N+1;j++) {
	  if ( (fl[j].IsVector() && ATOOLS::IsZero(fl[j].Mass())) ||
	       fl[j]==Flavour(kf_pol) ) {
	    n->MPolconvert(i+10+1,j+10);
	    break;
	  }
	}	      
	
      }
    }
  }
#else
  for (short int i=0;i<N;i++) {
    if (fl[i].IsVector() && ATOOLS::IsZero(fl[i].Mass()))  n->MPolconvert(i+masslessskip+1,99);
    if (fl[i].IsVector() && !ATOOLS::IsZero(fl[i].Mass())) n->MPolconvert(i+massiveskip,99);
  }  
#endif
}
  

