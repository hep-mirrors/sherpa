#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/String/String_Generator.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Kabbala Basic_Zfunc::Z(const int z1,const int z2)
{
  int sarg[8];

  for (short int i=0;i<4;++i) sarg[i]   = arg[4*z1+i];
  for (short int i=0;i<4;++i) sarg[i+4] = arg[4*z2+i];

  for (short int i=0;i<4;++i) {
    if (sarg[2*i]==99) {
      if (i<2) return X(z2,z1,1-i);
      else     return X(z1,z2,3-i);
    } 
  }

  Complex scoupl[4];

  scoupl[0] = coupl[2*z1];  scoupl[1] = coupl[2*z1+1];
  scoupl[2] = coupl[2*z2];  scoupl[3] = coupl[2*z2+1];

  return sgen->GetZnumber(sarg,scoupl,
			   Zcalc(arg[4*z1],arg[4*z1+1],arg[4*z1+2],arg[4*z1+3],
				 arg[4*z2],arg[4*z2+1],arg[4*z2+2],arg[4*z2+3],
				 coupl[2*z1],coupl[2*z1+1],coupl[2*z2],coupl[2*z2+1]));
}


Complex Basic_Zfunc::Zcalc(const int t1,const int sign1,const int t2,const int sign2,
			   const int t3,const int sign3,const int t4,const int sign4,
			   const Complex& cR1,const Complex& cL1,
			   const Complex& cR2,const Complex& cL2)
{
  Complex Zhelp;// = Complex(0.,0.);

  int sum = sign1+sign2+sign3+sign4;

  if (sum==4) {
    Zhelp  = BS->S0(t3,t1)*BS->S1(t4,t2)*cR1*cR2;
    Zhelp -= BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cR2*cL1;
    Zhelp -= BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cR1*cL2;
    return -2.*Zhelp;
  }
  if (sum==-4) {
    Zhelp  = BS->S1(t3,t1)*BS->S0(t4,t2)*cL1*cL2;
    Zhelp -= BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cL2*cR1;
    Zhelp -= BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cL1*cR2;
    return -2.*Zhelp;
  }
  if (sum==2) {
    if (sign4==-1) {
      if(t3==t4) if(ATOOLS::IsEqual(cL2,cR2)) return 0.;
      Zhelp = BS->Eta(t2)*cR1*( BS->S0(t4,t1)*BS->Mu(t3)*cL2 -
				BS->S0(t3,t1)*BS->Mu(t4)*cR2 );
      return -2.*Zhelp;
    }
    if (sign3==-1) {
      if(t3==t4) if(ATOOLS::IsEqual(cL2,cR2)) return 0.;
      Zhelp = BS->Eta(t1)*cR1*( BS->S1(t2,t3)*BS->Mu(t4)*cL2 -
				BS->S1(t2,t4)*BS->Mu(t3)*cR2 );
      return -2.*Zhelp;
    }
    if (sign2==-1) {
      if(t1==t2) if(ATOOLS::IsEqual(cL1,cR1)) return 0.;
      Zhelp = BS->Eta(t4)*cR2*( BS->S0(t3,t1)*BS->Mu(t2)*cR1 -
				BS->S0(t3,t2)*BS->Mu(t1)*cL1 );
      return -2.*Zhelp;
    }
    if (sign1==-1) {
      if(t1==t2) if(ATOOLS::IsEqual(cL1,cR1)) return 0.;
      Zhelp = BS->Eta(t3)*cR2*( BS->S1(t2,t4)*BS->Mu(t1)*cR1 -
				BS->S1(t1,t4)*BS->Mu(t2)*cL1 );
      return -2.*Zhelp;
    }
  }
  if (sum==-2) {
    if (sign4==1) {
      if(t3==t4) if(ATOOLS::IsEqual(cL2,cR2)) return 0.;
      Zhelp = BS->Eta(t2)*cL1*( BS->S1(t4,t1)*BS->Mu(t3)*cR2 -
				BS->S1(t3,t1)*BS->Mu(t4)*cL2 );
      return -2.*Zhelp;
    }
    if (sign3==1) {
      if(t3==t4) if(ATOOLS::IsEqual(cL2,cR2)) return 0.;
      Zhelp = BS->Eta(t1)*cL1*( BS->S0(t2,t3)*BS->Mu(t4)*cR2 -
				BS->S0(t2,t4)*BS->Mu(t3)*cL2 );
      return -2.*Zhelp;
    }
    if (sign2==1) {
      if(t1==t2) if(ATOOLS::IsEqual(cL1,cR1)) return 0.;
      Zhelp = BS->Eta(t4)*cL2*( BS->S1(t3,t1)*BS->Mu(t2)*cL1 -
				BS->S1(t3,t2)*BS->Mu(t1)*cR1 );
      return -2.*Zhelp;
    }
    if (sign1==1) {
      if(t1==t2) if(ATOOLS::IsEqual(cL1,cR1)) return 0.;
      Zhelp = BS->Eta(t3)*cL2*( BS->S0(t2,t4)*BS->Mu(t1)*cL1 -
				BS->S0(t1,t4)*BS->Mu(t2)*cR1 );
      return -2.*Zhelp;
    }
  }
  if (sum==0) {
    if ((sign2==-1) && (sign3==-1)) {
      if(t1==t2) if(ATOOLS::IsEqual(cL1,cR1)) return 0.;
      if(t3==t4) if(ATOOLS::IsEqual(cL2,cR2)) return 0.;
      Zhelp = BS->Mu(t1)*BS->Mu(t4)*BS->Eta(t2)*BS->Eta(t3)*cL2*cL1 +
	      BS->Mu(t2)*BS->Mu(t3)*BS->Eta(t1)*BS->Eta(t4)*cR2*cR1 -
	      BS->Mu(t1)*BS->Mu(t3)*BS->Eta(t2)*BS->Eta(t4)*cR2*cL1 -
	      BS->Mu(t2)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t3)*cL2*cR1;
      return -2.*Zhelp;
    }
    if ((sign2==1) && (sign3==1)) {
      if(t1==t2) if(ATOOLS::IsEqual(cL1,cR1)) return 0.;
      if(t3==t4) if(ATOOLS::IsEqual(cL2,cR2)) return 0.;
      Zhelp = BS->Mu(t1)*BS->Mu(t4)*BS->Eta(t2)*BS->Eta(t3)*cR2*cR1 +
  	      BS->Mu(t2)*BS->Mu(t3)*BS->Eta(t1)*BS->Eta(t4)*cL2*cL1 -
       	      BS->Mu(t1)*BS->Mu(t3)*BS->Eta(t2)*BS->Eta(t4)*cL2*cR1 -
	      BS->Mu(t2)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t3)*cR2*cL1;
      return -2.*Zhelp;
    }
    if ((sign2==-1) && (sign4==-1)) {
      Zhelp = Complex(0.,0.);
      return Zhelp;
    }
    if ((sign2==1) && (sign4==1)) {
      Zhelp = Complex(0.,0.);
      return Zhelp;
    }
    if ((sign3==-1) && (sign4==-1)) {
      Zhelp  = BS->S0(t1,t4)*BS->S1(t2,t3)*cL2*cR1;
      Zhelp -= BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cL2*cL1;
      Zhelp -= BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cR2*cR1;
      return -2.*Zhelp;
    }
    if ((sign3==1) && (sign4==1)) {
      Zhelp  = BS->S1(t1,t4)*BS->S0(t2,t3)*cR2*cL1;
      Zhelp -= BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cR2*cR1;
      Zhelp -= BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cL2*cL1;
      return -2.*Zhelp;
    }
  }
  return 0.;
}

int Basic_Zfunc::Zmassless(const int t1,const int sign1,const int t2,const int sign2,
			       const int t3,const int sign3,const int t4,const int sign4,
			       const Complex& cR1,const Complex& cL1,
			       const Complex& cR2,const Complex& cL2)
{
  int sum = sign1+sign2+sign3+sign4;

  if (sum==4) {
    Complex part1 = BS->S0(t3,t1)*BS->S1(t4,t2)*cR1*cR2;
    Complex part2 = 
       BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cR2*cL1
      +BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cR1*cL2;
    if (ATOOLS::IsZero(part2/(part1-part2))) return 1;
    return 0;
  }
  if (sum==-4) {
    Complex part1 = BS->S1(t3,t1)*BS->S0(t4,t2)*cL1*cL2;
    Complex part2 = 
       BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cL2*cR1
      +BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cL1*cR2;
    if (ATOOLS::IsZero(part2/(part1-part2))) return 1;
    return 0;
  }
  if (sum==2) return 0;
  if (sum==-2) return 0;
  if (sum==0) {
    if ((sign2==-1) && (sign3==-1)) return 0;
    if ((sign2==1) && (sign3==1))   return 0;
    if ((sign2==-1) && (sign4==-1)) return 0;
    if ((sign2==1) && (sign4==1)) return 0;
    if ((sign3==-1) && (sign4==-1)) {
      Complex part1 = BS->S0(t1,t4)*BS->S1(t2,t3)*cL2*cR1;
      Complex part2 = 
  	 BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cL2*cL1
	+BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cR2*cR1;
      if (ATOOLS::IsZero(part2/(part1-part2))) return 1;
      return 0;
    }
    if ((sign3==1) && (sign4==1)) {
      Complex part1 = BS->S1(t1,t4)*BS->S0(t2,t3)*cR2*cL1;
      Complex part2 = 
	 BS->Mu(t1)*BS->Mu(t2)*BS->Eta(t3)*BS->Eta(t4)*cR2*cR1
	+BS->Mu(t3)*BS->Mu(t4)*BS->Eta(t1)*BS->Eta(t2)*cL2*cL1;
      if (ATOOLS::IsZero(part2/(part1-part2))) return 1;
      return 0;
    }
  }
  return 0;
}





