#include "Basic_Func.H"
#include "String_Generator.H"
#include "Basic_Sfuncs.H"

using namespace AMEGIC;
using namespace std;

Kabbala Basic_Zfunc::Z(const int& z1,const int& z2)
{
  /*if (AMATOOLS::IsZero(coupl[2*z2]) && AMATOOLS::IsZero(coupl[2*z2+1])) {
    cerr<<"Error in Basic_Zfunc::Z()!!!"<<endl;
    abort();
  }
  
  if (AMATOOLS::IsZero(coupl[2*z1]) && AMATOOLS::IsZero(coupl[2*z1+1])) {
    cerr<<"Error in Basic_Zfunc::Z()!!!"<<endl;
    abort();
    }*/

  //for testing purpose -> still to be improved
  int sarg[8];

  for (short int i=0;i<4;i++) sarg[i]   = arg[4*z1+i].numb;
  for (short int i=0;i<4;i++) sarg[i+4] = arg[4*z2+i].numb;

  for (short int i=0;i<7;i+=2) {
    //Marker for -99
    if (sarg[i]==99) {
      int partner;
      switch (i) {
      case 0:partner = 1;break;
      case 2:partner = 0;break;
      case 4:partner = 3;break;
      case 6:partner = 2;break;
      }

      if (partner<2) return X(z2,z1,partner);
      if (partner>1) return X(z1,z2,partner-2);

      //int polnumb = BS->Get_Pol_Number(sarg[partner*2],sarg[partner*2+1],GetPMass(sarg[partner*2],sarg[partner*2+1]));
      //cout<<"Plnumber: "<<polnumb<<endl;
      //if (partner*2<4) return X(sarg[4],sarg[5],polnumb,sarg[6],sarg[7],coupl[2*z2],coupl[2*z2+1]);
      //if (partner*2>2) return X(sarg[0],sarg[1],polnumb,sarg[2],sarg[3],coupl[2*z1],coupl[2*z1+1]);
    }
  }

  Complex scoupl[4];

  scoupl[0] = coupl[2*z1];  scoupl[1] = coupl[2*z1+1];
  scoupl[2] = coupl[2*z2];  scoupl[3] = coupl[2*z2+1];

  //sarg --> &arg[4*z]

  return sgen->Get_Znumber(sarg,scoupl,
			   Zcalc(arg[4*z1].numb,arg[4*z1+1].numb,arg[4*z1+2].numb,arg[4*z1+3].numb,
				 arg[4*z2].numb,arg[4*z2+1].numb,arg[4*z2+2].numb,arg[4*z2+3].numb,
				 coupl[2*z1],coupl[2*z1+1],coupl[2*z2],coupl[2*z2+1]));
}


Complex Basic_Zfunc::Zcalc(const int& t1,const int& sign1,const int& t2,const int& sign2,
			   const int& t3,const int& sign3,const int& t4,const int& sign4,
			   const Complex& cR1,const Complex& cL1,
			   const Complex& cR2,const Complex& cL2)
{
  Complex Zhelp = Complex(0.,0.);

  int sum = sign1+sign2+sign3+sign4;

  if (sum==4) {
    Zhelp  = BS->S0(t3,t1)*BS->S1(t4,t2)*cR1*cR2;
    Zhelp -= BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cR2*cL1;
    Zhelp -= BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cR1*cL2;
    return -2.*Zhelp;
  }
  if (sum==-4) {
    Zhelp  = BS->S1(t3,t1)*BS->S0(t4,t2)*cL1*cL2;
    Zhelp -= BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cL2*cR1;
    Zhelp -= BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cL1*cR2;
    return -2.*Zhelp;
  }
  if (sum==2) {
    if (sign4==-1) {
      Zhelp = BS->eta(t2)*cR1*( BS->S0(t4,t1)*BS->mu(t3)*cL2 -
				BS->S0(t3,t1)*BS->mu(t4)*cR2 );
      return -2.*Zhelp;
    }
    if (sign3==-1) {
      Zhelp = BS->eta(t1)*cR1*( BS->S1(t2,t3)*BS->mu(t4)*cL2 -
				BS->S1(t2,t4)*BS->mu(t3)*cR2 );
      return -2.*Zhelp;
    }
    if (sign2==-1) {
      Zhelp = BS->eta(t4)*cR2*( BS->S0(t3,t1)*BS->mu(t2)*cR1 -
				BS->S0(t3,t2)*BS->mu(t1)*cL1 );
      return -2.*Zhelp;
    }
    if (sign1==-1) {
      Zhelp = BS->eta(t3)*cR2*( BS->S1(t2,t4)*BS->mu(t1)*cR1 -
				BS->S1(t1,t4)*BS->mu(t2)*cL1 );
      return -2.*Zhelp;
    }
  }
  if (sum==-2) {
    if (sign4==1) {
      Zhelp = BS->eta(t2)*cL1*( BS->S1(t4,t1)*BS->mu(t3)*cR2 -
				BS->S1(t3,t1)*BS->mu(t4)*cL2 );
      return -2.*Zhelp;
    }
    if (sign3==1) {
      Zhelp = BS->eta(t1)*cL1*( BS->S0(t2,t3)*BS->mu(t4)*cR2 -
				BS->S0(t2,t4)*BS->mu(t3)*cL2 );
      return -2.*Zhelp;
    }
    if (sign2==1) {
      Zhelp = BS->eta(t4)*cL2*( BS->S1(t3,t1)*BS->mu(t2)*cL1 -
				BS->S1(t3,t2)*BS->mu(t1)*cR1 );
      return -2.*Zhelp;
    }
    if (sign1==1) {
      Zhelp = BS->eta(t3)*cL2*( BS->S0(t2,t4)*BS->mu(t1)*cL1 -
				BS->S0(t1,t4)*BS->mu(t2)*cR1 );
      return -2.*Zhelp;
    }
  }
  if (sum==0) {
    if ((sign2==-1) && (sign3==-1)) {
      Zhelp = BS->mu(t1)*BS->mu(t4)*BS->eta(t2)*BS->eta(t3)*cL2*cL1 +
	BS->mu(t2)*BS->mu(t3)*BS->eta(t1)*BS->eta(t4)*cR2*cR1 -
	BS->mu(t1)*BS->mu(t3)*BS->eta(t2)*BS->eta(t4)*cR2*cL1 -
	BS->mu(t2)*BS->mu(t4)*BS->eta(t1)*BS->eta(t3)*cL2*cR1;
      return -2.*Zhelp;
    }
    if ((sign2==1) && (sign3==1)) {
      Zhelp = BS->mu(t1)*BS->mu(t4)*BS->eta(t2)*BS->eta(t3)*cR2*cR1 +
	BS->mu(t2)*BS->mu(t3)*BS->eta(t1)*BS->eta(t4)*cL2*cL1 -
	BS->mu(t1)*BS->mu(t3)*BS->eta(t2)*BS->eta(t4)*cL2*cR1 -
	BS->mu(t2)*BS->mu(t4)*BS->eta(t1)*BS->eta(t3)*cR2*cL1;
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
      Zhelp -= BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cL2*cL1;
      Zhelp -= BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cR2*cR1;
      return -2.*Zhelp;
    }
    if ((sign3==1) && (sign4==1)) {
      Zhelp  = BS->S1(t1,t4)*BS->S0(t2,t3)*cR2*cL1;
      Zhelp -= BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cR2*cR1;
      Zhelp -= BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cL2*cL1;
      return -2.*Zhelp;
    }
  }
  return 0.;
}

int Basic_Zfunc::Zmassless(const int& t1,const int& sign1,const int& t2,const int& sign2,
			       const int& t3,const int& sign3,const int& t4,const int& sign4,
			       const Complex& cR1,const Complex& cL1,
			       const Complex& cR2,const Complex& cL2)
{
  Complex Zhelp = Complex(0.,0.);

  int sum = sign1+sign2+sign3+sign4;

  if (sum==4) {
    Complex part1 = BS->S0(t3,t1)*BS->S1(t4,t2)*cR1*cR2;
    Complex part2 = 
       BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cR2*cL1
      +BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cR1*cL2;
    if (AMATOOLS::IsZero(part2/(part1-part2))) return 1;
    return 0;
  }
  if (sum==-4) {
    Complex part1 = BS->S1(t3,t1)*BS->S0(t4,t2)*cL1*cL2;
    Complex part2 = 
       BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cL2*cR1
      +BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cL1*cR2;
    if (AMATOOLS::IsZero(part2/(part1-part2))) return 1;
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
  	 BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cL2*cL1
	+BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cR2*cR1;
      if (AMATOOLS::IsZero(part2/(part1-part2))) return 1;
      return 0;
    }
    if ((sign3==1) && (sign4==1)) {
      Complex part1 = BS->S1(t1,t4)*BS->S0(t2,t3)*cR2*cL1;
      Complex part2 = 
	 BS->mu(t1)*BS->mu(t2)*BS->eta(t3)*BS->eta(t4)*cR2*cR1
	+BS->mu(t3)*BS->mu(t4)*BS->eta(t1)*BS->eta(t2)*cL2*cL1;
      if (AMATOOLS::IsZero(part2/(part1-part2))) return 1;
      return 0;
    }
  }
  return 0;
}





