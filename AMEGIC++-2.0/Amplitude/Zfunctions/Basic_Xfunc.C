#include "Basic_Func.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace AMATOOLS;

Kabbala Basic_Xfunc::X(const int& t1,const int& sign1,const int& t2,
		       const int& t3,const int& sign3,
		       const Complex& cR,const Complex& cL)
{

  //if(t1==99||t3==99)AORGTOOLS::msg.Error()<<"*****Xold: Alarm "<<cR<<","<<cL<<endl;


  //Marker for -99
  if(t1==99){
    if(t3==t2){
      //cout<<"X: warning, ep("<<t2<<") droped"<<endl;
      return sgen->Get_Enumber(Complex(0.,0.));
    }
    return Vcplx(BS->Get_Pol_Number(t3,sign3,GetPMass(t3,sign3)),t2);
  }
  //Marker for -99
  if(t3==99){
    if(t1==t2){
      //cout<<"X: warning, ep("<<t2<<") droped"<<endl;
      return sgen->Get_Enumber(Complex(0.,0.));
    }
    return Vcplx(BS->Get_Pol_Number(t1,sign1,GetPMass(t1,sign1)),t2);
  }

  int sarg[5];
  sarg[0] = t1;sarg[1] = sign1;
  sarg[2] = t2;
  sarg[3] = t3;sarg[4] = sign3;
  Complex scoupl[2];
  scoupl[0] = cR;scoupl[1] = cL;

  return sgen->Get_Xnumber(sarg,scoupl,Xcalc(t1,sign1,t2,t3,sign3,cR,cL));
}

Kabbala Basic_Xfunc::X(const int &a, const int &b)
{
  int sarg[5];

  for (short int i=0;i<2;i++) sarg[i] = arg[4*a+i].numb;
  sarg[2] = ps[iabs(b)].numb;
  for (short int i=3;i<5;i++) sarg[i] = arg[4*a+i-1].numb;
  
  int sign = Sign(b)*ps[iabs(b)].direction;
  if (ps[iabs(b)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(b)].numb);

  //Marker for -99
  if(sarg[0]==99){
    if(sarg[3]==sarg[2]){
      //cout<<"X: warning, ep("<<sarg[2]<<") droped"<<endl;
      return sgen->Get_Enumber(Complex(0.,0.));
    }
    return (sign>0) ? 
      Vcplx(BS->Get_Pol_Number(sarg[3],sarg[4],GetPMass(sarg[3],sarg[4])),sarg[2]) 
      :
      -Vcplx(BS->Get_Pol_Number(sarg[3],sarg[4],GetPMass(sarg[3],sarg[4])),sarg[2]); 
  }
  //Marker for -99
  if(sarg[3]==99){
    if(sarg[0]==sarg[2]){
      //cout<<"X: warning, ep("<<sarg[2]<<") droped"<<endl;
      return sgen->Get_Enumber(Complex(0.,0.));
    }
    return (sign>0) ? 
      Vcplx(BS->Get_Pol_Number(sarg[0],sarg[1],GetPMass(sarg[0],sarg[1])),sarg[2])
      :
      -Vcplx(BS->Get_Pol_Number(sarg[0],sarg[1],GetPMass(sarg[0],sarg[1])),sarg[2]);
  }

  return (sign>0) ? 
    sgen->Get_Xnumber(sarg,&coupl[2*a],
		      Xcalc(arg[4*a].numb  ,arg[4*a+1].numb,
			    ps[iabs(b)].numb,
			    arg[4*a+2].numb,arg[4*a+3].numb,
			    coupl[2*a],coupl[2*a+1]))
    :
    -sgen->Get_Xnumber(sarg,&coupl[2*a],
		       Xcalc(arg[4*a].numb  ,arg[4*a+1].numb,
			     ps[iabs(b)].numb,
			     arg[4*a+2].numb,arg[4*a+3].numb,
			     coupl[2*a],coupl[2*a+1]));
}

Kabbala Basic_Xfunc::X(const int &a, const int &b, const int &m)
{
  int sarg[5];

  for (short int i=0;i<2;i++) sarg[i] = arg[4*a+i].numb;
  sarg[2] = BS->Get_Pol_Number(arg[4*b+2*m].numb,arg[4*b+2*m+1].numb,
			       GetPMass(arg[4*b+2*m].numb,arg[4*b+2*m+1].numb));
  for (short int i=3;i<5;i++) sarg[i] = arg[4*a+i-1].numb;
  
  //Marker for -99
  if(sarg[0]==99){
    if(sarg[3]==sarg[2]){
      //cout<<"X: warning, ep("<<sarg[2]<<") droped"<<endl;
      return sgen->Get_Enumber(Complex(0.,0.));
    }
    return Vcplx(BS->Get_Pol_Number(sarg[3],sarg[4],GetPMass(sarg[3],sarg[4])),sarg[2]); 
  }
  //Marker for -99
  if(sarg[3]==99){
    if(sarg[0]==sarg[2]){
      //cout<<"X: warning, ep("<<sarg[2]<<") droped"<<endl;
      return sgen->Get_Enumber(Complex(0.,0.));
    }
    return Vcplx(BS->Get_Pol_Number(sarg[0],sarg[1],GetPMass(sarg[0],sarg[1])),sarg[2]);
  }

  return sgen->Get_Xnumber(sarg,&coupl[2*a],
		      Xcalc(arg[4*a].numb  ,arg[4*a+1].numb,
			    sarg[2],
			    arg[4*a+2].numb,arg[4*a+3].numb,
			    coupl[2*a],coupl[2*a+1]));
}


Complex Basic_Xfunc::Xcalc(const int& t1,const int& sign1,const int& t2,
			   const int& t3,const int& sign3,
			   const Complex& cR,const Complex& cL)
{
  
  double sum = sign1+sign3;

  if (sum==2)
    return cL*BS->eta(t2)*BS->eta(t2)*BS->mu(t1)*BS->mu(t3)+
           cR*(BS->eta(t1)*BS->eta(t3)*BS->mu(t2)*BS->mu(t2)+
	       BS->S0(t1,t2)*BS->S1(t2,t3));
  
  if (sum==-2)
    return cR*BS->eta(t2)*BS->eta(t2)*BS->mu(t1)*BS->mu(t3)+
           cL*(BS->eta(t1)*BS->eta(t3)*BS->mu(t2)*BS->mu(t2)+
	       BS->S1(t1,t2)*BS->S0(t2,t3));

  if (sign1==1) 
    return BS->eta(t2)*(cL*BS->mu(t1)*BS->S0(t2,t3)+
			cR*BS->mu(t3)*BS->S0(t1,t2));

  if (sign3==1) 
    return BS->eta(t2)*(cR*BS->mu(t1)*BS->S1(t2,t3)+
			cL*BS->mu(t3)*BS->S1(t1,t2));
  
  return 0.;
}









