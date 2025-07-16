#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;

Kabbala Basic_Epsilonfunc::Epsilon(const int a,const int b,const int c,const int d,int s)
{
  int arg[4];
  arg[0]=a;
  arg[1]=b;
  arg[2]=c;
  arg[3]=d;

  int rd=0;
  for (;rd==0;) {
    rd=1;
    for (size_t i=0;i<3;i++) {
      if (BS->IsComplex(arg[i])<BS->IsComplex(arg[i+1])) {
	rd=arg[i];arg[i]=arg[i+1];arg[i+1]=rd;
	rd=0;
	s=-s;
      }
    }
  }
  rd=4;
  for (size_t i=0;i<4;i++) {
    if (!BS->IsComplex(arg[i])) {
      rd=i;
      break;
    }
  }

  Complex res=EpsCalc(arg[0],arg[1],arg[2],arg[3],rd);
  
  return (s>0) ?
    sgen->GetEpsnumber(arg,rd,res)
    :
    -sgen->GetEpsnumber(arg,rd,res);
}

double Basic_Epsilonfunc::EC(const Vec4D* va,const Vec4D* vb,const Vec4D* vc,const Vec4D* vd)
{ 
  return (*va)[0]*( (*vb)[1]*((*vc)[2]*(*vd)[3]-(*vc)[3]*(*vd)[2])
	           +(*vb)[2]*((*vc)[3]*(*vd)[1]-(*vc)[1]*(*vd)[3])
	           +(*vb)[3]*((*vc)[1]*(*vd)[2]-(*vc)[2]*(*vd)[1]))
        +(*va)[1]*( (*vb)[0]*((*vc)[3]*(*vd)[2]-(*vc)[2]*(*vd)[3])
	           +(*vb)[2]*((*vc)[0]*(*vd)[3]-(*vc)[3]*(*vd)[0])
	           +(*vb)[3]*((*vc)[2]*(*vd)[0]-(*vc)[0]*(*vd)[2]))
        +(*va)[2]*( (*vb)[0]*((*vc)[1]*(*vd)[3]-(*vc)[3]*(*vd)[1])
	           +(*vb)[1]*((*vc)[3]*(*vd)[0]-(*vc)[0]*(*vd)[3])
	           +(*vb)[3]*((*vc)[0]*(*vd)[1]-(*vc)[1]*(*vd)[0]))
        +(*va)[3]*( (*vb)[0]*((*vc)[2]*(*vd)[1]-(*vc)[1]*(*vd)[2])
	           +(*vb)[1]*((*vc)[0]*(*vd)[2]-(*vc)[2]*(*vd)[0])
	           +(*vb)[2]*((*vc)[1]*(*vd)[0]-(*vc)[0]*(*vd)[1]));
}
  
Complex Basic_Epsilonfunc::EpsCalc(const int a,const int b,const int c,const int d,const int nc) {
  switch (nc) {
  case 0:
    return EpsCalc<0>(a,b,c,d);
  case 1:
    return EpsCalc<1>(a,b,c,d);
  case 2:
    return EpsCalc<2>(a,b,c,d);
  case 3:
    return EpsCalc<3>(a,b,c,d);
  case 4:
    return EpsCalc<4>(a,b,c,d);
  }
  return Complex(0.,0.);
}

template <>
Complex Basic_Epsilonfunc::EpsCalc<0>(const int a,const int b,const int c,const int d) {
  Vec4D* va = &(BS->Momentum(a));
  Vec4D* vb = &(BS->Momentum(b));
  Vec4D* vc = &(BS->Momentum(c));
  Vec4D* vd = &(BS->Momentum(d));
  return Complex(EC(va,vb,vc,vd),0.);
}

template <>
Complex Basic_Epsilonfunc::EpsCalc<1>(const int a,const int b,const int c,const int d) {
  Vec4D* va  = &(BS->Momentum(a));
  Vec4D* vai = &(BS->MomentumImg(a));
  Vec4D* vb  = &(BS->Momentum(b));
  Vec4D* vc  = &(BS->Momentum(c));
  Vec4D* vd  = &(BS->Momentum(d));
  return Complex(EC(va,vb,vc,vd),EC(vai,vb,vc,vd));
}

template <>
Complex Basic_Epsilonfunc::EpsCalc<2>(const int a,const int b,const int c,const int d) {
  Vec4D* va  = &(BS->Momentum(a));
  Vec4D* vai = &(BS->MomentumImg(a));
  Vec4D* vb  = &(BS->Momentum(b));
  Vec4D* vbi = &(BS->MomentumImg(b));
  Vec4D* vc  = &(BS->Momentum(c));
  Vec4D* vd  = &(BS->Momentum(d));
  return Complex(EC(va,vb,vc,vd)-EC(vai,vbi,vc,vd),EC(vai,vb,vc,vd)+EC(va,vbi,vc,vd));
}

template <>
Complex Basic_Epsilonfunc::EpsCalc<3>(const int a,const int b,const int c,const int d) {
  Vec4D* va  = &(BS->Momentum(a));
  Vec4D* vai = &(BS->MomentumImg(a));
  Vec4D* vb  = &(BS->Momentum(b));
  Vec4D* vbi = &(BS->MomentumImg(b));
  Vec4D* vc  = &(BS->Momentum(c));
  Vec4D* vci = &(BS->MomentumImg(c));
  Vec4D* vd  = &(BS->Momentum(d));
  return Complex(EC(va,vb,vc,vd)-EC(vai,vbi,vc,vd)-EC(vai,vb,vci,vd)-EC(va,vbi,vci,vd),
		 EC(vai,vb,vc,vd)+EC(va,vbi,vc,vd)+EC(va,vb,vci,vd)-EC(vai,vbi,vci,vd));
}

template <>
Complex Basic_Epsilonfunc::EpsCalc<4>(const int a,const int b,const int c,const int d) {
  Vec4D* va  = &(BS->Momentum(a));
  Vec4D* vai = &(BS->MomentumImg(a));
  Vec4D* vb  = &(BS->Momentum(b));
  Vec4D* vbi = &(BS->MomentumImg(b));
  Vec4D* vc  = &(BS->Momentum(c));
  Vec4D* vci = &(BS->MomentumImg(c));
  Vec4D* vd  = &(BS->Momentum(d));
  Vec4D* vdi = &(BS->MomentumImg(d));
  return Complex(EC(va,vb,vc,vd)+EC(vai,vbi,vci,vdi)
		 -EC(vai,vbi,vc,vd)-EC(vai,vb,vci,vd)-EC(va,vbi,vci,vd)
		 -EC(va,vbi,vci,vd)-EC(va,vbi,vc,vdi)-EC(va,vb,vci,vdi),
		 EC(vai,vb,vc,vd)+EC(va,vbi,vc,vd)+EC(va,vb,vci,vd)+EC(va,vb,vc,vdi)
		 -EC(va,vbi,vci,vdi)-EC(vai,vb,vci,vdi)
		 -EC(vai,vbi,vc,vdi)-EC(vai,vbi,vci,vd));
}




