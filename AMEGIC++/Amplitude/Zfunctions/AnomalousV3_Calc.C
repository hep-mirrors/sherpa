#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"
#include "ATOOLS/Org/Message.H"
#include "AMEGIC++/Amplitude/Zfunc_Generator.H"
#include "AMEGIC++/Main/Point.H"
#include "AMEGIC++/Amplitude/Zfunc.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

class AnomalousV3_Calc : public Zfunc_Calc, 
			 public Basic_Zfunc,
			 public Basic_Xfunc,
			 public Basic_Mfunc,
			 public Basic_Vfunc,
			 public Basic_Epsilonfunc,
                         public Unitarityfunc {
public:
  AnomalousV3_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
  ~AnomalousV3_Calc() {}
  ATOOLS::Kabbala Do();
  void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
	       Point *const p,Point *const pf,Point *&pb,
	       int *lfnumb,int *canumb);
};

DEFINE_ZFAGCCALC_GETTER(AnomalousV3_Calc,AnomalousV3Calc_Getter,"AV3","anomalous v3 calculator")

AnomalousV3_Calc::AnomalousV3_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Basic_Epsilonfunc(_sgen,_BS),
  Unitarityfunc(_sgen,_BS)
{ 
  type="AV3";
  ncoupl=14;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("AGauge3",LF_Key()));
  for (short int i=0;i<3;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[3]->SetParticleArg(0,1,2);     
}

Kabbala AnomalousV3_Calc::Do() 
{
  if (!IsZero(X(0,0)) || 
      !IsZero(X(1,1)) || 
      !IsZero(X(2,2))) {
    std::cerr<<"Error in AnomalousV3_Calc::Do(): not cutted vertex!"<<std::endl;
    abort();
  }

  Kabbala f1 = sgen->GetEnumber(coupl[6]);
  Kabbala f2 = sgen->GetEnumber(coupl[7]);
  Kabbala f3 = sgen->GetEnumber(coupl[8]);
  Kabbala f4 = sgen->GetEnumber(coupl[9]);
  Kabbala f5 = sgen->GetEnumber(coupl[10]);
  Kabbala f6 = sgen->GetEnumber(coupl[11]);
  Kabbala f7 = sgen->GetEnumber(coupl[12]);
  Kabbala fa = sgen->GetEnumber(coupl[13]);

  Kabbala t1 = Z(2,1)*(X(0,1)-X(0,2))+Z(2,0)*X(1,2)-Z(0,1)*X(2,1);
  Kabbala t2 = Z(0,1)*X(2,0)-Z(2,0)*X(1,0);
  Kabbala t3 = X(2,1)*X(0,2)*X(1,0)-X(0,1)*X(1,2)*X(2,0)
    +Z(0,1)*(X(2,0)*V(2,1)-X(2,1)*V(2,0))
    +Z(2,0)*(X(1,2)*V(1,0)-X(1,0)*V(1,2))
    +Z(1,2)*(X(0,1)*V(0,2)-X(0,2)*V(0,1));
  Kabbala t4 = Z(0,1)*X(2,0)+Z(2,0)*X(1,0);
  
  int sarg[6];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=ps[2].numb;
  sarg[3]=BS->GetPolNumber(arg[0],arg[1],GetPMass(arg[0],arg[1]));
  sarg[4]=BS->GetPolNumber(arg[4],arg[5],GetPMass(arg[4],arg[5]));
  sarg[5]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));

//   std::cout<<sarg[0]<<" "<<sarg[1]<<" "<<sarg[2]<<" "<<sarg[3]<<" "<<sarg[4]<<" "<<sarg[5]<<std::endl;

  Kabbala t5 = Epsilon(sarg[5],sarg[4],sarg[1],sarg[3],1)-
               Epsilon(sarg[5],sarg[4],sarg[2],sarg[3],1);
  Kabbala t6 = Epsilon(sarg[5],sarg[4],sarg[0],sarg[3],1);
  Kabbala t7 = 
     Epsilon(sarg[4],sarg[2],sarg[0],sarg[3],1)*X(2,1)+
     Epsilon(sarg[1],sarg[5],sarg[0],sarg[3],1)*X(1,2)-
     Epsilon(sarg[4],sarg[5],sarg[0],sarg[3],1)*V(1,2)-
     Epsilon(sarg[1],sarg[2],sarg[0],sarg[3],1)*Z(1,2);
//    std::cout<<f5.Value()<<" "<<f6.Value()<<" "<<f7.Value()<<std::endl;
  Kabbala uf = U();

  return fa*(t1+t2)+uf*((f1-fa)*t1+(f2-fa)*t2-f3*t3+f4*t4+t5*f5+t6*f6-t7*f7);
}

void AnomalousV3_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			       Point *const p,Point *const pf,Point *&pb,
			       int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  pb->cpl.resize(8);
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[2];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[3];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[4];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[5];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[6];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[7];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
}


class AnomalousZZZ_Calc : public Zfunc_Calc, 
			  public Basic_Zfunc,
			  public Basic_Xfunc,
			  public Basic_Mfunc,
			  public Basic_Vfunc,
			  public Basic_Epsilonfunc,
			  public Unitarityfunc  {
public:
  AnomalousZZZ_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
  ~AnomalousZZZ_Calc() {}
  ATOOLS::Kabbala Do();
  void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
	       Point *const p,Point *const pf,Point *&pb,
	       int *lfnumb,int *canumb);
};

DEFINE_ZFAGCCALC_GETTER(AnomalousZZZ_Calc,AnomalousZZZCalc_Getter,"AZZZ","anomalous ZZZ calculator")

AnomalousZZZ_Calc::AnomalousZZZ_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Basic_Epsilonfunc(_sgen,_BS),
  Unitarityfunc(_sgen,_BS)
{ 
  type="AZZZ";
  ncoupl=10;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("AZZZ",LF_Key()));
  for (short int i=0;i<3;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[3]->SetParticleArg(0,1,2);     
}

Kabbala AnomalousZZZ_Calc::Do() 
{
  if (!IsZero(X(0,0)) || 
      !IsZero(X(1,1)) || 
      !IsZero(X(2,2))) {
    std::cerr<<"Error in AnomalousZZZ_Calc::Do(): not cutted vertex!"<<std::endl;
    abort();
  }
  Kabbala mz2 = sgen->GetEnumber(sqr(Flavour(kf_Z).Mass()));

  Kabbala f4 = sgen->GetEnumber(coupl[6]);
  Kabbala f5 = sgen->GetEnumber(coupl[7]);

  Kabbala tf4 = (V(0,0)-mz2)*(X(1,0)*Z(0,2)+X(2,0)*Z(0,1))
    + (V(1,1)-mz2)*(X(2,1)*Z(0,1)+X(0,1)*Z(1,2))
    + (V(2,2)-mz2)*(X(0,2)*Z(1,2)+X(1,2)*Z(0,2));

  int sarg[6];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=ps[2].numb;
  sarg[3]=BS->GetPolNumber(arg[0],arg[1],GetPMass(arg[0],arg[1]));
  sarg[4]=BS->GetPolNumber(arg[4],arg[5],GetPMass(arg[4],arg[5]));
  sarg[5]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));

  Kabbala tf5 = (V(0,0)-V(2,2))*( Epsilon(sarg[3],sarg[4],sarg[5],sarg[1],1)
				 -Epsilon(sarg[3],sarg[4],sarg[5],sarg[2],1))
    + (V(1,1)-V(2,2))*( Epsilon(sarg[3],sarg[4],sarg[5],sarg[2],1)
		       -Epsilon(sarg[3],sarg[4],sarg[5],sarg[0],1));
  Kabbala uf = U();

  return uf*(f4*tf4+f5*tf5);

}

void AnomalousZZZ_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
				Point *const p,Point *const pf,Point *&pb,
				int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[2];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[3];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
}



class AnomalousZZG_Calc : public Zfunc_Calc, 
			  public Basic_Zfunc,
			  public Basic_Xfunc,
			  public Basic_Mfunc,
			  public Basic_Vfunc,
			  public Basic_Epsilonfunc,
			  public Unitarityfunc  {
public:
  AnomalousZZG_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
  ~AnomalousZZG_Calc() {}
  ATOOLS::Kabbala Do();
  void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
	       Point *const p,Point *const pf,Point *&pb,
	       int *lfnumb,int *canumb);
};

DEFINE_ZFAGCCALC_GETTER(AnomalousZZG_Calc,AnomalousZZGCalc_Getter,"AZZG","anomalous ZZG calculator")

AnomalousZZG_Calc::AnomalousZZG_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Basic_Epsilonfunc(_sgen,_BS),
  Unitarityfunc(_sgen,_BS)
{ 
  type="AZZG";
  ncoupl=14;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("AZZG",LF_Key()));
  for (short int i=0;i<3;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[3]->SetParticleArg(0,1,2);     
}

Kabbala AnomalousZZG_Calc::Do() 
{
  if (!IsZero(X(0,0)) || 
      !IsZero(X(1,1)) || 
      !IsZero(X(2,2))) {
    std::cerr<<"Error in AnomalousZZG_Calc::Do(): not cutted vertex!"<<std::endl;
    abort();
  }
  Kabbala mz2 = sgen->GetEnumber(sqr(Flavour(kf_Z).Mass()));

  Kabbala f4 = sgen->GetEnumber(coupl[6]);
  Kabbala f5 = sgen->GetEnumber(coupl[7]);
  Kabbala h1 = sgen->GetEnumber(coupl[8]);
  Kabbala h2 = sgen->GetEnumber(coupl[9]);
  Kabbala h3 = sgen->GetEnumber(coupl[10]);
  Kabbala h4 = sgen->GetEnumber(coupl[11]);

  int sarg[6];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=ps[2].numb;
  sarg[3]=BS->GetPolNumber(arg[0],arg[1],GetPMass(arg[0],arg[1]));
  sarg[4]=BS->GetPolNumber(arg[4],arg[5],GetPMass(arg[4],arg[5]));
  sarg[5]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));

  Kabbala tf4 = V(0,0)*(X(1,0)*Z(0,2)+X(2,0)*Z(0,1));

  Kabbala tf5 = V(0,0)*( Epsilon(sarg[3],sarg[4],sarg[5],sarg[1],1)
		        -Epsilon(sarg[3],sarg[4],sarg[5],sarg[2],1));

  Kabbala th1 = (V(2,2)-V(1,1))*(X(1,0)*Z(0,2)-X(2,0)*Z(0,1));
  
  Kabbala th2 = (V(2,2)-mz2)*X(1,2)*(X(2,0)*X(0,2)-V(0,2)*Z(0,2))
    + (V(1,1)-mz2)*X(2,1)*(X(1,0)*X(0,1)-V(0,1)*Z(0,1));

  Kabbala th3 = (V(1,1)-V(2,2))*Epsilon(sarg[3],sarg[4],sarg[5],sarg[0],1);
  Kabbala th4 = (V(2,2)-mz2)*X(1,2)*Epsilon(sarg[5],sarg[3],sarg[2],sarg[0],1)
    + (V(1,1)-mz2)*X(2,1)*Epsilon(sarg[4],sarg[3],sarg[1],sarg[0],1);
  Kabbala uf = U();

  return uf*(f4*tf4+f5*tf5+h1*th1+h2*th2+h3*th3+h4*th4);

}

void AnomalousZZG_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
				Point *const p,Point *const pf,Point *&pb,
				int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  pb->cpl.resize(8);
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[2];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[3];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[4];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[5];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[6];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[7];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
}


class AnomalousZGG_Calc : public Zfunc_Calc, 
			  public Basic_Zfunc,
			  public Basic_Xfunc,
			  public Basic_Mfunc,
			  public Basic_Vfunc,
			  public Basic_Epsilonfunc,
			  public Unitarityfunc  {
public:
  AnomalousZGG_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
  ~AnomalousZGG_Calc() {}
  ATOOLS::Kabbala Do();
  void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
	       Point *const p,Point *const pf,Point *&pb,
	       int *lfnumb,int *canumb);
};

DEFINE_ZFAGCCALC_GETTER(AnomalousZGG_Calc,AnomalousZGGCalc_Getter,"AZGG","anomalous ZGG calculator")

AnomalousZGG_Calc::AnomalousZGG_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Basic_Epsilonfunc(_sgen,_BS),
  Unitarityfunc(_sgen,_BS)
{ 
  type="AZGG";
  ncoupl=10;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("AZGG",LF_Key()));
  for (short int i=0;i<3;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[3]->SetParticleArg(0,1,2);     
}

Kabbala AnomalousZGG_Calc::Do() 
{
  if (!IsZero(X(0,0)) || 
      !IsZero(X(1,1)) || 
      !IsZero(X(2,2))) {
    std::cerr<<"Error in AnomalousZGG_Calc::Do(): not cutted vertex!"<<std::endl;
    abort();
  }
  Kabbala h1 = sgen->GetEnumber(coupl[6]);
  Kabbala h2 = sgen->GetEnumber(coupl[7]);
  Kabbala h3 = sgen->GetEnumber(coupl[8]);
  Kabbala h4 = sgen->GetEnumber(coupl[9]);

  int sarg[6];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=ps[2].numb;
  sarg[3]=BS->GetPolNumber(arg[0],arg[1],GetPMass(arg[0],arg[1]));
  sarg[4]=BS->GetPolNumber(arg[4],arg[5],GetPMass(arg[4],arg[5]));
  sarg[5]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));

  Kabbala th1 = V(0,0)*(X(1,2)*Z(0,2)-X(0,2)*Z(1,2))
    + V(2,2)*(X(1,0)*Z(0,2)-X(2,0)*Z(0,1));
  
  Kabbala th2 = (V(0,0)*X(1,0)+V(2,2)*X(1,2))*(X(2,0)*X(0,2)-V(0,2)*Z(0,2));

  Kabbala th3 = V(0,0)*Epsilon(sarg[3],sarg[4],sarg[5],sarg[2],1)
    + V(2,2)*Epsilon(sarg[5],sarg[4],sarg[3],sarg[0],1);

  Kabbala th4 = (V(0,0)*X(1,0)+V(2,2)*X(1,2))*Epsilon(sarg[3],sarg[5],sarg[0],sarg[2],1);
  Kabbala uf = U();

  return uf*(h1*th1+h2*th2+h3*th3+h4*th4);

}

void AnomalousZGG_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
				Point *const p,Point *const pf,Point *&pb,
				int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[2];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[3];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
}
