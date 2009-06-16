#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"
#include "AMEGIC++/Amplitude/Zfunc_Generator.H"
#include "AMEGIC++/Main/Point.H"
#include "AMEGIC++/Amplitude/Zfunc.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;


class AnomalousV4_Calc : public Zfunc_Calc, 
			 public Basic_Zfunc,
			 public Basic_Xfunc,
			 public Basic_Mfunc,
			 public Basic_Vfunc,
                         public Unitarityfunc {
public:
  AnomalousV4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
  ~AnomalousV4_Calc() {}
  ATOOLS::Kabbala Do();
  void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
	       Point *const p,Point *const pf,Point *&pb,
	       int *lfnumb,int *canumb);
};

DEFINE_ZFAGCCALC_GETTER(AnomalousV4_Calc,AnomalousV4Calc_Getter,"AV4","anomalous v4 calculator")

AnomalousV4_Calc::AnomalousV4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Unitarityfunc(_sgen,_BS) 
{ 
  type="AV4";
  ncoupl=10;narg=8;pn=4;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("AGauge4",LF_Key()));
  for (short int i=0;i<4;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[4]->SetParticleArg(0,1,2,3);     
}

Kabbala AnomalousV4_Calc::Do() 
{
  Kabbala factor1 = sgen->GetEnumber(coupl[8]);
  Kabbala factor2 = sgen->GetEnumber(coupl[9]);
  Kabbala uf = U();
 
  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3))) {
    return 2*factor1*Z(0,1)*Z(2,3)-factor2*(Z(0,2)*Z(1,3)+Z(0,3)*Z(1,2));
  }
  if (!IsZero(X(0,0)*M(0)) || 
      !IsZero(X(1,1)*M(1)) || 
      !IsZero(X(2,2)*M(2)) || 
      !IsZero(X(3,3)*M(3))) {
    std::cerr<<"Error in AnomalousV4_Calc::Do(): not cutted massive vertex!"<<std::endl;
    abort();
  }
  return uf*(2*factor1*Z(0,1)*Z(2,3)-factor2*(Z(0,2)*Z(1,3)+Z(0,3)*Z(1,2)));

}

void AnomalousV4_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			       Point *const p,Point *const pf,Point *&pb,
			       int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
}












