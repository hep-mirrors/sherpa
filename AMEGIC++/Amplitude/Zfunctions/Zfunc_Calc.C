#include "AMEGIC++/Amplitude/Zfunctions/Zfunc_Calc.H"
#include "AMEGIC++/Amplitude/Zfunc_Generator.H"
#include "AMEGIC++/Main/Point.H"
#include "AMEGIC++/Amplitude/Zfunc.H"
#include "ATOOLS/Org/Smart_Pointer.C"

using namespace AMEGIC;
using namespace ATOOLS;

#define COMPILE__Getter_Function
#define PARAMETER_TYPE AMEGIC::ZFCalc_Key
#define OBJECT_TYPE AMEGIC::Zfunc_Calc
#include "ATOOLS/Org/Getter_Function.C"

Zfunc_Calc::~Zfunc_Calc() 
{
  // for (size_t i(0);i<lorentzlist.size();++i) lorentzlist[i]->Delete();
}

Zfunc_Calc *Zfunc_Calc::GetCopy() const
{
  Zfunc_Calc *calc(ZFCalc_Getter::GetObject(type,ZFCalc_Key(sgen,BS,NULL)));
  if (calc==NULL) THROW(fatal_error,"Internal error.");
  return calc;
}

Kabbala Zfunc_Calc::Do() 
{ 
  std::cerr<<"Error: Virtual method Zfunc_Calc::Do() called!"<<std::endl; 
  return Kabbala();
}

int Zfunc_Calc::GetScalarNumb() 
{ 
  return 0; 
}

void Zfunc_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			 Point *const p,Point *const pf,Point *&pb,
			 int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
  if(GetScalarNumb()>0){
    int scnt(narg-GetScalarNumb());
    if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
    zfc->SetScalarArgs(zf,scnt,pb->left);
    zfc->SetScalarArgs(zf,scnt,pb->right);
    zfc->SetScalarArgs(zf,scnt,pb->middle);
  }
}
