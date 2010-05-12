#include "ATOOLS/Org/Message.H"
#include "AMEGIC++/Amplitude/Single_Amplitude.H"
#include "AMEGIC++/Amplitude/Prop_Generator.H"
#include "AMEGIC++/Amplitude/Color_Generator.H"
#include "ATOOLS/Math/Kabbala.H"
#include "AMEGIC++/Amplitude/Zfunc_Generator.H"
#include "AMEGIC++/Main/Process_Tags.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

#define Cut_Fermion_Prop

Single_Amplitude::Single_Amplitude(Point* _p,int _topnum, int _permnum,int* _b,int dep,int _no,
				   Topology* top,
				   Basic_Sfuncs* _BS,
				   Flavour* _fl,
				   String_Handler* _shand) 
  : Single_Amplitude_Base(_b,_no,_BS,_fl,_shand)
{

  static int first = 1;
  if (first) first = 0; 
  topnum  = _topnum;
  permnum = _permnum;
  
  icoul = 0;
  on = 1;
  Pointlist = new Point[dep];
  Next = 0;
  int ll = 0;
  top->Copy(_p,Pointlist,ll);

  CFlist  = NULL;
  CCFlist = NULL;
  spind   = NULL;
  SetStringOn(); 

  static int ampltotalnumber = 0;
  ampltotalnumber++;
  amplnumber = ampltotalnumber;
}


Single_Amplitude::Single_Amplitude(Point* _p,int* _b,int dep,int _no,
				   Topology* top,
				   Basic_Sfuncs* _BS,
				   Flavour* _fl,
				   String_Handler* _shand) 
  : Single_Amplitude_Base(_b,_no,_BS,_fl,_shand)
{
  topnum  = 0;
  permnum = 0;

  static int first = 1;

  if (first) first = 0;
  icoul = 0;
  on = 1;
  Pointlist = new Point[dep];
  Next = 0;
  int ll = 0;
  top->Copy(_p,Pointlist,ll);

  CFlist  = NULL;
  CCFlist = NULL;
  spind   = NULL;
  SetStringOn(); 

  static int ampltotalnumber = 0;
  ampltotalnumber++;
  amplnumber = ampltotalnumber;
}

Single_Amplitude::Single_Amplitude(int* _b,int _no,Process_Tags* pinfo,Single_Amplitude** sglist,
				   Basic_Sfuncs* _BS,ATOOLS::Flavour* _fl,String_Handler* _shand) : 
  Single_Amplitude_Base(_b,_no,_BS,_fl,_shand)
{
  topnum  = 0;
  permnum = 0;

  static int first = 1;
  if (first) first = 0;
  icoul = 0;
  on = 1;
  int nin = 1;
  if (_b[1]==-1) nin=2;
  int ndecays = pinfo->Ndecays();
  Pointlist = new Point[2*N-3+ndecays];
  Point** pl = new Point*[ndecays+1];
  for (int i=0; i<ndecays+1;i++) pl[i]=sglist[i]->GetPointlist();
  pinfo->MergePointList(pl,Pointlist,nin);
  Pointlist->ResetProps();

  Next = 0;

  CFlist  = NULL;
  CCFlist = NULL;
  spind   = NULL;
  SetStringOn(); 

  static int ampltotalnumber = 0;
  ampltotalnumber++;
  amplnumber = ampltotalnumber;
}

Single_Amplitude::~Single_Amplitude() 
{
  delete[] Pointlist;
  //Zlist,Clist,Plist

  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) delete (*pit);
  
  SpinorDirection* sd;
  SpinorDirection* sd2;
  sd = spind;
  while (sd) {
    sd2 = sd;
    sd = sd->Next;
    delete sd2;
  }
  
  if (CFlist)  delete CFlist;
  if (CCFlist) delete CCFlist;

}

Point* Single_Amplitude::GetPointlist() {
  return Pointlist;
} 

void Single_Amplitude::AddSpinorDirection(const int& from,const int& to)
{
  SpinorDirection* sd = new SpinorDirection;
  sd->from = from;
  sd->to   = to;
  sd->Next = 0;

  if (spind) {
    SpinorDirection* help;
    help = spind;
    while (help->Next) help = help->Next; 
    help->Next = sd;  
  }
  else spind = sd;
}



void Single_Amplitude::PrintGraph() 
{
  if (!msg_LevelIsTracking()) return;
  
  msg_Out()<<"--------"<<amplnumber+1<<". Amplitude----------"<<endl;

  Single_Amplitude_Base::PrintGraph();

  Color_Function* c;
  c = CFlist;
  msg_Out()<<"Color-matrix: ";
  while(c) {
    switch (c->Type()) {
      case  0: {
	msg_Out()<<"T("<<c->ParticleArg(0)<<" "<<c->ParticleArg(1)
			    <<" "<<c->ParticleArg(2)<<") ";
	break;
      }
      case  1: {
	msg_Out()<<"F("<<c->ParticleArg(0)<<" "<<c->ParticleArg(1)
			    <<" "<<c->ParticleArg(2)<<") ";
	break;
      }
      case 10: {
	msg_Out()<<"TP("<<c->ParticleArg(0)<<" "<<c->ParticleArg(1)
			    <<" "<<c->ParticleArg(2)<<") ";
	break;
      }
      default : break;
    }
    c = c->Next();     
  }
  msg_Out()<<endl<<"Color-string: "<<Colstring<<endl<<endl<<"Spinflow:"<<endl;
  SpinorDirection* sd = spind;
  while(sd) {
    msg_Out()<<sd->from<<" -> "<<sd->to<<endl;
    sd = sd->Next;     
  }
  msg_Out()<<"Overall sign "<<sign<<endl;
}

void Single_Amplitude::Zprojecting(Flavour* fl,int ngraph,bool gc)
{
  CFlist  = NULL;
  CCFlist = NULL;
  
  if (gc){
    Color_Generator cgen;
    int dummy = 0;
    cgen.CFConvert(N,dummy,Pointlist);  
    cgen.CFKill();
    cgen.CFBuildString(N);  
    CFlist  = cgen.Get_CF();
    CCFlist = cgen.Get_CCF();
    CFColstring = cgen.CF2String(CFlist);
    CFColstringC = cgen.CF2String(CCFlist);
  }        

  Zfunc_Generator zgen(BS);
  zgen.BuildZlist(shand->Get_Generator(),BS,ngraph);
  zgen.LorentzConvert(Pointlist);
  zgen.MarkCut(Pointlist,0);
  zgen.Convert(Pointlist);
  zgen.SetDirection(N,spind);
  zgen.Get(*zlist);
  
  Prop_Generator pgen;
  pgen.Convert(Pointlist);
  pgen.Fill();
  pgen.Kill(*zlist);
  pgen.Get(plist);
}

void Single_Amplitude::FillCoupling(String_Handler* _shand) 
{
  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
    Zfunc* z = (*zit);
    for (short int i=0;i<z->m_ncoupl;i++) {
      (_shand->Get_Generator())->GetCnumber(z->p_couplings[i]);
    }
  }
}

void Single_Amplitude::MPolconvert(int alt,int neu)
{
  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
    Zfunc* z = (*zit);
    for (int i=0;i<z->m_narg;i++) {
      if (z->p_arguments[i]==alt) z->p_arguments[i]=neu;
    } 
  }  
}

void Single_Amplitude::Prop_Replace(Flavour falt,int alt,int neu1,int neu2)
{
  Pfunc* Ph = new Pfunc;

  Ph->on = 0;
  Ph->fl = falt;
  Ph->arg = new int[3];
  Ph->argnum = 3;
  Ph->arg[0] = alt;
  Ph->arg[1] = neu1;
  Ph->arg[2] = neu2;

  plist.push_back(Ph);
}

void Single_Amplitude::SetOrderQED() {
  int oqed =0;
  m_oqed = GetPointlist()->FindQEDOrder(oqed);
}
void Single_Amplitude::SetOrderQCD() {
  int oqcd =0;
  m_oqcd = GetPointlist()->FindQCDOrder(oqcd);
}

int Single_Amplitude::GetOrderQED() { return m_oqed; }

int Single_Amplitude::GetOrderQCD() { return m_oqcd; }






