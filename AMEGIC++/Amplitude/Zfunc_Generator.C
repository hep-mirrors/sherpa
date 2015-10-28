#include "AMEGIC++/Amplitude/Zfunc_Generator.H"
#include "AMEGIC++/String/String_Generator.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/Amplitude/Zfunctions/Zfunc_Calc.H"
#include "AMEGIC++/Amplitude/Zfunc.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

#define lorentz_type

extern int iabs(int&);

Zfunc_Generator::~Zfunc_Generator() 
{
}

ZF_Vector Zfunc_Generator::zcalc;

void Zfunc_Generator::BuildZlist(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS, int ngraph)
{
  if (ngraph!=1) return;
  zcalc.clear();
  ZFCalc_Key key(_sgen,_BS,MODEL::s_model);
  ZFCalc_Getter::Getter_List zfclist(ZFCalc_Getter::GetGetters());
  for (ZFCalc_Getter::Getter_List::const_iterator 
	 git(zfclist.begin());git!=zfclist.end();++git) {
    Zfunc_Calc *nc((*git)->GetObject(key));
    if (nc!=NULL) zcalc.push_back(nc);
  }
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n   Implemented calculators:\n\n";
    ZFCalc_Getter::PrintGetterInfo(msg_Out(),15);
    msg_Out()<<"\n   Implemented Lorentz functions:\n\n";
    LF_Getter::PrintGetterInfo(msg_Out(),15);
    msg_Out()<<"\n}\n";
  }
}

void Zfunc_Generator::LorentzConvert(Point* p)
{
  if (p==0) return;
  Lorentz_Function* l = p->Lorentz;

  int partarg[4]={-1,-1,-1,-1};
  for (short int i=0;i<l->NofIndex();i++) {
    switch (l->ParticleArg(i)) {
    case 0: partarg[i] = p->number;break;
    case 1: partarg[i] = p->left->number;break;
    case 2: partarg[i] = p->right->number;break;
    case 3: partarg[i] = p->middle->number;break;
    }
  }
  l->SetParticleArg(partarg[0],partarg[1],partarg[2],partarg[3]);
  LorentzConvert(p->right);
  LorentzConvert(p->left);  
  LorentzConvert(p->middle);
}

void Zfunc_Generator::MarkCut(Point* p,int notcut,bool fromfermion,bool cutvectors)
{
  if (p==0) return; 

  if (cutvecprop && p->fl.IsVector() && p->number>99){
    p->m = 1;
    notcut++;
    if(fromfermion && p->left->fl.IsFermion()){
      p->m=0;
    }
    if (p->fl.Mass()==0.0) p->m=0;

    if (p->Lorentz->CutVectors()||cutvectors) p->m=1;
  }
  else p->m = 0;

  // spin 2 particles must be cutted
  if (p->fl.IsTensor() && p->number>99) p->m = 1;

  // "new gauge test" cut all massless propagators
  if (p->fl.IsVector() && p->number>99  && rpa->gen.CutScheme()==1) {
    if(p->fl.Mass()==0.0) {
      p->m=1;
    }	
  }
  MarkCut(p->right,notcut,p->fl.IsFermion(),p->Lorentz->CutVectors());
  MarkCut(p->left,notcut,p->fl.IsFermion(),p->Lorentz->CutVectors());
  MarkCut(p->middle,notcut,p->fl.IsFermion(),p->Lorentz->CutVectors()); 
}

void Zfunc_Generator::Convert(Point* p)
{
  Zfunc* Zh = 0;
  if ((p->left==0) && (p->right==0)) return;
 
  if ( p->fl.IsFermion() || p->fl.IsScalar() || (p->fl.IsTensor() && p->number<99) ||
       (p->fl.IsVector() && p->number<99) || p->m==1) {
    Zh = new Zfunc;
    Point* pb;
    Point* pf;
    pb = pf = 0;
    
    if (p->fl.IsFermion()) {
      if (p->left->fl.IsBoson()) {
	pb = p->left;
	pf = p->right;
      }
      if (p->right->fl.IsBoson()) {
	pb = p->right;
	pf = p->left;
      }
      if(pb->fl.IsTensor()||pb->fl.IsScalar()){
	if(p->middle)pb =p->middle; 
      }
    }
    else {
      //Incoming Boson
      pb = p;
    }
    if(!LFDetermine_Zfunc(Zh,p,pf,pb)){
      Point* ph1 = pb;
      if (ph1->left->fl.Kfcode()==kf_shgluon) { 
	if (ph1->right->fl.IsScalar() || ph1->right->m==1 || !ph1->right->left) ph1=ph1->left;
	else if (ph1->right->left->fl.IsFermion()) ph1=ph1->left;
      }
      if (ph1->right->fl.Kfcode()==kf_shgluon) { 
	if (ph1->left->fl.IsScalar() || ph1->left->m==1 || !ph1->left->left) ph1=ph1->right;
	else if (ph1->left->left->fl.IsFermion()) ph1=ph1->right;
      }
      Point* ph=ph1->right;
      if (!( ph->fl.IsFermion() || ph->fl.IsScalar() || 
	     (ph->fl.IsVector() && ph->number<99) || ph->m==1 || ph->fl.Kfcode()==kf_shgluon)
	  &&ph->left)
	if(!(ph->left->fl.IsFermion())||ph->middle){
	  ph->m=1;
	  Convert(p); 
	  return;
	}
      ph=ph1->left;
      if (!( ph->fl.IsFermion() || ph->fl.IsScalar() || 
	     (ph->fl.IsVector() && ph->number<99) || ph->m==1 || ph->fl.Kfcode()==kf_shgluon)
	  &&ph->left)
	if(!(ph->left->fl.IsFermion())||ph->middle){
	  ph->m=1;
	  Convert(p); 
	  return;
	}
      if(ph1->middle){
	ph=ph1->middle;
	if (!( ph->fl.IsFermion() || ph->fl.IsScalar() || 
	       (ph->fl.IsVector() && ph->number<99) || ph->m==1 || ph->fl.Kfcode()==kf_shgluon)
	    &&ph->left)
	  if(!(ph->left->fl.IsFermion())||ph->middle){
	    ph->m=1;
	    Convert(p); 
	    return;
	  }
      }
      if(pf!=0 && ((p->Lorentz)->Type()=="FFVT" || 
                   (p->Lorentz)->Type()=="FFVGS"  ) ){
	pb->m=1;
	Convert(p);
	return;
      }
      msg_Error()<<"Zfunc_Generator::Convert(Point* p) : Cutting Error, abort the run."<<endl;
      abort();
    }
  }
  if (Zh) zlist.push_back(Zh);

  Convert(p->right);
  Convert(p->left);
  if (p->middle) Convert(p->middle);
}

void Zfunc_Generator::Lorentz_Sequence(Point* pb,vector<Lorentz_Function*> &lflist)
{ 
  if (pb->left==0 && (pb->fl.IsScalar()||pb->fl.IsTensor())) return;
 
  lflist.push_back(pb->Lorentz->GetCopy());

  //Endpoints are scalar or non-bosons
  int skal,vec;
  IsGaugeV(pb,skal,vec);      
  if (skal+vec>=2) {
    if (pb->left->fl.IsVector() && pb->left->m==0)  Lorentz_Sequence(pb->left,lflist);
    if (pb->right->fl.IsVector() && pb->right->m==0) Lorentz_Sequence(pb->right,lflist);
    if (pb->middle!=0) {
            if (pb->middle->fl.IsVector() && pb->middle->m==0) Lorentz_Sequence(pb->middle,lflist);
      if (pb->middle->m==1 && !pb->middle->fl.IsTensor()) {
	Lorentz_Function *lf(LF_Getter::GetObject("Pol",LF_Key()));
	lf->SetParticleArg(pb->middle->number);
	lflist.push_back(lf);
      }
    }

    if (pb->left->m==1 && !pb->left->fl.IsTensor()) {
      Lorentz_Function *lf(LF_Getter::GetObject("Pol",LF_Key()));
      lf->SetParticleArg(pb->left->number);
      lflist.push_back(lf);
    }
    if (pb->right->m==1 && !pb->right->fl.IsTensor()) {
      Lorentz_Function *lf(LF_Getter::GetObject("Pol",LF_Key()));
      lf->SetParticleArg(pb->right->number);
      lflist.push_back(lf);
    }
  }
}

void Zfunc_Generator::LFPrint(const vector<Lorentz_Function*> &lflist)
{
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<"LorentzList: "<<endl;
  for (size_t i=0;i<lflist.size();i++)
    msg_Out()<<lflist[i]->String(1)<<endl;
  msg_Out()<<endl;
}

std::string Zfunc_Generator::LFEff(const std::string &type)
{ 
  return (type=="Pol") ? "FFV" : type;
}

int Zfunc_Generator::LFDetermine_Zfunc(Zfunc* Zh,Point* p,Point* pf,Point* pb)
{
  Zh->m_type = "";
  vector<Lorentz_Function*> lflist;
  if (pf!=0) lflist.push_back(p->Lorentz->GetCopy());

  if (pf==0 && pb->fl.IsVector()) {
    Lorentz_Function *lf(LF_Getter::GetObject("Pol",LF_Key()));
    lf->SetParticleArg(pb->number);
    lflist.push_back(lf);
  }
  
  if (!( (pb->fl.IsScalar() || pb->m==1) && pf!=0)) Lorentz_Sequence(pb,lflist);
  else {
    if (pb->m==1 && pf!=0 && !pb->fl.IsTensor()) {
      Lorentz_Function *lf(LF_Getter::GetObject("Pol",LF_Key()));
      lf->SetParticleArg(pb->number);
      lflist.push_back(lf);
    }
  }
  //LFPrint(lflist);

  for (size_t i=0;i<zcalc.size();i++) {
    if (lflist.size()==(zcalc[i]->lorentzlist).size()) {
      int hit = 1; 
      vector<std::string> typerem;
      
      for (size_t j=0;j<lflist.size();j++) {
	int hit2 = 1;
	for (size_t k=0;k<typerem.size();k++) {
	  if (typerem[k]==LFEff(lflist[j]->Type())) {
	    hit2 = 0;
	    break;
	  }
	}
	if (hit2) {
	  //counting
	  int type1 = 0;
	  for (size_t k=j;k<lflist.size();k++) {
	    if (LFEff(lflist[j]->Type())==LFEff(lflist[k]->Type())) type1++;
	  }
	  int type2 = 0;
	  for (size_t k=0;k<(zcalc[i]->lorentzlist).size();k++) {
	    if (LFEff(lflist[j]->Type())==LFEff((zcalc[i]->lorentzlist[k])->Type())) type2++;
	  }  
	  if (type1!=type2) {
	    hit = 0;
	    break;
	  }
	  else typerem.push_back(LFEff(lflist[j]->Type()));
	}
      }
      if (hit) {
	Zh->m_type       = zcalc[i]->type;
	Zh->p_calculator = zcalc[i];
	break;
      }
    }
  }
  if (Zh->m_type=="") {
    for (size_t i(0);i<lflist.size();++i) delete lflist[i];
    return 0;
    msg_Error()<<METHOD<<"(): Invalid Lorentz function."<<endl;
    LFPrint(lflist);  
    abort();
  }

  LFFill_Zfunc(Zh,lflist,p,pf,pb);
  for (size_t i(0);i<lflist.size();++i) delete lflist[i];
  return 1;
}


void Zfunc_Generator::CopyOrder(vector<Lorentz_Function*> &lflist,vector<Lorentz_Function*> & lfpointer)
{
  for (size_t i=0;i<lflist.size();i++) lfpointer.push_back(lflist[i]);

  for (size_t i=0;i<lfpointer.size();i++) 
    for (size_t j=i+1;j<lfpointer.size();j++) {
      if (lfpointer[i]->NofIndex()<lfpointer[j]->NofIndex()) {
	Lorentz_Function* help;
	help         = lfpointer[i];
	lfpointer[i] = lfpointer[j];
	lfpointer[j] = help;
      }
    }  
}

int Zfunc_Generator::Compare(int Nargs,
			     const vector<Lorentz_Function*> &lfpointer,
			     int *lfnumb,
			     const vector<Lorentz_Function*> &capointer,
			     int* canumb)
{
  for (short int i=0;i<Nargs;i++) {
    lfnumb[i] = -1;
    canumb[i] = -1;
  }
  
  int numbcount = 0;
  
  for (size_t i=0;i<lfpointer.size();i++) {
    for (int k=0;k<lfpointer[i]->NofIndex();k++) {
      int lfarg = abs(lfpointer[i]->ParticleArg(k));
      int caarg = abs(capointer[i]->ParticleArg(k));
      
      int hit = 1;
      
      for (int j=0;j<numbcount;j++) {
	if (lfnumb[j]==lfarg) {
	  if (canumb[j]==caarg) {
	    hit = 0;
	    break;
	  }
	  else return i;
	}
      }

      if (hit) {
	lfnumb[numbcount] = abs(lfarg);
	canumb[numbcount] = abs(caarg);
	numbcount++;
      }
    }
  }

  return lfpointer.size();
}

void Zfunc_Generator::LFFill_Zfunc(Zfunc* Zh,vector<Lorentz_Function*> &lflist,
				   Point* p,Point* pf,Point* pb)
{
  vector<Lorentz_Function*> lfpointer;
  CopyOrder(lflist,lfpointer);
  vector<Lorentz_Function*> capointer;
  CopyOrder(Zh->p_calculator->lorentzlist,capointer);

  vector<Lorentz_Function*> permpointer;

  for (size_t j=0;j<lfpointer.size();j++) {
    permpointer.push_back(lfpointer[j]);
    permpointer[j]->InitPermutation();
  }

  int* lfnumb = new int[Zh->p_calculator->pn];
  int* canumb = new int[Zh->p_calculator->pn];

  //loop over all permutations
  for (;;) {
    int i = Compare(Zh->p_calculator->pn,lfpointer,lfnumb,capointer,canumb);
    if (i==(int)lfpointer.size()) break;

    //loop over previous permutations
    for (;;) {
      int typecount = 0;
      int typemin   = 1000; 
      int typemax   = 0;

      for (int j=0;j<(int)lfpointer.size();j++) {
	if (LFEff(lfpointer[j]->Type())==LFEff(lfpointer[i]->Type())) {
	  if (typemin>j) typemin = j;
	  if (typemax<j) typemax = j;
	  typecount++;
	}
      }
      
      vector<Lorentz_Function*> copypointer;
      for (size_t j=0;j<lfpointer.size();j++) copypointer.push_back(lfpointer[j]);

      if (typecount>1) {
	int over = 0;
	int* ii = new int[lfpointer.size()]; 
	for (short int j=typemin;j<=typemax;j++) ii[j] = typemin;
	//loop over permutation of one type
	for (;;) {
	  int hit = 1;
	  for (short int j=typemin;j<=typemax;j++) {
	    for (short int k=j+1;k<=typemax;k++) {
	      if (ii[j]==ii[k]) {hit = 0;break;}
	    }
	  }
	  
	  if (hit) {
	    for (short int j=typemin;j<=typemax;j++) copypointer[j] = lfpointer[ii[j]];
	    i = Compare(Zh->p_calculator->pn,copypointer,lfnumb,capointer,canumb);
	    if (i>typemax) {
	      for (short int j=typemin;j<=typemax;j++) lfpointer[j] = copypointer[j];
	      break;
	    }
	  }
	  //next permutation
	  for (short int j=typemax;j>=typemin;j--) {
	    if (ii[j]<typemax) {
	      ii[j]++;            
	      break;
	    }
	    else {
	      ii[j] = typemin;
	      if (j==typemin) over = 1;
	    }
	  }
	  if (over) break;
	}
	delete[] ii;
      }
      else i = Compare(Zh->p_calculator->pn,lfpointer,lfnumb,capointer,canumb);
      if (i==(int)lfpointer.size()) break;
	
      if (i<=typemax) {
	int sn  = 1;
	int max = typemax;
	do {
	  sn = permpointer[max]->NextPermutation();
	  if (sn==0) {
	    permpointer[max]->ResetPermutation();
	    max--;
	    if (max==-1) break;
	  }	  	  
	}
	while (sn==0);

	if (max<typemin) i = typemin-1;
	if (i<0) {
	  msg_Error()<<"ERROR in Zfunc_Generator::LFFill_Zfunc() : abort the run."<<endl;
	  abort();
	}
	//LFPrint(lfpointer);
      }
    }
    if (i==(int)lfpointer.size()) break;
  }

  //Total Sign......
  Zh->m_sign = 1;

  for (size_t j=0;j<permpointer.size();j++) Zh->m_sign *= permpointer[j]->GetSign();

  SetPropDirection(Zh->p_calculator->pn,pb->number,lfpointer,lfnumb,capointer,canumb);


  //Setting the arguments

  Zh->m_narg   = Zh->p_calculator->narg;
  Zh->m_ncoupl = Zh->p_calculator->ncoupl;
  Zh->m_nprop  = Zh->p_calculator->pn;

  Zh->p_arguments   = NULL;
  Zh->p_propagators = NULL;

  if (Zh->m_narg>0)   Zh->p_arguments   = new int[Zh->m_narg];
  if (Zh->m_ncoupl>0) Zh->p_couplings   = new Complex[Zh->m_ncoupl];
  if (Zh->m_nprop>0)  Zh->p_propagators = new Argument[Zh->m_nprop];
  
  //Set all Couplings Zero
  for (short int i=0;i<Zh->m_ncoupl;i++) Zh->p_couplings[i] = Complex(0.,0.);
  for (short int i=0;i<Zh->p_calculator->pn;i++) {
    if (lfnumb[i]==pb->number) {
      Set_In(Zh,canumb[i],p,pf,pb);
      break;
    }
  }

  
  //Special cases
  Zh->p_calculator->SetArgs(this,Zh,p,pf,pb,lfnumb,canumb);

  delete[] lfnumb;
  delete[] canumb;
}

void Zfunc_Generator::SetPropDirection(int Nargs,int incoming,
				       const vector<Lorentz_Function*> &lfpointer,
				       int *lfnumb,
				       const vector<Lorentz_Function*> &capointer,
				       int* canumb)
{
  //Search Incoming
  int start = -1;
  //works only for incoming vectors!!!!
  for (size_t i=0;i<lfpointer.size();i++) {
    if (LFEff(lfpointer[i]->Type())=="FFV") {
      for (short int k=0;k<lfpointer[i]->NofIndex();k++) {
	if (lfpointer[i]->ParticleArg(k)==incoming) {
	  start = i;
	  break;
	}
      }
      if (start!=-1) break;
    }
  }
  
  //pseudo solution for scalars
  if (start!=-1) SearchNextProp(Nargs,lfpointer,lfnumb,capointer,canumb,incoming,start);
}

void Zfunc_Generator::SearchNextProp(int Nargs,
				     const vector<Lorentz_Function*> &lfpointer,
				     int *lfnumb,
				     const vector<Lorentz_Function*> &capointer,
				     int* canumb,
				     int incoming,
				     int position)
{
  int start= -1;
  for (size_t i=0;i<lfpointer.size();i++) {
    if ((int)i!=position) {
      for (short int k=0;k<lfpointer[i]->NofIndex();k++) {
	if (lfpointer[i]->ParticleArg(k)==incoming) {
	  start = i;
	  break;
	}
      }
      if (start!=-1) break;
    }
  }

  if (start==-1) return;

  for (short int k=0;k<lfpointer[position]->NofIndex();k++) {
    if (lfpointer[position]->ParticleArg(k)==incoming) {
      if (capointer[position]->ParticleArg(k)<0) {
	for (short int i=0;i<Nargs;i++) {
	  if (lfnumb[i]==incoming) {
	    canumb[i] = -canumb[i];
	    break;
	  }
	}
      }
      break;
    }
  }
  
  for (short int k=0;k<lfpointer[start]->NofIndex();k++) {
    if (lfpointer[start]->ParticleArg(k)!=incoming) 
      SearchNextProp(Nargs,lfpointer,lfnumb,capointer,canumb,lfpointer[start]->ParticleArg(k),start);
  }
}

void Zfunc_Generator::SetArgs(Zfunc* Zh,int* lfnumb,int* canumb,Point* pb,Point* p,int& icoupl)
{
  if (pb==0||pb->fl.IsTensor()) return;
  
  for (short int i=0;i<Zh->p_calculator->pn;i++) {
    if (lfnumb[i]==pb->number) {
      if (pb->number<99 || pb->fl.IsScalar() || pb->m==1) Set_Out(Zh,canumb[i],pb,p);
      else {
	if  (!pb->left->fl.IsVector() && !pb->right->fl.IsVector()) Set_Out(Zh,canumb[i],pb,p);
	else {
	  Zh->p_couplings[icoupl]                   = pb->cpl[1];icoupl++;
	  Zh->p_propagators[abs(canumb[i])].numb      = pb->number;
	  Zh->p_propagators[abs(canumb[i])].kfcode    = (pb->fl).Kfcode();
	  Zh->p_propagators[abs(canumb[i])].direction = Direction::Outgoing;
	  if (canumb[i]<0) 
	    Zh->p_propagators[abs(canumb[i])].direction = Direction::Incoming;
	  SetArgs(Zh,lfnumb,canumb,pb->left,p,icoupl);
	  SetArgs(Zh,lfnumb,canumb,pb->right,p,icoupl);
	  SetArgs(Zh,lfnumb,canumb,pb->middle,p,icoupl);
	}
      }
      break;
    }
  }
}

void Zfunc_Generator::SetScalarArgs(Zfunc* Zh,int &scnt,Point* pb)
{
  if (pb==0) return;
  
  if (scnt==Zh->m_narg) return;
  if(pb->fl.IsScalar()){
    if  (scnt<Zh->m_narg)  Zh->p_arguments[scnt]=pb->number;
    else{
      Zh->Print();
      msg_Error()<<"ERROR in Zfunc_Generator::SetScalarArgs : "<<std::endl
		 <<"   scnt : "<<scnt<<" Zh->m_narg : "<<Zh->m_narg<<", will abort."<<endl;
      abort();
    }
    scnt++;
    return;
  }
  if (pb->number<99 || pb->m==1) {}
  else {
    if  (!pb->left->fl.IsVector() && !pb->right->fl.IsVector() && !pb->middle) {}
    else {
      SetScalarArgs(Zh,scnt,pb->left);
      SetScalarArgs(Zh,scnt,pb->right);
      SetScalarArgs(Zh,scnt,pb->middle);
    }
  }
}


void Zfunc_Generator::Set_In(Zfunc* Zh,int number, Point* p, Point* pf,Point* pb)
{
  if(p->fl.IsTensor())return;
  int nb=number;
  if (Zh->m_type=="FFVT"||Zh->m_type=="FFVGS")nb--;
  if (Zh->m_nprop>nb && nb>=0) {
    Zh->p_propagators[nb].numb      = pb->number;
    Zh->p_propagators[nb].kfcode    = (pb->fl).Kfcode();
    Zh->p_propagators[nb].direction = Direction::Incoming;
  }
  if (pf!=0) {
    if (pb->m==1) {
      if (Zh->m_nprop>nb && nb>=0)Zh->p_propagators[nb].direction = Direction::Outgoing;
    }
    if (pb->number<99) {
      if (BS->Sign(pb->number)==1) {
	if (Zh->m_nprop>nb && nb>=0)Zh->p_propagators[nb].direction = Direction::Outgoing;
      }
    }

    if ((p->fl).IsAnti()) {
      Zh->p_arguments[number*2+1] = p->number;
      Zh->p_arguments[number*2]   = pf->number;
    }
    else {
      Zh->p_arguments[number*2]   = p->number;
      Zh->p_arguments[number*2+1] = pf->number;
    }
    Zh->p_couplings[number*2]   = p->cpl[0];
    Zh->p_couplings[number*2+1] = p->cpl[1];
  }
  else {
    if (pb->number<99) {
      if (BS->Sign(pb->number)==-1) {
	if (Zh->m_nprop>nb && nb>=0) Zh->p_propagators[nb].direction = Direction::Outgoing;
      }
    }

    if (p->m==1) {
      Zh->p_arguments[number*2]   = p->number;
      if (!(p->fl).IsTensor()) {
      	Zh->p_arguments[number*2+1] = 99;
	Zh->p_couplings[number*2]   = Complex(1.,0.);
	Zh->p_couplings[number*2+1] = Complex(1.,0.);
	}
	else {
	  Zh->p_arguments[number*2+1]   = p->number;
	  }
    }
    else {
      //incoming boson
      Zh->p_arguments[number*2+1]   = p->number;
      
      if ((p->fl).IsScalar()) {
	Zh->p_arguments[number*2]     = p->number;
	Zh->p_couplings[number*2]   = Complex(0.,0.);
	Zh->p_couplings[number*2+1] = Complex(0.,0.);
      }
      else {
	if (p->fl.IsVector() && p->fl.Mass()!=0.0) Zh->p_arguments[number*2] = p->number+massiveskip;
	else Zh->p_arguments[number*2] = p->number+masslessskip+1;
	Zh->p_couplings[number*2]   = Complex(1.,0.);
	Zh->p_couplings[number*2+1] = Complex(1.,0.);
      }
    }
  }
}

void Zfunc_Generator::Set_Out(Zfunc* Zh,int number,Point* pg,Point* p)
{
  int nb=number;
  if (Zh->m_type=="FFVT"||Zh->m_type=="FFVGS")nb--;
  if (Zh->m_nprop>nb && nb>=0) {
    Zh->p_propagators[nb].numb      = pg->number;
    Zh->p_propagators[nb].kfcode    = (pg->fl).Kfcode();
    Zh->p_propagators[nb].direction = Direction::Outgoing;
  }

  if ((pg->fl).IsScalar() && (pg->left==0 || Zh->m_type=="SSV" || p!=pg)) {
    Zh->p_arguments[number*2]     = pg->number;
    Zh->p_arguments[number*2+1]   = pg->number;
    Zh->p_couplings[number*2]   = Complex(0.,0.);
    Zh->p_couplings[number*2+1] = Complex(0.,0.);
  }
  else {
    if (pg->left!=0) {
      if (pg->m==1 && pg!=p) {
	Zh->p_arguments[number*2]     = pg->number;
	if (!(pg->fl).IsTensor()) {
	  Zh->p_arguments[number*2+1]   = 99;
	  
	  Zh->p_couplings[number*2]   = Complex(1.,0.);
	  Zh->p_couplings[number*2+1] = Complex(1.,0.);
	}
      }
      else {
	Zh->p_arguments[number*2]     = pg->left->number;
	Zh->p_arguments[number*2+1]   = pg->right->number;
	if(pg->middle)if(pg->middle->fl.IsFermion()){
	  if(!(pg->left->fl.IsFermion())) Zh->p_arguments[number*2]   = pg->middle->number;
	  if(!(pg->right->fl.IsFermion()))Zh->p_arguments[number*2+1] = pg->middle->number;
	}
	Zh->p_couplings[number*2]   = pg->cpl[0];
	Zh->p_couplings[number*2+1] = pg->cpl[1];
      }
    }
    else {
      Zh->p_arguments[number*2]     = pg->number;
      if (BS->Sign(pg->number)==-1) {
	Zh->p_arguments[number*2+1]     = pg->number;
	if (pg->fl.IsVector() && pg->fl.Mass()!=0.0) Zh->p_arguments[number*2] = pg->number+massiveskip;
	else Zh->p_arguments[number*2] = pg->number+masslessskip+1;
      }
      else {
	if (pg->fl.IsVector() && pg->fl.Mass()!=0.0) Zh->p_arguments[number*2+1] = pg->number+massiveskip;
	else Zh->p_arguments[number*2+1] = pg->number+masslessskip+1;
      }
      Zh->p_couplings[number*2]   = Complex(1.,0.);
      Zh->p_couplings[number*2+1] = Complex(1.,0.);	
    }
  }
}	
void Zfunc_Generator::Set_Tensor(Zfunc* Zh,Point* p)
{
  Point *pb=p,*pt;
  if(p->fl.IsFermion()){
    if(p->left->fl.IsBoson())pb=p->left;
    if(!(pb->fl.IsTensor())&&p->right->fl.IsBoson())pb=p->right;
    if(!(pb->fl.IsTensor())&&p->middle)
      if(p->middle->fl.IsBoson())pb=p->middle;
  }
  pt=pb;
  if(!(pb->fl.IsTensor())){
    if(pb->left==0)return;
    if(pb->left->fl.IsTensor())pt=pb->left;
    else if(pb->right->fl.IsTensor())pt=pb->right;
    else if(pb->middle)if(pb->middle->fl.IsTensor())pt=pb->middle;
  }else pb=p;
  if(!(pt->fl.IsTensor()))return;
  Zh->p_propagators[Zh->m_nprop-1].numb      = pt->number;
  Zh->p_propagators[Zh->m_nprop-1].kfcode    = (pt->fl).Kfcode();
  Zh->p_propagators[Zh->m_nprop-1].direction = Direction::Outgoing;
  
  int narg=Zh->m_narg - Zh->p_calculator->GetScalarNumb();
  Zh->p_arguments[narg-2]     = pt->number;
  Zh->p_arguments[narg-1]     = pt->number;
  int ic=narg-2;
  if (Zh->m_type=="FFT") Zh->p_couplings[2]=pb->cpl[2];
  if (Zh->m_type=="FFT" || Zh->m_type=="FFVT") ic=0;
  Zh->p_couplings[ic]       = pb->cpl[0];
  Zh->p_couplings[ic+1]     = pb->cpl[1];  
}	

void Zfunc_Generator::Set_FermionProp(Zfunc* Zh,Point* p,Point* pf)
{
  if(Zh->m_nprop!=3)return;
  if(pf){
    int i1=1,i2=2;
    if((p->fl).IsAnti()){i1=2;i2=1;}
    Zh->p_propagators[i2].numb      = p->number;
    Zh->p_propagators[i2].kfcode    = (p->fl).Kfcode();
    Zh->p_propagators[i2].direction = Direction::Incoming;
    if(p->number==0)Zh->p_propagators[i2].direction = Direction::Outgoing;
    Zh->p_propagators[i1].numb      = pf->number;
    Zh->p_propagators[i1].kfcode    = (pf->fl).Kfcode();
    Zh->p_propagators[i1].direction = Direction::Outgoing;
  }
  else{
    Zh->p_propagators[2].numb      = p->left->number;
    Zh->p_propagators[2].kfcode    = (p->left->fl).Kfcode();
    Zh->p_propagators[2].direction = Direction::Outgoing;
    Zh->p_propagators[1].numb      = p->right->number;
    Zh->p_propagators[1].kfcode    = (p->right->fl).Kfcode();
    Zh->p_propagators[1].direction = Direction::Outgoing;
  }
}

void Zfunc_Generator::IsGaugeV(Point* p,int& skal,int& vec)
{
  skal = 0;
  vec  = 0;
  if (p->left!=0) {
    if ((p->fl).IsScalar())        skal++;
    if ((p->left->fl).IsScalar())  skal++;
    if ((p->right->fl).IsScalar()) skal++;
    if ((p->fl).IsTensor())        skal++;
    if ((p->left->fl).IsTensor())  skal++;
    if ((p->right->fl).IsTensor()) skal++;

    if ((p->fl).IsVector())        vec++;
    if ((p->left->fl).IsVector())  vec++;
    if ((p->right->fl).IsVector()) vec++;

    if (p->middle!=0) {
      if ((p->middle->fl).IsScalar()) skal++;
      if ((p->middle->fl).IsTensor()) skal++;
      if ((p->middle->fl).IsVector()) vec++;
    }
  }
}

void Zfunc_Generator::SetDirection(int N,SpinorDirection* spind)
{
  int partner,h,oldp;
  short int i,j;
  short int ip,jp;
  ip = jp = -1;
  for (Zfunc_Iterator zit=zlist.begin();zit!=zlist.end();++zit) {
    Zfunc* z = (*zit);
    for (int pos=0;pos<(z->m_narg - z->p_calculator->GetScalarNumb());pos+=2) {
      int swchange = 0;
      i = pos;
      int first = 0;
      if (z->p_arguments[i]<N) {
	SpinorDirection* sd = spind;
	while (sd) {
	  if (sd->to==z->p_arguments[i]) {
	    first = 1;
	    swchange = 1;
	    ip = i;
	    jp = i+1;
	    // ip<jp
	    break;
	  }
	  sd = sd->Next;
	}
      }
      if (!swchange) {
	i=pos+1;
	if (z->p_arguments[i]<N) {
	  SpinorDirection* sd = spind;
	  while (sd) {
	    if (sd->from==z->p_arguments[i]) {
	      first = 0;
	      swchange = 1;
	      ip = i;
	      jp = i-1;
	      // jp<ip
	      break;
	    }
	    sd = sd->Next;
	  }
	}
      }
      if (swchange) {
	//change
	partner = z->p_arguments[jp];
	oldp    = z->p_arguments[ip];
	//
	h          = z->p_arguments[ip];
	z->p_arguments[ip] = z->p_arguments[jp];
	z->p_arguments[jp] = h;
	if (partner>99) {
	  //Fermionline
	  for (;;) {
	    int end = 0;
	    for (Zfunc_Iterator zit=zlist.begin();zit!=zlist.end();++zit) {
	      Zfunc* zh = (*zit);
	      if (zh!=z) {
		end = 0;
		for (int pos2=0;pos2<(zh->m_narg - zh->p_calculator->GetScalarNumb());pos2+=2) {
		  j = pos2;
		  swchange = 0;
		  if (zh->p_arguments[j]==partner) {
		    if (first==0) {end = 1;break;}
		    ip = j+1;
		    swchange = 1;
		  }
		  else {
		    j = pos2+1;
		    if (zh->p_arguments[j]==partner) {
		      if (first==1) {end = 1;break;}
		      ip = j-1;
		      swchange = 1;
		    }
		  }
		  jp = j;
		  if (swchange) {
		    if (zh->p_arguments[ip] == oldp) break;
		    partner = zh->p_arguments[ip];
		    oldp = zh->p_arguments[jp];
		    if (ip>jp) first = 1;
		          else first = 0; 
		    
		    h           = zh->p_arguments[jp];
		    zh->p_arguments[jp] = zh->p_arguments[ip];
		    zh->p_arguments[ip] = h;
		    break;
		  }
		}
		if (partner<99 || end) break;
	      }
	    }
	    if (partner<99 || end) break;
	  }
	}
      }
    }
  }
}
