#include "Zfunc_Generator.H"
#include "Vector.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace std;

#define lorentz_type

extern int iabs(int&);

void Zfunc_Generator::BuildZlist(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS)
{
  zcalc.push_back(new Y_Calc(_sgen,_BS));
  zcalc.push_back(new Z_Calc(_sgen,_BS));
  zcalc.push_back(new VVV_Calc(_sgen,_BS));
  zcalc.push_back(new VVVV_Calc(_sgen,_BS));
  zcalc.push_back(new V4_Calc(_sgen,_BS));
  zcalc.push_back(new G4_Calc(_sgen,_BS));
  zcalc.push_back(new V5_Calc(_sgen,_BS));
  zcalc.push_back(new V4V3_Calc(_sgen,_BS));
  zcalc.push_back(new VVS_Calc(_sgen,_BS));
  zcalc.push_back(new SSV_Calc(_sgen,_BS));
  zcalc.push_back(new SSS_Calc(_sgen,_BS));
  zcalc.push_back(new VVSS_Calc(_sgen,_BS));
  zcalc.push_back(new ZNC_Calc(_sgen,_BS));
  zcalc.push_back(new VVVNC_Calc(_sgen,_BS));  
}

void Zfunc_Generator::LorentzConvert(Point* p)
{
  if (p==0) return;
  Lorentz_Function* l = p->Lorentz;

  for (short int i=0;i<l->NofIndex();i++) {
    switch (l->partarg[i]) {
    case 0: l->partarg[i] = p->number;break;
    case 1: l->partarg[i] = p->left->number;break;
    case 2: l->partarg[i] = p->right->number;break;
    case 3: l->partarg[i] = p->middle->number;break;
    }
  }
  LorentzConvert(p->right);
  LorentzConvert(p->left);  
  LorentzConvert(p->middle);
}

void Zfunc_Generator::MarkCut(Point* p,int notcut,bool fromfermion)
{
  if (p==0) return; 

  if (p->fl.isvector() && p->number>99){
    p->m = 1;
    notcut++;
    //cout<<"+++++++++Cut:"<< p->number<<endl;
    if(fromfermion&&p->left->fl.isfermion()){
      //      cout<<"MarkCut: "<<p->number<<" not cutted"<<endl;
      p->m=0;
      }
    if(AMATOOLS::IsZero(p->fl.mass())){
      p->m=0;
      }	
    /*if(!(p->left->fl.kfcode()==kf::gluon&&p->right->fl.kfcode()==kf::gluon&&!fromfermion)){
      cout<<"MarkCut: "<<p->number<<" not cutted (1)"<<endl;
      p->m=0;
      //notcut=true;
      }*/
    /* if(!(p->fl.kfcode()==kf::photon||p->fl.kfcode()==kf::Z)){
      cout<<"MarkCut: "<<p->number<<" not cutted (1)"<<endl;
      p->m=0;
      }*/
  }
  else p->m = 0;

  // to aktivate Cut of Boson lines comment out the following line
  //p->m = 0;

  // "new gauge test" cut all massless propagators
  if (p->fl.isvector() && p->number>99  && rpa.me.CutScheme()==1) {
    if(AMATOOLS::IsZero(p->fl.mass())){
      p->m=1;
      }	

  }


  MarkCut(p->right,notcut,p->fl.isfermion());
  MarkCut(p->left,notcut,p->fl.isfermion());
  MarkCut(p->middle,notcut,p->fl.isfermion()); 
}

void Zfunc_Generator::Convert(Point* p)
{
  Zfunc* Zh = 0;
  if ((p->left==0) && (p->right==0)) return;
 
  if ( p->fl.isfermion() || p->fl.isscalar() || 
       (p->fl.isvector() && p->number<99) || p->m==1) {
    Zh = new Zfunc;
    Point* pb;
    Point* pf;
    pb = pf = 0;
    
    if (p->fl.isfermion()) {
      if (p->left->fl.isboson()) {
	pb = p->left;
	pf = p->right;
      }
      if (p->right->fl.isboson()) {
	pb = p->right;
	pf = p->left;
      }
    }
    else {
      //Incoming Boson
      pb = p;
    }
    if(!LFDetermine_Zfunc(Zh,p,pf,pb)){
      cout<<"Unknown lorentz sequence found! Cutting...";
      Point* ph=pb->right;
      //cout<<"right "<<ph->fl<<endl;
      if (!( ph->fl.isfermion() || ph->fl.isscalar() || 
	     (ph->fl.isvector() && ph->number<99) || ph->m==1)&&ph->left)
	if(!(ph->left->fl.isfermion()))
	  {ph->m=1;cout<<"right"<<endl;}
      else {
	ph=pb->left;
	//cout<<"left "<<ph->fl<<endl;
	if (!( ph->fl.isfermion() || ph->fl.isscalar() || 
	       (ph->fl.isvector() && ph->number<99) || ph->m==1)&&ph->left)
	  if(!(ph->left->fl.isfermion()))
	    {ph->m=1;cout<<"left"<<endl;}
	else if(p->middle){
	  ph=pb->middle;
	  //cout<<"middle "<<ph->fl<<endl;
	  if (!( ph->fl.isfermion() || ph->fl.isscalar() || 
		 (ph->fl.isvector() && ph->number<99) || ph->m==1)&&ph->left)
	    if(!(ph->left->fl.isfermion()))
	      {ph->m=1;cout<<"middle"<<endl;}
	}
      }
      Convert(p); 
      return;
    }
  }
  if (Zh) zlist.push_back(Zh);

  Convert(p->right);
  Convert(p->left);
  if (p->middle) Convert(p->middle);
}

void Zfunc_Generator::Lorentz_Sequence(Point* pb,vector<Lorentz_Function> &lflist)
{ 
  if (pb->left==0 && pb->fl.isscalar()) return;
 
  lflist.push_back(*(pb->Lorentz));

  //Endpoints are scalar or non-bosons
  int skal,vec;
  IsGaugeV(pb,skal,vec);      
  
  if (skal+vec>=3) {
    if (pb->left->fl.isvector() && pb->left->m==0)  Lorentz_Sequence(pb->left,lflist);
    if (pb->right->fl.isvector() && pb->right->m==0) Lorentz_Sequence(pb->right,lflist);
    if (pb->middle!=0) {
            if (pb->middle->fl.isvector() && pb->middle->m==0) Lorentz_Sequence(pb->middle,lflist);
      if (pb->middle->m==1) {
	Lorentz_Function lf(lf::Pol);
	lf.SetParticleArg(pb->middle->number);
	lflist.push_back(lf);
      }
    }

    if (pb->left->m==1) {
      Lorentz_Function lf(lf::Pol);
      lf.SetParticleArg(pb->left->number);
      lflist.push_back(lf);
    }
    if (pb->right->m==1) {
      Lorentz_Function lf(lf::Pol);
      lf.SetParticleArg(pb->right->number);
      lflist.push_back(lf);
    }
  }
}

void Zfunc_Generator::LFPrint(const vector<Lorentz_Function> &lflist)
{
  AORGTOOLS::msg.Out()<<"LorentzList: "<<endl;
  for (short int i=0;i<lflist.size();i++)
    AORGTOOLS::msg.Out()<<lflist[i].String(1)<<endl;
  AORGTOOLS::msg.Out()<<endl;
}

void Zfunc_Generator::LFPrint(const vector<Lorentz_Function*> &lflist)
{
  AORGTOOLS::msg.Out()<<"LorentzList: "<<endl;
  for (short int i=0;i<lflist.size();i++)
    AORGTOOLS::msg.Out()<<lflist[i]->String(1)<<endl;
  AORGTOOLS::msg.Out()<<endl;
}

lf::code Zfunc_Generator::LFEff(lf::code type)
{ 
  if (AORGTOOLS::rpa.me.Model()==AORGTOOLS::Model_Type::NCQED) return (type==lf::Pol) ? lf::Gamma_NC : type;
  return (type==lf::Pol) ? lf::Gamma : type;
}

int Zfunc_Generator::LFDetermine_Zfunc(Zfunc* Zh,Point* p,Point* pf,Point* pb)
{
  Zh->type = -10;
  vector<Lorentz_Function> lflist;

  if (pf!=0) lflist.push_back(*(p->Lorentz));

  if (pf==0 && pb->fl.isvector()) {
    Lorentz_Function lf(lf::Pol);
    lf.SetParticleArg(pb->number);
    lflist.push_back(lf);
  }
  
  if (!( (pb->fl.isscalar() || pb->m==1) && pf!=0)) Lorentz_Sequence(pb,lflist);
  else {
    if (pb->m==1 && pf!=0) {
      Lorentz_Function lf(lf::Pol);
      lf.SetParticleArg(pb->number);
      lflist.push_back(lf);
    }
  }
  //LFPrint(lflist);  

  for (short int i=0;i<zcalc.size();i++) {
    if (lflist.size()==(zcalc[i]->lorentzlist).size()) {
      int hit = 1; 
      vector<int> typerem;
      
      for (short int j=0;j<lflist.size();j++) {
	int hit2 = 1;
	for (short int k=0;k<typerem.size();k++) {
	  if (typerem[k]==LFEff(lflist[j].type)) {
	    hit2 = 0;
	    break;
	  }
	}
	if (hit2) {
	  //counting 
	  int type1 = 0;
	  for (short int k=j;k<lflist.size();k++) {
	    if (LFEff(lflist[j].type)==LFEff(lflist[k].type)) type1++;
	  }
	  int type2 = 0;
	  for (short int k=0;k<(zcalc[i]->lorentzlist).size();k++) {
	    if (LFEff(lflist[j].type)==LFEff((zcalc[i]->lorentzlist[k]).type)) type2++;
	  }  
	  if (type1!=type2) {
	    hit = 0;
	    break;
	  }
	  else typerem.push_back(LFEff(lflist[j].type));
	}
      }
      if (hit) {
	Zh->type       = zcalc[i]->type;
	Zh->calculator = zcalc[i];
	break;
      }
    }
  }
  if (Zh->type==-10) {
    return 0;
    AORGTOOLS::msg.Error()<<"No Lorentzfunction found!"<<endl;
    LFPrint(lflist);  
    abort();
  }

  LFFill_Zfunc(Zh,lflist,p,pf,pb);
  return 1;
}


void Zfunc_Generator::CopyOrder(vector<Lorentz_Function> &lflist,vector<Lorentz_Function*> & lfpointer)
{
  for (short int i=0;i<lflist.size();i++) lfpointer.push_back(&lflist[i]);

  for (short int i=0;i<lfpointer.size();i++) 
    for (short int j=i+1;j<lfpointer.size();j++) {
      if (lfpointer[i]->NofIndex()<lfpointer[j]->NofIndex()) {
	Lorentz_Function* help;
	help         = lfpointer[i];
	lfpointer[i] = lfpointer[j];
	lfpointer[j] = help;
      }
    }  
  //LFPrint(lfpointer);
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
  
  for (short int i=0;i<lfpointer.size();i++) {
    for (short int k=0;k<lfpointer[i]->NofIndex();k++) {
      int lfarg = abs(lfpointer[i]->partarg[k]);
      int caarg = abs(capointer[i]->partarg[k]);
      
      int hit = 1;
      
      for (short int j=0;j<numbcount;j++) {
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

void Zfunc_Generator::LFFill_Zfunc(Zfunc* Zh,vector<Lorentz_Function> &lflist,Point* p,Point* pf,Point* pb)
{
  vector<Lorentz_Function*> lfpointer;
  CopyOrder(lflist,lfpointer);
  //LFPrint(lfpointer);
  vector<Lorentz_Function*> capointer;
  CopyOrder(Zh->calculator->lorentzlist,capointer);
  //LFPrint(capointer);

  vector<Lorentz_Function*> permpointer;

  for (short int j=0;j<lfpointer.size();j++) {
    permpointer.push_back(lfpointer[j]);
    permpointer[j]->InitPermutation();
  }
  
  int* lfnumb = new int[Zh->calculator->pn];
  int* canumb = new int[Zh->calculator->pn];

  //loop over all permutations
  for (;;) {
    int i = Compare(Zh->calculator->pn,lfpointer,lfnumb,capointer,canumb);
    if (i==lfpointer.size()) break;

    //loop over previous permutations
    for (;;) {
      int typecount = 0;
      int typemin   = 1000; 
      int typemax   = 0;

      for (short int j=0;j<lfpointer.size();j++) {
	if (LFEff(lfpointer[j]->type)==LFEff(lfpointer[i]->type)) {
	  if (typemin>j) typemin = j;
	  if (typemax<j) typemax = j;
	  typecount++;
	}
      }
      
      vector<Lorentz_Function*> copypointer;
      for (short int j=0;j<lfpointer.size();j++) copypointer.push_back(lfpointer[j]);

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
	    i = Compare(Zh->calculator->pn,copypointer,lfnumb,capointer,canumb);
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
      else i = Compare(Zh->calculator->pn,lfpointer,lfnumb,capointer,canumb);

      if (i==lfpointer.size()) break;
	
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
	  AORGTOOLS::msg.Error()<<"Error in Zfunc_Generator::LFFill_Zfunc()!!"<<endl;
	  abort();
	}
	//LFPrint(lfpointer);
      }
    }
    if (i==lfpointer.size()) break;
  }

  //Total Sign......
  Zh->sign = 1;

  for (short int j=0;j<permpointer.size();j++) Zh->sign *= permpointer[j]->GetSign();
  
  SetPropDirection(Zh->calculator->pn,pb->number,lfpointer,lfnumb,capointer,canumb);


  //Setting the arguments

  Zh->narg   = Zh->calculator->narg;
  Zh->ncoupl = Zh->calculator->ncoupl;
  Zh->pn     = Zh->calculator->pn;

  Zh->arg    = NULL;
  Zh->psnew  = NULL;

  if (Zh->narg>0)   Zh->arg = new int[Zh->narg];
  if (Zh->ncoupl>0) Zh->coupl = new Complex[Zh->ncoupl];
  if (Zh->pn>0)     Zh->psnew = new Argument[Zh->pn];
  
  //Set all Couplings Zero
  for (short int i=0;i<Zh->ncoupl;i++) Zh->coupl[i] = Complex(0.,0.);

  for (short int i=0;i<Zh->calculator->pn;i++) {
    if (lfnumb[i]==pb->number) {
      Set_In(Zh,canumb[i],p,pf,pb);
      break;
    }
  }
  //Special cases
  switch (Zh->type) {
    
  case zl::Y:
        if (pf==0) Set_Out(Zh,0,pb,p);
              else Set_In(Zh,0,p,pf,pb);
    break;
  case zl::Z:Set_Out(Zh,1,pb,p);break;
  case zl::ZNC:Set_Out(Zh,1,pb,p);break;
  default:
    int icoupl        = Zh->narg;
    Zh->coupl[icoupl] = pb->cpl[1];icoupl++;

    SetArgs(Zh,lfnumb,canumb,pb->left,p,icoupl);
    SetArgs(Zh,lfnumb,canumb,pb->right,p,icoupl);
    SetArgs(Zh,lfnumb,canumb,pb->middle,p,icoupl);
  }

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
  for (short int i=0;i<lfpointer.size();i++) {
    if (LFEff(lfpointer[i]->type)==lf::Gamma) {
      for (short int k=0;k<lfpointer[i]->NofIndex();k++) {
	if (lfpointer[i]->partarg[k]==incoming) {
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
  for (short int i=0;i<lfpointer.size();i++) {
    if (i!=position) {
      for (short int k=0;k<lfpointer[i]->NofIndex();k++) {
	if (lfpointer[i]->partarg[k]==incoming) {
	  start = i;
	  break;
	}
      }
      if (start!=-1) break;
    }
  }

  if (start==-1) return;

  for (short int k=0;k<lfpointer[position]->NofIndex();k++) {
    if (lfpointer[position]->partarg[k]==incoming) {
      if (capointer[position]->partarg[k]<0) {
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
    if (lfpointer[start]->partarg[k]!=incoming) 
      SearchNextProp(Nargs,lfpointer,lfnumb,capointer,canumb,lfpointer[start]->partarg[k],start);
  }
}

void Zfunc_Generator::SetArgs(Zfunc* Zh,int* lfnumb,int* canumb,Point* pb,Point* p,int& icoupl)
{
  if (pb==0) return;
  //if (!pb->fl.isvector()) return;


  for (short int i=0;i<Zh->calculator->pn;i++) {
    if (lfnumb[i]==pb->number) {
      if (pb->number<99 || pb->fl.isscalar() || pb->m==1) Set_Out(Zh,canumb[i],pb,p);
      else {
	if  (!pb->left->fl.isvector() && !pb->right->fl.isvector()) Set_Out(Zh,canumb[i],pb,p);
	else {
	  Zh->coupl[icoupl]                   = pb->cpl[1];icoupl++;
	  Zh->psnew[abs(canumb[i])].numb      = pb->number;
	  Zh->psnew[abs(canumb[i])].kfcode    = (pb->fl).kfcode();
	  Zh->psnew[abs(canumb[i])].direction = Direction::Outgoing;
	  if (canumb[i]<0) 
	    Zh->psnew[abs(canumb[i])].direction = Direction::Incoming;
	  SetArgs(Zh,lfnumb,canumb,pb->left,p,icoupl);
	  SetArgs(Zh,lfnumb,canumb,pb->right,p,icoupl);
	  SetArgs(Zh,lfnumb,canumb,pb->middle,p,icoupl);
	}
      }
      break;
    }
  }
}

void Zfunc_Generator::Set_In(Zfunc* Zh,int number, Point* p, Point* pf,Point* pb)
{
  //  cout<<"Incoming: "<<number<<" : "<<pb->number<<endl;
  Zh->psnew[number].numb      = pb->number;
  Zh->psnew[number].kfcode    = (pb->fl).kfcode();
  Zh->psnew[number].direction = Direction::Incoming;
  
  if (pf!=0) {
    if (pb->m==1) {
      //      cout<<"Changing Directions...."<<number<<endl;
      Zh->psnew[number].direction = Direction::Outgoing;
    }
    if (pb->number<99) {
      if (BS->Sign(pb->number)==1) {
	Zh->psnew[number].direction = Direction::Outgoing;
      }
    }

    if ((p->fl).isanti()) {
      Zh->arg[number*2+1] = p->number;
      Zh->arg[number*2]   = pf->number;
    }
    else {
      Zh->arg[number*2]   = p->number;
      Zh->arg[number*2+1] = pf->number;
    }
    Zh->coupl[number*2]   = p->cpl[0];
    Zh->coupl[number*2+1] = p->cpl[1];
  }
  else {
    if (pb->number<99) {
      if (BS->Sign(pb->number)==-1) {
	Zh->psnew[number].direction = Direction::Outgoing;
      }
    }

    if (p->m==1) {
      Zh->arg[number*2]   = p->number;
      //Marker for -99
      Zh->arg[number*2+1] = 99;
      Zh->coupl[number*2]   = Complex(1.,0.);
      Zh->coupl[number*2+1] = Complex(1.,0.);
    }
    else {
      //incoming boson
      Zh->arg[number*2+1]   = p->number;
      
      if ((p->fl).isscalar()) {
	Zh->arg[number*2]     = p->number;
	Zh->coupl[number*2]   = Complex(0.,0.);
	Zh->coupl[number*2+1] = Complex(0.,0.);
      }
      else {
	if (p->fl.isvector() && !AMATOOLS::IsZero(p->fl.mass())) Zh->arg[number*2] = p->number+20;
	                                                    else Zh->arg[number*2] = p->number+10+1;
	Zh->coupl[number*2]   = Complex(1.,0.);
	Zh->coupl[number*2+1] = Complex(1.,0.);
      }
    }
  }
}

void Zfunc_Generator::Set_Out(Zfunc* Zh,int number,Point* pg,Point* p)
{
  //setting an outgoing in gauge vertices
  if (Zh->pn>number) {
    Zh->psnew[number].numb      = pg->number;
    Zh->psnew[number].kfcode    = (pg->fl).kfcode();
    Zh->psnew[number].direction = Direction::Outgoing;
  }

  if ((pg->fl).isscalar() && pg->left==0) {
    Zh->arg[number*2]     = pg->number;
    Zh->arg[number*2+1]   = pg->number;
    Zh->coupl[number*2]   = Complex(0.,0.);
    Zh->coupl[number*2+1] = Complex(0.,0.);
  }
  else {
    if (pg->left!=0) {
      if (pg->m==1 && pg!=p) {
	Zh->arg[number*2]     = pg->number;
	//Marker for -99
	Zh->arg[number*2+1]   = 99;
	Zh->coupl[number*2]   = Complex(1.,0.);
	Zh->coupl[number*2+1] = Complex(1.,0.);
      }
      else {
	Zh->arg[number*2]     = pg->left->number;
	Zh->arg[number*2+1]   = pg->right->number;
	Zh->coupl[number*2]   = pg->cpl[0];
	Zh->coupl[number*2+1] = pg->cpl[1];
      }
    }
    else {
      Zh->arg[number*2]     = pg->number;
      if (BS->Sign(pg->number)==-1) {
	Zh->arg[number*2+1]     = pg->number;
	if (pg->fl.isvector() && !AMATOOLS::IsZero(pg->fl.mass())) Zh->arg[number*2] = pg->number+20;
	                                                      else Zh->arg[number*2] = pg->number+10+1;
      }
      else {
	if (pg->fl.isvector() && !AMATOOLS::IsZero(pg->fl.mass())) Zh->arg[number*2+1] = pg->number+20;
	else Zh->arg[number*2+1] = pg->number+10+1;
      }
      Zh->coupl[number*2]   = Complex(1.,0.);
      Zh->coupl[number*2+1] = Complex(1.,0.);	
    }
  }
}	

void Zfunc_Generator::IsGaugeV(Point* p,int& skal,int& vec)
{
  skal = 0;
  vec  = 0;
  if (p->left!=0) {
    if ((p->fl).isscalar())        skal++;
    if ((p->left->fl).isscalar())  skal++;
    if ((p->right->fl).isscalar()) skal++;

    if ((p->fl).isvector())        vec++;
    if ((p->left->fl).isvector())  vec++;
    if ((p->right->fl).isvector()) vec++;

    if (p->middle!=0) {
      if ((p->middle->fl).isscalar()) skal++;
      if ((p->middle->fl).isvector()) vec++;
    }
  }
}

void Zfunc_Generator::SetDirection(int N,SpinorDirection* spind)
{
  int partner,h,oldp;
  short int i,j;
  short int ip,jp;
  ip = jp = -1;
  for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
    Zfunc* z = (*zit);
    for (int pos=0;pos<z->narg;pos+=2) {
      int swchange = 0;
      i = pos;
      int first;
      if (z->arg[i]<N) {
	SpinorDirection* sd = spind;
	while (sd) {
	  if (sd->to==z->arg[i]) {
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
	if (z->arg[i]<N) {
	  SpinorDirection* sd = spind;
	  while (sd) {
	    if (sd->from==z->arg[i]) {
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
	partner = z->arg[jp];
	oldp    = z->arg[ip];
	//
	h          = z->arg[ip];
	z->arg[ip] = z->arg[jp];
	z->arg[jp] = h;
	if (partner>99) {
	  //Fermionline
	  for (;;) {
	    int end;
	    for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
	      Zfunc* zh = (*zit);
	      if (zh!=z) {
		end = 0;
		for (int pos2=0;pos2<zh->narg;pos2+=2) {
		  j = pos2;
		  swchange = 0;
		  if (zh->arg[j]==partner) {
		    if (first==0) {end = 1;break;}
		    ip = j+1;
		    swchange = 1;
		  }
		  else {
		    j = pos2+1;
		    if (zh->arg[j]==partner) {
		      if (first==1) {end = 1;break;}
		      ip = j-1;
		      swchange = 1;
		    }
		  }
		  jp = j;
		  if (swchange) {
		    if (zh->arg[ip] == oldp) break;
		    partner = zh->arg[ip];
		    oldp = zh->arg[jp];
		    if (ip>jp) first = 1;
		          else first = 0; 
		    
		    h           = zh->arg[jp];
		    zh->arg[jp] = zh->arg[ip];
		    zh->arg[ip] = h;
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
