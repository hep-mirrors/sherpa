#include "Single_Amplitude_Base.H"
#include "Message.H"
#include "Run_Parameter.H"

using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMEGIC;
using namespace std;

#define Cut_Fermion_Prop

void Single_Amplitude_Base::PrintGraph() 
{
  for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) 
    (*zit)->Print(); 

  msg.Out()<<endl<<endl<<"Propagators: "<<endl;
  for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc* p = *pit;
    msg.Out()<<p->fl<<"("<<p->arg[0]<<")\t --> ";
    for (int i=1;i<p->argnum;i++) msg.Out()<<p->arg[i]<<",";
    msg.Out()<<"on = "<<p->on<<endl;
  }
  msg.Out()<<endl;
}

void Single_Amplitude_Base::ClearCalcList()
{
  //clear z's 
  for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
    Zfunc* z = (*zit);
    /*
    if (!z->value.empty()) {
      z->value.clear();
      for (short int i=0;i<z->calclist.size();i++) {
	delete[] z->calclist[i];
      }
      z->calclist.clear();
    }
    */

    z->ClearCalcList();
  }
}

#define new_calc


int Single_Amplitude_Base::Fill_Args(Zfunc* z, Argument* args, vector<int>* iz, vector<int>* iargs)  
{
  int k=-1;
  for (short int i=0;i<z->narg;i++) {
    args[2*i].numb = z->arg[i];
    int hit=0,j;
    if(iz!=0){
      for(j=0;j<iz->size();j++){
	if(iabs((*iz)[j])==z->arg[i]){
	  hit=1;
	  break;
	}
      }
    }
    if(hit){
      if((*iargs)[2*j+1]>100){
	//cout<<"Tensor: "<<(*iargs)[2*j+1]<<":"<<j<<","<<i<<endl;
	//z->Print();
	k= j;           //spin2 tensor dummies have to be replaced
      }
      if((z->arg[i]>99)&&(z->arg[i]<200)){
	//fermions
	if((*iz)[j]<0) args[2*i].numb=(*iargs)[2*j]; // new cut fermion prop  
	else{ 
	  if ((*iargs)[2*j]<0) {                  // old cut fermion prop
	    if (i%2==0) args[2*i].spinortype = Spinor::vbar;
	           else args[2*i].spinortype = Spinor::v;
	  }
	}      
      }

      if(z->arg[i]<99) 
	if (fl[z->arg[i]].Majorana() && i%2!=0) args[2*i].spinortype = Spinor::v;
      args[2*i+1].numb = (*iargs)[2*j+1];

      if(i<z->narg-1)
	if(z->arg[i]==z->arg[i+1]&&(*iz)[j]==(*iz)[j+1]){ // spin2 boson
	  i++;j++;
	  args[2*i].numb = z->arg[i];
	  args[2*i+1].numb = (*iargs)[2*j+1];
	}

    }
    else{
      if (z->arg[i]<20) {                                //old external massless Vector Boson treatment
	args[2*i].numb = z->arg[i]-10;
	for(short int j=0;j<iz->size();j++){
	  if(iabs((*iz)[j])==z->arg[i]-10){
	    args[2*i+1].numb = (*iargs)[2*j+1];
	    break;
	  }
	}
      }
      else{
	if (z->arg[i]>=20 && z->arg[i]<99) {
	  //old  massive Vector bosons
	  args[2*i].numb   = z->arg[i]-20;
	  args[2*i+1].numb = -1;
	}
        //else args[2*i+1].numb = 0;
      }
    }
  }
  return k;
}

Kabbala Single_Amplitude_Base::Single_ZvalueTensor(Zfunc* z,vector<int>* iz, vector<int>* iargs,int k)
{
  //cout<<"in Single_ZvalueTensor"<<endl;
  Kabbala value;
  int narg=z->narg - z->calculator->GetScalarNumb();
  if(z->arg[narg-2]!=z->arg[narg-1]){
    cout<<"Single_Amplitude_Base::Single_ZvalueTensor: Unexpected tensor sign! "<<(*iargs)[2*k+1]<<" "<<k<<endl;
    z->Print();
    abort();
  }
  vector<vector<int> > pol;
  vector<int> sign;
  tensor_struc ts;
  int tensor_type=(*iargs)[2*k+1];
  ts.Get_Pol_Combos(tensor_type,&pol,&sign);
  //for(int i=0;i<pol.size();i++)
  //cout<<sign[i]<<"("<<pol[i][0]<<","<<pol[i][1]<<") ";cout<<endl;
  //cout<<"TViarg:"<<(*iargs)[2*k+1]<<","<<(*iargs)[2*k+3]<<endl;
  for(int i=0;i<pol.size();i++){
    (*iargs)[2*k+1]=pol[i][0];
    (*iargs)[2*k+3]=pol[i][1];
    if(sign[i]==-1) value -= Single_Zvalue(z,iz,iargs);
    else            value += Single_Zvalue(z,iz,iargs);
  }
  (*iargs)[2*k+1]=tensor_type;
  (*iargs)[2*k+3]=tensor_type;
  return value;
}

Kabbala Single_Amplitude_Base::Single_ZGroup_value(Zfunc* z,
					     vector<int>* iz,              //list of indizes(propagators to sum)
					     vector<int>* iargs,int last)
{ 
  Kabbala value;
  if(z->GetOp()=='+') { 
     for (int i=0;i<z->GetSize();i++) {
      Kabbala hlp = Single_Zvalue((*z)[i],iz,iargs);
      if ((*z)[i]->pn>0) hlp*= Get_Prop((*z)[i]->pn,(*z)[i]->psnew,(*z)[i]->GetOp());

      if (z->GetSign(i)==-1) value -= hlp;
                        else value += hlp;
    }
  }  

  if(z->GetOp()=='*'){
    Kabbala hlp;

    if(z->GetSize()!=2){
      cout<<"Invalid Zfunc_ProdGroup!"<<endl;
      abort();
    }
    vector<int> iz_s;
    iz_s.push_back(z->GetSumIndex());

    vector<int> dummy;
    vector<vector<int> > iargs_s;

    iargs_s.reserve(2);iargs_s.resize(2,dummy);

    SetLoopVar(iz_s,iargs_s);
    iz->push_back(iz_s[0]);

    bool tensor=Get_Pflav(z->GetSumIndex())->IsTensor();
    if(tensor) iz->push_back(iz_s[0]);

    if(!tensor){                                 //fermions and cutted vector bosons
      for(int i1=0;i1<iargs_s[0].size();i1++)
	for(int i2=0;i2<iargs_s[1].size();i2++){
	  iargs->push_back(iargs_s[0][i1]);
	  iargs->push_back(iargs_s[1][i2]);
	  
	  hlp= Single_Zvalue((*z)[0],iz,iargs)*
	       Single_Zvalue((*z)[1],iz,iargs);
	  
	  if ((z->GetSumIndex()>99) && (z->GetSumIndex()<199))
	    hlp*= Single_Mass_Terms((*iz)[iz->size()-1],(*iargs)[iargs->size()-2]);
	  
	  value+=hlp;
	  iargs->pop_back();iargs->pop_back();
	}
    }
    else {                                              //cutted spin2 bosons notation: hep-ph/9811350
      for(int i1=0;i1<iargs_s[1].size();i1++)
	for(int i2=i1;i2<iargs_s[1].size();i2++){
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i1]);
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i2]);

	  hlp= Single_Zvalue((*z)[0],iz,iargs)*
	       Single_Zvalue((*z)[1],iz,iargs);
#ifdef Kabbala_on
	  if(i1!=i2)hlp*= Kabbala(string("2"),Complex(2.,0.));
#else
	  if(i1!=i2)hlp*2.;
#endif
	  value+=hlp;
 
	  iargs->pop_back();iargs->pop_back();
	  iargs->pop_back();iargs->pop_back();
	}
      Kabbala hlp1,hlp2;
      for(int i1=0;i1<iargs_s[1].size();i1++){
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i1]);
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i1]);

	  hlp1+= Single_Zvalue((*z)[0],iz,iargs);
	  hlp2+= Single_Zvalue((*z)[1],iz,iargs);

	  iargs->pop_back();iargs->pop_back();
	  iargs->pop_back();iargs->pop_back();
      }
#ifdef Kabbala_on
      Kabbala n33;
      if(buildstring) n33=(shand->Get_Generator())->Get_Enumber(Complex(1./3.,0.));
      else n33=Kabbala(string(""),Complex(1./3.,0.));
      value -= n33*hlp1*hlp2;

#else
      value-=(1./3.)*hlp1*hlp2;
#endif

      iz->pop_back();
    }

    iz->pop_back();

  }
  return value;
}

Kabbala Single_Amplitude_Base::Single_Zvalue(Zfunc* z,vector<int>* iz, vector<int>* iargs,int last)
{
  if(last && z->GetSize()>1) return Single_ZGroup_value(z,iz,iargs,last);
   
  Argument* args = new Argument[2*z->narg];
  int k=Fill_Args(z,args,iz,iargs);
  if(k!=-1&&z->GetSize()==1){
    delete[] args;
    return Single_ZvalueTensor(z,iz,iargs,k); //produce tensors for extern spin2 bosons
  }

#ifdef new_calc
  if (z->Equal!=z) {
    //cout<<"In not equal part!!!!"<<endl;
    Zfunc* ze = z->Equal;
    for (short int i=0;i<ze->value.size();i++) {
      //cout<<i<<": "<<ze->value[i].String()<<endl;
      int hit = 1;
      //Signs Equal ??
      for (short int j=1;j<2*ze->narg;j+=2) {
	if (ze->calclist[i][j]!=args[j]) {
	  hit = 0;
	  break;
	}
      }
      if (hit) {
	//Arguments equal ?
	for (short int j=0;j<2*ze->narg;j+=2) {
	  if (ze->calclist[i][j]!=args[j]) {
	    hit = 0;
	    break;
	  }
	}
      }

      if (hit) {
	delete[] args;
	if (z->sign==-1) return -ze->value[i];
	            else return ze->value[i];
      }
    }
  }
  else {
    //cout<<"In else part: "<<z->value.size()<<endl;
    //Check signlist
    for (short int i=0;i<z->value.size();i++) {
      int hit = 1;
      for (short int j=0;j<2*z->narg;j++) {
	if (z->calclist[i][j]!=args[j]) {
	  hit = 0;
	  break;
	}
      }
      if (hit) {
	delete[] args;
	if (z->sign==-1) return -z->value[i];
	            else return z->value[i];
      }
    }
  }
#endif
  Argument* newargs = new Argument[2*z->narg];
  for (short int i=0;i<2*z->narg;i++) newargs[i] = args[i];

  Kabbala value;

  if(z->GetSize()>1)value=Single_ZGroup_value(z,iz,iargs,last);
  else {
    z->calculator->SetArgCouplProp(2*z->narg,args,z->coupl,z->pn,z->psnew,&plist);

    value = z->calculator->Do();
  }

#ifdef Kabbala_on
  if (buildstring) {
    if ( (value.String()).find(string("+"))!=-1 ||
	 (value.String()).find(string("-"))!=-1 ||
	 (value.String()).find(string("*"))!=-1 )
      value = (shand->Get_Generator())->Get_CZnumber(value.Value(),value.String());
  }
  else value.SetString("");

#endif
#ifdef new_calc  
  if (z->Equal!=z) {
    z->Equal->value.push_back(value);
    z->Equal->calclist.push_back(newargs);
  }
  else {
    z->value.push_back(value);
    z->calclist.push_back(newargs);
  }

#endif
  delete[] args;
  if (z->sign==-1) return -value;

  return value;
}

void Single_Amplitude_Base::SetLoopVar(vector<int>& iz,vector<vector<int> >& iargs)
{
  for (short int i=0;i<iz.size();i++) if(iz[i]>99){
    //cout<<"New Prop: "<<iz[i]<<endl;
    if (iz[i]<199) {
#ifdef Cut_Fermion_Prop
      Pfunc* p;
      for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
	p = *pit;
	if (p->arg[0]==iz[i]) break;	
      }
      double mass = 0.;
      //cout<<iz[i]<<";"<<p->fl<<" --> ";
      for (short j=1;j<p->argnum;j++) {
	//cout<<p->arg[j]<<","<<fl[p->arg[j]]<<","<<b[p->arg[j]]<<" ; ";
	if (fl[p->arg[j]].IsAnti()) mass -= b[p->arg[j]] * fl[p->arg[j]].Mass();
	                       else mass += b[p->arg[j]] * fl[p->arg[j]].Mass();
      }
      //cout<<endl;

      double particlemass = p->fl.Mass();

      if (p->fl.IsAnti()) particlemass = -particlemass;

      if (AMATOOLS::IsEqual(particlemass,mass) && (p->fl).Width()==0.) {
	for (short j=1;j<p->argnum;j++)	iargs[2*i].push_back(p->arg[j]);
	//mark this special propagator as negative
	iz[i] = -iz[i];
	//cout<<"Propagator "<<iz[i]<<" cut new"<<endl;
      } 
      else {
	//fermions
	iargs[2*i].push_back(-1);
	iargs[2*i].push_back(1);
	//cout<<"Propagator "<<iz[i]<<" cut old"<<endl;
      }
#else
      //fermions
      iargs[2*i].push_back(-1);
      iargs[2*i].push_back(1);
#endif
      //helicity of fermions
      iargs[2*i+1].push_back(-1);
      iargs[2*i+1].push_back(1);
    }
    else {
      //bosons
      iargs[2*i].push_back(0);
      //helicity of bosons
      BS->PropPolarisation(iz[i],plist,iargs[2*i+1]);
    }
  }
  /*
  for (short int i=0;i<2*iz.size();i++) { 
    cout<<"0: ";
    if (!i%2) cout<<iz[i/2];
    cout<<" --> ";
    for (short int j=0;j<iargs[i].size();j++) 
      cout<<iargs[i][j]<<";";
    cout<<endl;
  }
  */
}

Complex Single_Amplitude_Base::Zvalue(String_Handler * sh, int ihel) 
{
  if (sh) return sh->Zvalue(amplnumber,ihel);
  msg.Error()<<" ERROR in Single_Amplitude_Base::Zvalue(String_Handler * sh, int ihel) "<<endl;
  shand->Zvalue(amplnumber,ihel); 
}

Complex Single_Amplitude_Base::Zvalue(int ihel,int* signlist) 
{
  if (signlist==0) return shand->Zvalue(amplnumber,ihel);

    //cout<<"Amplitude number "<<amplnumber;
    //cout<<"============================================================================"<<endl;
    
    return Zvalue_sum(ihel,signlist);
}

Flavour* Single_Amplitude_Base::Get_Pflav(int pn)
{
  for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc*  p = *pit;
    if(pn==p->arg[0])return &(p->fl);
  }
  msg.Error()<<"Single_Amplitude_Base::Get_Pflav: Propagator not found!"<<endl;
  abort();
}

Kabbala Single_Amplitude_Base::Get_Prop(int pn, Argument *ps, char op)
{
#ifdef Kabbala_on
  Kabbala Pols(string("1"),Complex(1.,0.));
#else
  Kabbala Pols(1.,0.);
#endif
  Basic_Pfunc bpf(shand->Get_Generator(),BS);

  int cnt=0,sign=1;
  for(int i=0;i<pn;i++){

    for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
      Pfunc*  p = *pit;
      int pc;

      if (!ps[i].maped) pc = p->arg[0];
      else {
	pc = p->momnum;
	if((p->fl).Kfcode()!=ps[i].kfcode) pc = -1;
      }
      //cout<<ps[i].numb<<"/"<<pc<<" OP:"<<op<<" On:"<<p->on<<endl;

      if (pc==iabs(ps[i].numb)) {
	if(op==0 && p->on==1)break;
	if (p->haspol && p->fl.IsVector()) sign*=-1;
        Pols*= bpf.P(p);
	cnt++;
	break;
      } 
    }
  }
  if(sign<0)Pols=-Pols;

#ifdef Kabbala_on
  if (buildstring) {
    if (cnt>1)
      Pols = (shand->Get_Generator())->Get_CZnumber(Pols.Value(),Pols.String());
  }
#endif

  return Pols;
}


Kabbala Single_Amplitude_Base::Single_Mass_Terms(int iz,int iarg)
{  
#ifdef Kabbala_on
  Kabbala factor(string("1"),Complex(1.,0.));
#else
  Kabbala factor(1.,0.);
#endif    
  if (iabs(iz)>=199)return factor;
  if (iz<0) if (b[iarg]<0) {
    factor= -factor;
    //cout<<"Factor "<<factor.String()<<endl;
    return factor;
  }    

  Basic_MassTermfunc bmtf(shand->Get_Generator(),BS);

  bmtf.SetArgCouplProp(0,0,0,0,0,&plist); 
  //cout<<"Single_Mass_Terms_new:"<<iz<<" "<<iarg;
  if(iz>0) {
    factor = bmtf.MassTerm(iz*iarg);
#ifdef Kabbala_on
    factor*= Kabbala(string("0.5"),Complex(1./2.,0.));
#else
    factor*= 0.5;
#endif
  }
  return factor;
}



Complex Single_Amplitude_Base::Zvalue_sum(int ihel, int* signlist) 
{
  vector<int> iz;

  Zfunc* z;

  if(zlist.size()>1){
    for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
      Pfunc* p = *pit;
      if (p->on) {
	if (p->arg[0]>99) iz.push_back(p->arg[0]);
      }
    }
  }

  //cout<<"Zvalue_new: iz"<<endl;
  
  /*cout<<"Number of propagators: "<<iz.size()<<endl;
  for(int i=0;i<iz.size();i++) cout<<iz[i]<<" ";
   cout<<endl;
  */
  
  if (!iz.empty()) {
    vector<int> dummy;
    vector<vector<int> > iargs;iargs.reserve(2*iz.size());iargs.resize(2*iz.size(),dummy);

    SetLoopVar(iz,iargs);
    
    vector<vector<int> > indexlist;
    int cnt=0;
    for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
      Zfunc* z = (*zit);
      indexlist.push_back(dummy);
      for (short int i=0;i<z->narg;i++) {
	if (z->arg[i]>99) {
	  for (short int j=0;j<iz.size();j++) {
	    if (iabs(iz[j])==z->arg[i]) {
	      indexlist[cnt].push_back(j);
	      break;
	    }
	  }
	}
      }
      cnt++;
    }
   int over;

    do{
      over=1;
      int min=1000000;
      int imin;
      for(short int i=0;i<iz.size();i++){
	int adds=0;
	for(cnt=0;cnt<indexlist.size();cnt++)
	  for(short int j=0;j<indexlist[cnt].size();j++) {
	    over=0;
	    if(indexlist[cnt][j]==i){
	      if(adds==0)adds=1;
	      for(short int k=0;k<indexlist[cnt].size();k++)
		if(indexlist[cnt][k]!=i)
		  {
		    adds*=iargs[2*indexlist[cnt][k]].size()*iargs[2*indexlist[cnt][k]+1].size();
		  }
	      break;
	    }
	  }      
	if (adds>0&&adds<min) {min=adds;imin=i;}
      }
      int ia=0;
      if(over==0){
	Zfunc* zh[2]; 
	indexlist.push_back(dummy);

	vector<vector<int> >::iterator ilt=indexlist.begin();	
	for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();) {
	  int hit=0;
	  for(short int j=0;j<(*ilt).size();j++) {
	    //cout<<j<<":"<<(*ilt)[j]<<endl;
	    if((*ilt)[j]==imin){
	      //cout<<"ini "<<ia<<endl;
	      zh[ia]=(*zit);
	      ia++;
	      for(short int k=0;k<(*ilt).size();k++)
		if((*ilt)[k]!=imin)
		  {
		    indexlist[indexlist.size()-1].push_back((*ilt)[k]);
		  }
	      ilt=indexlist.erase(ilt);
	      zit=zlist.erase(zit);
	      hit=1;
	      break;
	    }
	  }
	  //if(ia==2)break;
	  if(!hit){++ilt;++zit;}
	}
	if(ia!=2){
	  cout<<"Error Zvalue_new: index appeared "<<ia<<" times!"<<endl;
	  zh[0]->Print();zh[1]->Print();
	  abort();
	}
	
        Zfunc_Group *sf=new Zfunc_Group(*zh[0],*zh[1],iabs(iz[imin]),&plist);
	zlist.push_back(sf);

      }
    }while(over==0);
  }

#ifdef Kabbala_on
  Kabbala value(string("1"),Complex(1.,0.));
#else
  Kabbala value(1.,0.);
#endif
  vector<int> iarg,ize;
  for(int j=0;j<N;j++){
    ize.push_back(j);
    iarg.push_back(0);
    iarg.push_back(signlist[j]);
    if(signlist[j]>100){
      ize.push_back(j);
      iarg.push_back(0);
      iarg.push_back(signlist[j]);
    }
  }

  for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
    Zfunc* z = (*zit);

    value*= Single_Zvalue(z,&ize,&iarg,1);
    if (z->pn>0) value*= Get_Prop(z->pn,z->psnew,z->GetOp());
  }
  //cout<<"Amplitude vor prop:"<<amplnumber<<": "<<value.Value()<<endl;

  Basic_Pfunc bpf(shand->Get_Generator(),BS);

  //   scalar propagators:
#ifndef Scalar_Args
  for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc*  p = *pit;
    if ((p->fl).IsScalar()) {
      value *= bpf.P(p);
    } 
  }
#endif

  if (sign<0) value = -value; 
 
  //if (buildstring) cout<<"Amplitude nach prop:"<<amplnumber<<": "<<value.String()<<endl;
#ifdef Kabbala_on
  if (buildstring) shand->Set_String(amplnumber,ihel,value.String());
  return value.Value();
#else
  return value;
#endif      
}
