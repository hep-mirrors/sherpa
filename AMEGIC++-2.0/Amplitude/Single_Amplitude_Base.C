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
    for (int i=0;i<z->GetSize();i++) {
      Zfunc* z1 = (*z)[i];
      if (!z1->value.empty()) {
	z1->value.clear();
	for (short int i=0;i<z1->calclist.size();i++) {
	  delete[] z1->calclist[i];
	}
	z1->calclist.clear();
      }
    }
  }
}

#define new_calc

Kabbala Single_Amplitude_Base::Mass_Terms(vector<int>& iz,vector<int>& ii,
				     vector<vector<int> >& iargs)
{
#ifdef Kabbala_on
  Kabbala factor(string("1"),Complex(1.,0.));
#else
  Kabbala factor(1.,0.);
#endif    

  Basic_MassTermfunc bmtf(shand->Get_Generator(),BS);

  bmtf.SetArgCouplProp(0,0,0,0,0,&plist); 

  for (short int i=0;i<iz.size();i++) {
    if (iz[i]<199 && iz[i]>0) factor *= bmtf.MassTerm(iz[i]*iargs[2*i][ii[2*i]]);
  }

  return factor;
}

void Single_Amplitude_Base::Fill_Args(Zfunc* z,Argument* args,int* signlist,
				      vector<int>* iz,              //list of indizes(propagators to sum)
				      vector<int>* ii,              //current positions of indizes 
				      vector<vector<int> >* iargs)  //list of indizes with different index values
                                      //e.g (*iargs)[3][(*ii)[3]] is the current value of the third index 
{
  int nargs = 0;
  for (short int i=0;i<z->narg;i++) {
    //cout<<"Setting Argument : "<<z->arg[i]<<endl;
    if (z->arg[i]>99) {
      if (z->arg[i]>199) {
	//boson
	args[nargs].numb   = z->arg[i];
	//scalar bosons have no helicity, result should be independent of special choice....
	nargs++;
	args[nargs].numb = 1;
	if (iz!=0) {
	  for (short int j=0;j<iz->size();j++) {
	    if ((*iz)[j]==z->arg[i]) {
	      args[nargs].numb = (*iargs)[2*j+1][(*ii)[2*j+1]];
	      break;
	    }
	  }
	}
      }
      else {
	//fermion prop
	for (short int j=0;j<iz->size();j++) {
	  if (iabs((*iz)[j])==z->arg[i]) {
	    if ((*iz)[j]>0) {
	      args[nargs].numb = (*iz)[j];
	      if ((*iargs)[2*j][(*ii)[2*j]]<0) {
		if (i%2==0) args[nargs].spinortype = Spinor::vbar;
                       else args[nargs].spinortype = Spinor::v;
	      }
	    }
	    else {
	      //new case
	      args[nargs].numb = (*iargs)[2*j][(*ii)[2*j]];
	    }
	    args[++nargs].numb = (*iargs)[2*j+1][(*ii)[2*j+1]];
	    break;
	  }
	}
      }
    }
    else {
      args[nargs].numb   = z->arg[i];
      if ((args[nargs].numb>N-1) && (args[nargs].numb<20)) {
	for (short int j=z->arg[i]-10-1;j>=0;j--) {
	  if ( (fl[j].IsVector() && AMATOOLS::IsZero(fl[j].Mass())) ||
	       fl[j]==Flavour(kf::pol) ) {
	    args[++nargs].numb = signlist[j];	    
	    break;
	  }
	}	      
	args[nargs-1].numb = z->arg[i]-10;
      }
      else {
	if (args[nargs].numb>=20 && args[nargs].numb<99) {
	  //massive Vector bosons
	  args[nargs].numb   = z->arg[i]-20;
	  args[++nargs].numb = -1;
	}
	else {
	  if (z->arg[i]!=99) {
	    if (fl[z->arg[i]].Majorana() && i%2!=0) args[nargs].spinortype = Spinor::v;
	    args[++nargs].numb = signlist[z->arg[i]];
	  }
	  else nargs++; 
	}
      }
    }
    nargs++;
  }
}

Kabbala Single_Amplitude_Base::Single_Zvalue(Argument* args,Zfunc* z)
{ 
#ifdef new_calc
  if (z->Equal!=z) {
    //cout<<"In not equal part!!!!"<<endl;
    Zfunc* ze = z->Equal;
    for (short int i=0;i<ze->value.size();i++) {
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
	if (z->sign==-1) return -z->value[i];
	            else return z->value[i];
      }
    }
  }
#endif
  Argument* newargs = new Argument[2*z->narg];
  for (short int i=0;i<2*z->narg;i++) newargs[i] = args[i];

  z->calculator->SetArgCouplProp(2*z->narg,args,z->coupl,z->pn,z->psnew,&plist);
  Kabbala value = z->calculator->Do();

#ifdef Kabbala_on
  if (buildstring) {
    if ( (value.String()).find(string("+"))!=-1 ||
	 (value.String()).find(string("-"))!=-1 ||
	 (value.String()).find(string("*"))!=-1 )
      value = (shand->Get_Generator())->Get_CZnumber(value.Value(),value.String());
  }

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
  if (z->sign==-1) return -value;

  return value;
}

void Single_Amplitude_Base::SetLoopVar(vector<int>& iz,vector<vector<int> >& iargs)
{
  for (short int i=0;i<iz.size();i++) {
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
	if (fl[p->arg[j]].IsAnti()) mass -= b[p->arg[j]]*rpa.consts.Mass(fl[p->arg[j]],sqr(rpa.gen.Ecms()));
	                       else mass += b[p->arg[j]]*rpa.consts.Mass(fl[p->arg[j]],sqr(rpa.gen.Ecms()));
      }
      //cout<<endl;

      double particlemass = rpa.consts.Mass(p->fl,sqr(rpa.gen.Ecms()));

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

  // use new sum scheme
  return Zvalue_new_sum(ihel,signlist);

  Kabbala Zval;

  vector<int> iz;

  for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc* p = *pit;
    if ((p->on && (p->fl).IsFermion()) ||
	p->haspol) {
      if (p->arg[0]>99) iz.push_back(p->arg[0]);
    }
  }

  //cout<<"Number of fermions: "<<iz.size()<<endl;

  Kabbala sum;
  Argument* args;
  
  if (iz.empty()) {
#ifdef Kabbala_on
    Kabbala factor(string("1"),Complex(1.,0.));
#else
    Kabbala factor(1.,0.);
#endif
    for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
      Zfunc* z = (*zit);
      args = new Argument[2*z->narg];
      Fill_Args(z,args,signlist);      
      factor *= Single_Zvalue(args,z);
      delete[] args;
    }
    sum += factor;
  }
  else {
    vector<int> ii;ii.reserve(2*iz.size());ii.resize(2*iz.size(),0);
    vector<int> dummy;
    vector<vector<int> > iargs;iargs.reserve(2*iz.size());iargs.resize(2*iz.size(),dummy);

    SetLoopVar(iz,iargs);

    int over = 0;
    int sw1 = 0;
    do {
#ifdef Kabbala_on
      Kabbala factor(string("1"),Complex(1.,0.));
#else
      Kabbala factor(1.,0.);
#endif
      int minus = 0;
      for (int j=0;j<iz.size();j++) {
	if (iabs(iz[j])<199 && iz[j]<0) {
	  int currval = iargs[2*j][ii[2*j]];
	  if (b[currval]<0) minus++;
	}
      }

      /*
      cout<<"Next: ";
      for (short int i=0;i<2*iz.size();i++)
	cout<<iargs[i][ii[i]]; 
      cout<<endl;
      */
      for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
	Zfunc* z = (*zit);
	args = new Argument[2*z->narg];
	Fill_Args(z,args,signlist,&iz,&ii,&iargs);
	factor *= Single_Zvalue(args,z);
	delete[] args;
      }
      factor *= Mass_Terms(iz,ii,iargs);
      if (minus%2) sum -= factor;
              else sum += factor;
      for (short int j=2*iz.size()-1;j>=0;j--) {
	ii[j]++;
	if (ii[j]==iargs[j].size()) {
	  ii[j]=0;
	  if (j==0) over = 1;
	}
	else break;
      }
    } while (over==0);

#ifdef Kabbala_on
    for (short int k=0;k<iz.size();k++) {
      if (iz[k]<199 && iz[k]>0) sum *= Kabbala(string("0.5"),Complex(1./2.,0.));
    }
#else
    for (short int k=0;k<iz.size();k++) {
      if (iz[k]<199 && iz[k]>0) sum *= 0.5;
    }
#endif    
  }
  
  //Propagators...
  Basic_Pfunc bpf(shand->Get_Generator(),BS);

#ifdef Kabbala_on
  Kabbala Pols(string("1"),Complex(1.,0.));
#else
  Kabbala Pols(1.,0.);
#endif
  for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc* p = *pit;
    if (p->on) {
      //cout<<"Multiply outer with "<<p->arg[0]<<" = "<<p->fl<<endl;
      if (p->haspol) Pols *= -bpf.P(p);
                else Pols *= bpf.P(p); 
    }
  }

#ifdef Kabbala_on
  if (buildstring) {
    if ( (Pols.String()).find(string("+"))!=-1 ||
	 (Pols.String()).find(string("-"))!=-1 ||
	 (Pols.String()).find(string("*"))!=-1 )
      Pols = (shand->Get_Generator())->Get_CZnumber(Pols.Value(),Pols.String());
  }
#endif

  sum *= Pols;

  if (sign<0) sum = -sum; 
    
#ifdef Kabbala_on
  if (buildstring) shand->Set_String(amplnumber,ihel,sum.String());

  return sum.Value();
#else
  return sum;
#endif      
}

Kabbala Single_Amplitude_Base::Single_Mass_Terms_new(int iz,int iarg)
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
  if(iz>0) factor = bmtf.MassTerm(iz*iarg);
  return factor;
}

Kabbala Single_Amplitude_Base::Generate_Block(vector<Block_Info>& BlockList,int bi,
					 vector<int>& iz,vector<int>& ii,vector<vector<int> >& iargs,
					 int* signlist)
{
  Kabbala sum,hlp;
  if(BlockList[bi].l>-1){
    //cout<<"Generate_Block: index:"<<BlockList[bi].index<<" #:"<<bi<<" l:"<<BlockList[bi].l<<" r:"<<BlockList[bi].r<<endl;
    vector<int> iii(ii);
    for(short int i1=0;i1<iargs[2*BlockList[bi].index].size();i1++)
      for(short int i2=0;i2<iargs[2*BlockList[bi].index+1].size();i2++){
	iii[2*BlockList[bi].index]=i1;
	iii[2*BlockList[bi].index+1]=i2;

	hlp= Generate_Block(BlockList,BlockList[bi].l,iz,iii,iargs,signlist)
	    *Generate_Block(BlockList,BlockList[bi].r,iz,iii,iargs,signlist);
	//hlp = (shand->Get_Generator())->Get_CZnumber(hlp.Value(),hlp.String());

	if(BlockList[bi].index<iz.size())
	  sum+=hlp*Single_Mass_Terms_new(iz[BlockList[bi].index],iargs[2*BlockList[bi].index][i1]);    
	else
	  sum+=hlp;    
      }
    //cout<<"Generate_Block: index:"<<BlockList[bi].index<<" #:"<<bi<<" l:"<<BlockList[bi].l<<" r:"<<BlockList[bi].r<<endl;
    //cout<<sum.String()<<endl;
    if (buildstring) sum = (shand->Get_Generator())->Get_CZnumber(sum.Value(),sum.String());
  }
  else{
    Argument* args = new Argument[2*BlockList[bi].Z->narg];
    Fill_Args(BlockList[bi].Z,args,signlist,&iz,&ii,&iargs);
    sum = Single_Zvalue(args,BlockList[bi].Z);
    delete[] args;
    //cout<<"Generate_Block: Z:"<<bi<<" >"<<sum.String()<<";";
    //for(int i=0;i<ii.size();i++)cout<<" "<<ii[i];
    //cout<<endl;
  }
  return sum;
}

Complex Single_Amplitude_Base::Zvalue_new_sum(int ihel, int* signlist) 
{
  Kabbala Zval;

  vector<int> iz;

  Zfunc* z;

  for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc* p = *pit;
    if ((p->on && (p->fl).IsFermion()) ||
	p->haspol ) {
      if (p->arg[0]>99) iz.push_back(p->arg[0]);
    }
  }

  //cout<<"Zvalue_new: iz"<<endl;

  //cout<<"Number of propagators: "<<iz.size()<<endl;

  Kabbala sum;
  Argument* args;
  
  if (iz.empty()) {
#ifdef Kabbala_on
    Kabbala factor(string("1"),Complex(1.,0.));
#else
    Kabbala factor(1.,0.);
#endif
    for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
      Zfunc* z = (*zit);
      args = new Argument[2*z->narg];
      Fill_Args(z,args,signlist);      
      factor *= Single_Zvalue(args,z);
      delete[] args;
    }
    sum += factor;
  }
  else {
    vector<int> dummy;
    vector<vector<int> > iargs;iargs.reserve(2*iz.size());iargs.resize(2*iz.size(),dummy);

    SetLoopVar(iz,iargs);
    
    vector<vector<int> > indexlist;
    vector<Zfunc*> ScalarZ;
    int cnt=0;
    Block_Info bdummy;
    bdummy.l = -1;bdummy.r = -1;bdummy.index = -1;
    vector<Block_Info> BlockList;
    for (list<Zfunc*>::iterator zit=zlist.begin();zit!=zlist.end();++zit) {
      Zfunc* z = (*zit);
      indexlist.push_back(dummy);
      bdummy.Z=z;
      BlockList.push_back(bdummy);
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
      if(indexlist[cnt].size()==0){
	//	msg.Error()<<"Zvalue_new: Z-function without sum index found"<<endl;
	ScalarZ.push_back(z);
      }
      cnt++;
    }
    /*  cout<<"Zvalue_new: indexlist"<<endl;
    for(cnt=0;cnt<indexlist.size();cnt++){
      cout<<cnt<<": ";
      for(int j=0;j<indexlist[cnt].size();j++)cout<<indexlist[cnt][j]<<" ";
      cout<<endl;
    }
    cout<<"Zvalue_new: iargs"<<endl;
    for(cnt=0;cnt<iargs.size();cnt++){
      cout<<cnt<<": ";
      for(int j=0;j<iargs[cnt].size();j++)cout<<iargs[cnt][j]<<" ";
      cout<<endl;
      }*/
    int over;
    bdummy.Z=0;
    vector<int> TreeList;

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
	//cout<<"Zvalue_new: sum over prop "<<imin<<endl;
	indexlist.push_back(dummy);
	bdummy.index=imin;
	BlockList.push_back(bdummy);
	//if(min==1)TreeList.push_back(BlockList.size()-1);
	for(cnt=0;cnt<indexlist.size()-1;cnt++){
	  for(short int j=0;j<indexlist[cnt].size();j++) {
	    if(indexlist[cnt][j]==imin){
	      ia++;
	      for(short int k=0;k<indexlist[cnt].size();k++)
		if(indexlist[cnt][k]!=imin)
		  {
		    indexlist[indexlist.size()-1].push_back(indexlist[cnt][k]);
		  }
	      indexlist[cnt].clear();
	      switch(ia){
	      case 1:BlockList[BlockList.size()-1].l=cnt;break;
	      case 2:BlockList[BlockList.size()-1].r=cnt;break;		
	      }
	      break;
	    }
	  }
	  //if(ia==2)break;
	}
	if(indexlist[indexlist.size()-1].empty())TreeList.push_back(BlockList.size()-1);
	if(ia!=2) msg.Error()<<"Error Zvalue_new: index appeared "<<ia<<" times!"<<endl;
      }
    }while(over==0);

    int l=TreeList.size();
    //if(l!=1)cout<<"TreeList.size:"<<l<<endl;
    if(l>1){
      dummy.push_back(0);
      iargs.push_back(dummy);
      iargs.push_back(dummy);
    }
    while(l>1){
      BlockList.push_back(bdummy);
      BlockList[BlockList.size()-1].l=TreeList[l-1];
      BlockList[BlockList.size()-1].r=TreeList[l-2];
      BlockList[BlockList.size()-1].index=iargs.size()/2-1;
      TreeList[l-2]=BlockList.size()-1;
      l--;
    }

    /*   if(amplnumber>83){
 cout<<"Zvalue_new: BlockList"<<endl;
    for(cnt=0;cnt<BlockList.size();cnt++){
      cout<<cnt<<": ";
      if(BlockList[cnt].Z)cout<<"Z"<<endl;
      else cout<<" l:"<<BlockList[cnt].l<<" r:"<<BlockList[cnt].r<<" index:"<<BlockList[cnt].index<<endl;
      }
      }*/
    vector<int> ii;ii.reserve(iargs.size());ii.resize(iargs.size(),0);
    int bi=BlockList.size()-1;
    if(BlockList[bi].l>-1){
      for(short int i1=0;i1<iargs[2*BlockList[bi].index].size();i1++)
	for(short int i2=0;i2<iargs[2*BlockList[bi].index+1].size();i2++){
	  ii[2*BlockList[bi].index]=i1;
	  ii[2*BlockList[bi].index+1]=i2;

	  Kabbala hlp=
	    Generate_Block(BlockList,BlockList[bi].l,iz,ii,iargs,signlist)*
	    Generate_Block(BlockList,BlockList[bi].r,iz,ii,iargs,signlist);

	  //hlp = (shand->Get_Generator())->Get_CZnumber(hlp.Value(),hlp.String());
	  if(BlockList[bi].index<iz.size())
	    sum+=hlp*Single_Mass_Terms_new(iz[BlockList[bi].index],iargs[2*BlockList[bi].index][i1]);    
	  else
	    sum+=hlp;    
		}
    }
    //cout<<sum.String()<<endl;    

#ifdef Kabbala_on
      Kabbala factor(string("1"),Complex(1.,0.));
#else
      Kabbala factor(1.,0.);
#endif

    for(short int i=0;i<ScalarZ.size();i++){
      args = new Argument[2*ScalarZ[i]->narg];
      Fill_Args(ScalarZ[i],args,signlist);      
      factor *= Single_Zvalue(args,ScalarZ[i]);
      delete[] args;
    }
    if(ScalarZ.size()>1) {
      if (buildstring) factor = (shand->Get_Generator())->Get_CZnumber(factor.Value(),factor.String());
    }

    sum*=factor;

#ifdef Kabbala_on
    for (short int k=0;k<iz.size();k++) {
      if (iz[k]<199 && iz[k]>0) sum *= Kabbala(string("0.5"),Complex(1./2.,0.));
    }
#else
    for (short int k=0;k<iz.size();k++) {
      if (iz[k]<199 && iz[k]>0) sum *= 0.5;
    }
#endif    
  }
  
  //Propagators...
  Basic_Pfunc bpf(shand->Get_Generator(),BS);

#ifdef Kabbala_on
  Kabbala Pols(string("1"),Complex(1.,0.));
#else
  Kabbala Pols(1.,0.);
#endif
  for (list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc*  p = *pit;
    if (p->on) {
      if (p->haspol) Pols *= -bpf.P(p);
      else Pols *= bpf.P(p);
    } 
    /*
    if (p->haspol) Pols *= -1;
    if (p->on) Pols *= bpf.P(p);
    */
  }

#ifdef Kabbala_on
  if (buildstring) {
    if ( (Pols.String()).find(string("+"))!=-1 ||
	 (Pols.String()).find(string("-"))!=-1 ||
	 (Pols.String()).find(string("*"))!=-1 )
      Pols = (shand->Get_Generator())->Get_CZnumber(Pols.Value(),Pols.String());
  }
#endif

  sum *= Pols;

  if (sign<0) sum = -sum; 
    
#ifdef Kabbala_on
  if (buildstring) shand->Set_String(amplnumber,ihel,sum.String());
  return sum.Value();
#else
  return sum;
#endif      

  /*
  p = Plist;
  
  Basic_Pfunc bpf(shand->Get_Generator(),BS);

#ifdef Kabbala_on
  Kabbala Pols(string("1"),Complex(1.,0.));
#else
  Kabbala Pols(1.,0.);
#endif

  while (p) {
    if (p->haspol) Pols *= -bpf.P(p);
              else Pols *= bpf.P(p);
    p = p->Next;
  }

#ifdef Kabbala_on
  if (buildstring) {
    if ( (Pols.String()).find(string("+"))!=-1 ||
	 (Pols.String()).find(string("-"))!=-1 ||
	 (Pols.String()).find(string("*"))!=-1 )
      Pols = (shand->Get_Generator())->Get_CZnumber(Pols.Value(),Pols.String());
  }
#endif

  sum *= Pols;

  if (sign<0) sum = -sum; 
 
  //cout<<"done! "<<ihel<<" "<<igraph<<endl;
   
#ifdef Kabbala_on
  if (buildstring) shand->Set_String(igraph,ihel,sum.String());

  return sum.Value();
#else
  return sum;
#endif
  */      
}
