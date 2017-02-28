#include "AMEGIC++/Amplitude/Single_Amplitude_Base.H"
#include "AMEGIC++/Main/Tools.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include "AMEGIC++/Amplitude/Zfunc.H"
#include "AMEGIC++/Amplitude/Zfunctions/Zfunc_Calc.H"
#include "AMEGIC++/String/String_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <cassert>

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

#define Cut_Fermion_Prop

Single_Amplitude_Base::Single_Amplitude_Base(int* _b,int _n, Basic_Sfuncs* _BS,
					     ATOOLS::Flavour* _fl,
					     String_Handler* _shand) 
  : b(_b), N(_n), shand(_shand), BS(_BS), fl(_fl) 
{
  zlist= new Zfunc_List;
}

Single_Amplitude_Base::Single_Amplitude_Base(String_Handler* _shand, int an) 
  : shand(_shand) 
{
  zlist = NULL;
  amplnumber = an;
}

Zfunc_List* Single_Amplitude_Base::GetZlist() {
  return zlist;
}

Pfunc_List* Single_Amplitude_Base::GetPlist() {
  return &plist;
}

int Single_Amplitude_Base::GetSign() {
  return sign;
}

void Single_Amplitude_Base::SetSign(int s) {
  sign=s;
}

void Single_Amplitude_Base::SetNumber(int& i) {
  amplnumber = i;i++;
}

Amplitude_Base* Single_Amplitude_Base::GetAmplitude(const int n) {
  return (n==amplnumber) ? this : 0;
}
 
Single_Amplitude_Base::~Single_Amplitude_Base() 
{
  if (zlist){
    for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) delete (*zit);
    delete zlist;
  }
}

void Single_Amplitude_Base::PrintGraph() 
{
  if (!msg_LevelIsTracking()) return;
  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) 
    (*zit)->Print(); 

  msg_Out()<<endl<<endl<<"Propagators: "<<endl;
  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc* p = *pit;
    msg_Out()<<p->fl<<"("<<p->arg[0]<<")\t --> ";
    for (int i=1;i<p->argnum;i++) msg_Out()<<p->arg[i]<<",";
    msg_Out()<<"on = "<<p->on<<endl;
  }
  msg_Out()<<endl;
}

void Single_Amplitude_Base::ClearCalcList()
{
  //clear z's 
  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
    Zfunc* z = (*zit);

    z->ClearCalcList();
  }
}

void Single_Amplitude_Base::KillZList()
{
  if (!zlist) return;    
  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
    (*zit)->KillZList();
    delete (*zit);
  }
  delete zlist;
  zlist = NULL;
}

int Single_Amplitude_Base::FillArgs(Zfunc* z, int* args, vector<int>* iz, vector<int>* iargs)  
{
  int k=-1;
  for (short int i=0;i<z->m_narg;i++) {
    args[2*i] = z->p_arguments[i];
    int hit=0,j;
    if(iz!=0){
      for(j=0;j<(int)iz->size();j++){
	if(iabs((*iz)[j])==z->p_arguments[i]){
	  hit=1;
	  break;
	}
      }
    }
    if(hit){
      if((*iargs)[2*j+1]>100){
	k= j;           //spin2 tensor dummies have to be replaced
      }
      args[2*i+1] = (*iargs)[2*j+1];
      if((z->p_arguments[i]>99)&&(z->p_arguments[i]<200)){
	//fermions
	if((*iz)[j]<0) args[2*i]=(*iargs)[2*j]; // new cut fermion prop  
	else{ 
	  if ((*iargs)[2*j]<0) {                  // old cut fermion prop
	    args[2*i]   *= -1;
	    args[2*i+1] *= -1;
	  }
	}      
      }
    
      if(z->p_arguments[i]<99) {
	if (((fl[z->p_arguments[i]].Majorana() || 
	      //final-state line
	      (IsChargino(fl[z->p_arguments[i]]) && b[z->p_arguments[i]]==1 && !fl[z->p_arguments[i]].IsAnti()) ||
	      //initial-state line
	      (IsChargino(fl[z->p_arguments[i]]) && b[z->p_arguments[i]]==-1 && fl[z->p_arguments[i]].IsAnti())
	      )
	     && i%2!=0 ) || 
	    ((
	      //final-state line
	      (IsChargino(fl[z->p_arguments[i]]) && b[z->p_arguments[i]]==1 && fl[z->p_arguments[i]].IsAnti()) || 
	      //initial state line
	      (IsChargino(fl[z->p_arguments[i]]) && b[z->p_arguments[i]]==-1 && !fl[z->p_arguments[i]].IsAnti()))
	     
	     && i%2==0)) //args[2*i].spinortype = Spinor::v;
	  {
	    args[2*i]   *= -1;args[2*i+1] *= -1;
	  }
      }
      if(i<z->m_narg-1)
	if(z->p_arguments[i]==z->p_arguments[i+1]&&(*iz)[j]==(*iz)[j+1]){ // spin2 boson
	  i++;j++;
	  args[2*i] = z->p_arguments[i];
	  args[2*i+1] = (*iargs)[2*j+1];
	}
    }
    else{
      if (z->p_arguments[i]<massiveskip) {                                //old external massless Vector Boson treatment
	assert(iz!=NULL);
	args[2*i] = z->p_arguments[i]-masslessskip;
	for(size_t j=0;j<iz->size();j++){
	  if(iabs((*iz)[j])==z->p_arguments[i]-masslessskip){
	    args[2*i+1] = (*iargs)[2*j+1];
	    break;
	  }
	}
      }
      else{
	if (z->p_arguments[i]>=massiveskip && z->p_arguments[i]<99) {
	  //old  massive Vector bosons
	  args[2*i]   = z->p_arguments[i]-massiveskip;
	  args[2*i+1] = -1;
	}
	else args[2*i+1] = 0;
      }
    }
  }
  return k;
}

Kabbala Single_Amplitude_Base::SingleZvalueTensor(Zfunc* z,vector<int>* iz, vector<int>* iargs,int k)
{
  Kabbala value;
  int narg=z->m_narg - z->p_calculator->GetScalarNumb();
  if(z->p_arguments[narg-2]!=z->p_arguments[narg-1]){
    msg_Error()<<"ERROR in Single_Amplitude_Base::SingleZvalueTensor: "<<std::endl
	       <<"   Unexpected tensor sign! "<<(*iargs)[2*k+1]<<" "<<k<<endl;
    z->Print();
    Abort();
  }
  vector<vector<int> > pol;
  vector<int> sign;
  Tensor_Struc ts;
  int tensor_type=(*iargs)[2*k+1];
  ts.GetPolCombos(tensor_type,&pol,&sign);
  for(size_t i=0;i<pol.size();i++){
    (*iargs)[2*k+1]=pol[i][0];
    (*iargs)[2*k+3]=pol[i][1];
    if(sign[i]==-1) value -= SingleZvalue(z,iz,iargs);
    else            value += SingleZvalue(z,iz,iargs);
  }
  (*iargs)[2*k+1]=tensor_type;
  (*iargs)[2*k+3]=tensor_type;
  return value;
}

Kabbala Single_Amplitude_Base::SingleZGroupvalue(Zfunc* z,
						 vector<int>* iz,//list of indizes(propagators to sum)
						 vector<int>* iargs,int last)
{ 
  Kabbala value;
  if(z->GetOp()=='+') { 
     for (int i=0;i<z->GetSize();i++) {
      Kabbala hlp = SingleZvalue((*z)[i],iz,iargs);
      if ((*z)[i]->m_nprop>0) hlp*= GetProp((*z)[i]);

      if (z->GetSign(i)==-1) value -= hlp;
                        else value += hlp;
    }
  }  
  if(z->GetOp()=='*'){
    Kabbala hlp;

    if(z->GetSize()!=2){
      msg_Error()<<"ERROR in Single_Amplitude_Base::SingleZGroupvalue : "<<std::endl
		 <<"   Invalid Zfunc_ProdGroup. Abort the run."<<endl;
      Abort();
    }
    vector<int> iz_s;
    iz_s.push_back(z->GetSumIndex());

    vector<int> dummy;
    vector<vector<int> > iargs_s;

    iargs_s.reserve(2);iargs_s.resize(2,dummy);

    SetLoopVar(iz_s,iargs_s);
    iz->push_back(iz_s[0]);

    bool tensor=GetPflav(z->GetSumIndex())->IsTensor();
    if(tensor) iz->push_back(iz_s[0]);

    if(!tensor){                                 //fermions and cutted vector bosons
      for(size_t i1=0;i1<iargs_s[0].size();i1++)
	for(size_t i2=0;i2<iargs_s[1].size();i2++){
	  iargs->push_back(iargs_s[0][i1]);
	  iargs->push_back(iargs_s[1][i2]);
	  
	  hlp= SingleZvalue((*z)[0],iz,iargs)*
	       SingleZvalue((*z)[1],iz,iargs);

	  if ((z->GetSumIndex()>99) && (z->GetSumIndex()<199))
	    hlp*= SingleMassTerms((*iz)[iz->size()-1],(*iargs)[iargs->size()-2]);
	  value+=hlp;
	  iargs->pop_back();iargs->pop_back();
	}
    }
    else {                                              //cutted spin2 bosons notation: hep-ph/9811350
      for(size_t i1=0;i1<iargs_s[1].size();i1++)
	for(size_t i2=i1;i2<iargs_s[1].size();i2++){
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i1]);
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i2]);

	  hlp= SingleZvalue((*z)[0],iz,iargs)*
	       SingleZvalue((*z)[1],iz,iargs);

	  if(i1!=i2)hlp*= Kabbala(string("2"),Complex(2.,0.));
	  value+=hlp;
 
	  iargs->pop_back();iargs->pop_back();
	  iargs->pop_back();iargs->pop_back();
	}
      Kabbala hlp1,hlp2;
      for(size_t i1=0;i1<iargs_s[1].size();i1++){
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i1]);
	  iargs->push_back(iargs_s[0][0]);
	  iargs->push_back(iargs_s[1][i1]);

	  hlp1+= SingleZvalue((*z)[0],iz,iargs);
	  hlp2+= SingleZvalue((*z)[1],iz,iargs);

	  iargs->pop_back();iargs->pop_back();
	  iargs->pop_back();iargs->pop_back();
      }

      Kabbala n33;
      if(buildstring) n33=(shand->Get_Generator())->GetEnumber(Complex(1./3.,0.));
      else n33=Kabbala(string(""),Complex(1./3.,0.));
      value -= n33*hlp1*hlp2;

      iz->pop_back();
    }

    iz->pop_back();

  }
  return value;
}

Kabbala Single_Amplitude_Base::SingleZvalue(Zfunc* z,vector<int>* iz, vector<int>* iargs,int last)
{
  if(last && z->GetSize()>1) return SingleZGroupvalue(z,iz,iargs,last);
   
  int* args = new int[2*z->m_narg];
  int k=FillArgs(z,args,iz,iargs);
  if(k!=-1&&z->GetSize()==1){
    delete[] args;
    return SingleZvalueTensor(z,iz,iargs,k); //produce tensors for extern spin2 bosons
  }

  Zfunc* ze = z->p_equal;
  //if (z->p_equal!=z) ze = z->p_equal;

  for (CL_Iterator clit=ze->m_calclist.begin();clit!=ze->m_calclist.end();++clit) {
     if ((*clit).Compare(args,2*z->m_narg)) {
      delete[] args;
      if (z->m_sign==-1) return -(*clit).m_value;
                    else return  (*clit).m_value;
     }
  }

  int* newargs = new int[2*z->m_narg];
  for (short int i=0;i<2*z->m_narg;i++) newargs[i] = args[i];
  Kabbala value;

  if(z->GetSize()>1)value=SingleZGroupvalue(z,iz,iargs,last);
  else {
    z->p_calculator->SetArgCouplProp(2*z->m_narg,args,z->p_couplings,z->m_nprop,z->p_propagators,&plist);

    value = z->p_calculator->Do();
  }
  
  if (buildstring) {
    if ( (value.String()).find(string("+"))!=string::npos ||
	 (value.String()).find(string("-"))!=string::npos ||
	 (value.String()).find(string("*"))!=string::npos )
      value = (shand->Get_Generator())->GetCZnumber(value.Value(),value.String());
  }
  else value.SetString("");

  ze->m_calclist.push_back(CValue(newargs,value));

  delete[] args;
  if (z->m_sign==-1) return -value;

  return value;
}

void Single_Amplitude_Base::SetLoopVar(vector<int>& iz,vector<vector<int> >& iargs)
{
  for (size_t i=0;i<iz.size();i++) if(iz[i]>99){
    if (iz[i]<199) {
#ifdef Cut_Fermion_Prop
      Pfunc* p = 0;
      for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
	p = *pit;
	if (p->arg[0]==iz[i]) break;	
      }
      double mass = 0.;
      for (short j=1;j<p->argnum;j++) {
	if (fl[p->arg[j]].IsAnti()) mass -= b[p->arg[j]] * fl[p->arg[j]].Mass();
	                       else mass += b[p->arg[j]] * fl[p->arg[j]].Mass();
      }

      double particlemass = p->fl.Mass();

      if (p->fl.IsAnti()) particlemass = -particlemass;

      if (ATOOLS::IsEqual(particlemass,mass) && (p->fl).Width()==0.) {
	for (short j=1;j<p->argnum;j++)	iargs[2*i].push_back(p->arg[j]);
	//mark this special propagator as negative
	iz[i] = -iz[i];
      } 
      else {
	//fermions
	iargs[2*i].push_back(-1);
	iargs[2*i].push_back(1);
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
}

Flavour* Single_Amplitude_Base::GetPflav(int pn)
{
  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc*  p = *pit;
    if(pn==p->arg[0])return &(p->fl);
  }

  MyStrStream stream;
  stream << "ERROR in Single_Amplitude_Base::GetPflav:\n";
  stream << "  Propagator " << pn << " not found. Abort the run.";
  THROW(fatal_error, stream.str());
}

Kabbala Single_Amplitude_Base::GetProp(Zfunc* z)
{
  Kabbala Pols(string("1"),Complex(1.,0.));
  Basic_Pfunc bpf(shand->Get_Generator(),BS);
  int cnt=0,sign=1;
  for(int i=0;i<z->m_nprop;i++){

    for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
      Pfunc*  p = *pit;
      int pc;

      if (!(z->p_propagators[i].maped)) pc = p->arg[0];
      else {
	pc = p->momnum;
	if((p->fl).Kfcode()!=z->p_propagators[i].kfcode) pc = -1;
      }
 
      if (pc==iabs(z->p_propagators[i].numb) && p->fl.Kfcode()!=kf_shgluon) {
	if(z->GetOp()==0 && p->on==1)break;
	if (p->haspol && p->fl.IsVector()) sign*=-1;
	Pols*= bpf.P(p);
	cnt++;
	break;
      } 
    }
  }
  if(sign<0)Pols=-Pols;

  if (buildstring) {
    if (cnt>1)
      Pols = (shand->Get_Generator())->GetCZnumber(Pols.Value(),Pols.String());
  }
  return Pols;
}


Kabbala Single_Amplitude_Base::SingleMassTerms(int iz,int iarg)
{  
  Kabbala factor(string("1"),Complex(1.,0.));
  if (iabs(iz)>=199)return factor;
  if (iz<0) if (b[iarg]<0) {
    factor= -factor;
    return factor;
  }    

  Basic_MassTermfunc bmtf(shand->Get_Generator(),BS);

  bmtf.SetArgCouplProp(0,0,0,0,0,&plist); 
  if(iz>0) {
    factor = bmtf.MassTerm(iz*iarg);
    factor*= Kabbala(string("0.5"),Complex(1./2.,0.));
  }
  return factor;
}


void Single_Amplitude_Base::GroupZfuncs() 
{
  vector<int> iz;

  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc* p = *pit;
    if (p->on) {
      if (p->arg[0]>99) iz.push_back(p->arg[0]);
    }
  }
  if (iz.empty()) return;
  vector<int> dummy;
  vector<vector<int> > iargs;iargs.reserve(2*iz.size());iargs.resize(2*iz.size(),dummy);
  
  SetLoopVar(iz,iargs);
  
  vector<vector<int> > indexlist;
  int cnt=0;
  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
    Zfunc* z = (*zit);
    indexlist.push_back(dummy);
    for (short int i=0;i<z->m_narg;i++) {
      if (z->p_arguments[i]>99) {
	for (size_t j=0;j<iz.size();j++) {
	  if (iabs(iz[j])==z->p_arguments[i]) {
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
    int imin=0;
    for(int i=0;i<(int)iz.size();i++){
      int adds=0;
      for(cnt=0;cnt<(int)indexlist.size();cnt++)
	for(size_t j=0;j<indexlist[cnt].size();j++) {
	  over=0;
	  if(indexlist[cnt][j]==i){
	    if(adds==0)adds=1;
	    for(size_t k=0;k<indexlist[cnt].size();k++)
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
      zh[0]=0;zh[1]=0;
      indexlist.push_back(dummy);
      
      vector<vector<int> >::iterator ilt=indexlist.begin();	
      for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();) {
	int hit=0;
	for(size_t j=0;j<(*ilt).size();j++) {
	  if((*ilt)[j]==imin){
	    zh[ia]=(*zit);
	    ia++;
	    for(size_t k=0;k<(*ilt).size();k++)
	      if((*ilt)[k]!=imin)
		{
		  indexlist[indexlist.size()-1].push_back((*ilt)[k]);
		}
	    ilt=indexlist.erase(ilt);
	    zit=zlist->erase(zit);
	    hit=1;
	    break;
	  }
	}
	if(!hit){++ilt;++zit;}
      }
      if(ia!=2){
	msg_Error()<<"ERROR Single_Amplitude_Base::GroupZfuncs: "<<std::endl
		   <<"   Index appeared "<<ia<<" times, will abort the run."<<endl;
	zh[0]->Print();zh[1]->Print();
	Abort();
      }
      
      Zfunc_Group *sf=new Zfunc_Group(*zh[0],*zh[1],iabs(iz[imin]),&plist);
      zlist->push_back(sf);
      
    }
  } while(over==0);
}

Complex Single_Amplitude_Base::Zvalue(String_Handler * sh, int ihel) 
{
  if (sh) return sh->Zvalue(amplnumber,ihel);
  msg_Error()<<"ERROR in Single_Amplitude_Base::Zvalue(String_Handler * sh, int ihel) : "<<endl
	     <<"   Will try to circumvent this and continue the run."<<std::endl; 
  shand->Zvalue(amplnumber,ihel);   
  return Complex(0.,0.);
}

Complex Single_Amplitude_Base::Zvalue(int ihel) 
{
  return shand->Zvalue(amplnumber,ihel);
}

Complex Single_Amplitude_Base::Zvalue(int ihel,int* signlist) 
{
  if (signlist==0) return shand->Zvalue(amplnumber,ihel);
  
  if(zlist->size()>1) GroupZfuncs();
    
  Kabbala value(string("1"),Complex(1.,0.));
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
  for (Zfunc_Iterator zit=zlist->begin();zit!=zlist->end();++zit) {
    Zfunc* z = (*zit);

    value*= SingleZvalue(z,&ize,&iarg,1);
    if (z->m_nprop>0) value*= GetProp(z);
  }

  //   scalar propagators:
#ifndef Scalar_Args
  Basic_Pfunc bpf(shand->Get_Generator(),BS);

  for (Pfunc_Iterator pit=plist.begin();pit!=plist.end();++pit) {
    Pfunc*  p = *pit;
    if ((p->fl).IsScalar()) {
      value *= bpf.P(p);
    } 
  }
#endif

  if (sign<0 && value.Value()!=Complex(0.,0.)) value = -value; 
  if (buildstring) shand->Set_String(amplnumber,ihel,value.String());
  return value.Value();
}

void Single_Amplitude_Base::DefineOrder(const std::vector<int> &o)
{
  m_order=o;
}

const std::vector<int> &Single_Amplitude_Base::GetOrder()
{
  return m_order;
}
