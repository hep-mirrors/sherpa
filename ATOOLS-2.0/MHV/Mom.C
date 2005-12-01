#include "Mom.H"
#include "Color.H"
#include "Exception.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"
#include <fstream.h>
#include <iostream.h>
#include "MathTools.H"


using namespace ATOOLS;
using namespace AMEGIC;

// class Mom

// constructor
Mom::Mom(Vec4D& m,int h,int p):hel(h),part(p) {
  for (size_t i=0;i<4;i++) (*this)[i]=m[i];
  if (!IsZero(Abs2()))  THROW(fatal_error,"Massive particles are not implemented yet");
} 

Mom::Mom(Vec4D& m):hel(0),part(0) {
  for (size_t i=0;i<4;i++) (*this)[i]=m[i];
} 


inline int Mom::GetHelicity() {
  return hel;
}

inline int Mom::GetPartType() {
  return part;
}

void Mom::Print() {
  msg_Info()<<(*this)<<"    p^2=";
  if (IsZero(Abs2())) {msg_Info()<<0;}
  else msg_Info()<<Abs2(); 
  msg_Info()<<"    helicity="<<hel<<"    particle ("<<part<<")"<<endl;
  }


// class MomentumList


// constructor
MomentumList::MomentumList(const char*  file) {
  etha[0]=0; etha[1]=0;
  if (!Get(file))  THROW(fatal_error,"Error in initializing class MomentumList");
  MakeHlist(); 
  MakePlist();
  Vec4D m;
  Mom* momentum = new Mom(m) ; 
  this->push_back(momentum);
  if (!etha[0] && !etha[1]) {
    etha[0]=1; 
    etha[1]=1;
  }
  cout<<"etha = ("<<etha[0]<<","<<etha[1]<<")"<<endl; //print
} 

//destructor
MomentumList::~MomentumList() {
  for (size_t i=0;i<size();i++) delete (*this)[i];
  delete hlist;
  delete plist;
}

bool MomentumList::Get(const char* file) {
  std::string  line;
  ifstream dat(file);
  if (dat==NULL)
    THROW(fatal_error,"Error in opening a data file");
  while (getline(dat,line)) {
    if (line.find("%")!=std::string::npos) line=line.substr(0,line.find("%"));
    if (line.find("p_[")!=std::string::npos) {
      size_t stpos(line.find("p_[")+2);
      size_t endpos(line.find("]"));
      if (stpos>endpos)
	THROW(fatal_error,"Momentum is not properly defined");
      size_t c1pos(line.find(',',stpos+1));
      if (c1pos==std::string::npos)
	THROW(fatal_error,"Invalid number of momentum indices");
      size_t c2pos(line.find(',',c1pos+1));
      if (c2pos==std::string::npos)
	THROW(fatal_error,"Invalid number of momentum indices");
      size_t c3pos(line.find(',',c2pos+1));
      if (c3pos==std::string::npos || 
	  line.find(',',c3pos+1)<endpos || c3pos>endpos)
	THROW(fatal_error,"Invalid number of momentum indices");
      Vec4D m;
      m[0]= ToType<double>(line.substr(stpos+1,c1pos-stpos-1));
      m[1]= ToType<double>(line.substr(c1pos+1,c2pos-c1pos-1));
      m[2]= ToType<double>(line.substr(c2pos+1,c3pos-c2pos-1));
      m[3]= ToType<double>(line.substr(c3pos+1,endpos-c3pos-1));
      if (line.substr(endpos+1,2)=="_t") m[0]=m.PSpat();
      if (line.substr(endpos+1,2)=="_s") {
	double sp(m.PSpat());
	for (int i=1;i<4;i++) m[i]=m[i]*abs(m[0])/sp;
      }
      if (line.find("h")==std::string::npos)	
	THROW(fatal_error,"No helicity specified");
      int h=0;
      if (line.substr(line.find("h")+1,1)=="+") h=1;
      if (line.substr(line.find("h")+1,1)=="-") h=-1;
      if (!h) THROW(fatal_error,"No helicity specified");
      stpos=line.find("(");
      endpos=line.find(")");
      if (stpos+2>endpos)
	THROW(fatal_error,"Particle type is not properly defined");
      int p=ToType<int>(line.substr(stpos+1,endpos-stpos-1));
      Mom* momentum = new Mom(m,h,p) ; 
      this->push_back(momentum);
    } 
    if (line.find("etha=[")!=std::string::npos) {
      size_t stpos(line.find("etha=[")+5);
      size_t endpos(line.find("]"));
      if (stpos>endpos || endpos==std::string::npos)
	THROW(fatal_error,"etha is not properly defined");
      size_t c1pos(line.find(',',stpos+1));
      if (c1pos==std::string::npos || line.find(',',c1pos+1)<endpos || c1pos>endpos)
	THROW(fatal_error,"Definition of etha requires exectly two numbers");
      etha[0]= ToType<double>(line.substr(stpos+1,c1pos-stpos-1));
      etha[1]= ToType<double>(line.substr(c1pos+1,endpos-c1pos-1));
    }
  }
  dat.close();
  m_size=size();
  if (size()>0) return 1;
  return 0;
}

void MomentumList::MakeHlist() {
  hlist = new int[Size()];
  for (size_t i=0;i<Size();i++) hlist[i]=(*this)[i]->GetHelicity();
}

void MomentumList::MakePlist() {
  plist = new int[Size()];
  int aqpos1=-1, aqtype1=0;
  int aqpos2=-1, aqtype2=0;
  for (size_t i=0;i<Size();i++) {
    plist[i]=(*this)[i]->GetPartType();
    if (plist[i]<0 && plist[i]>-10) {
      if (aqpos1==-1) {aqpos1=i; aqtype1=plist[i];}
      else if (aqpos2==-1) {aqpos2=i; aqtype2=plist[i];}
      else THROW(fatal_error,"To many quark-antiquark pairs");
    }
  }
  if (aqpos1>-1 && plist[(aqpos1+1)%Size()]!=-aqtype1) THROW(fatal_error,"Wrong quark-antiquark structure"); 
  if (aqpos2>-1 && plist[(aqpos2+1)%Size()]!=-aqtype2) THROW(fatal_error,"Wrong quark-antiquark structure");
}

int* MomentumList::GetHList() {
  return hlist;
}

int* MomentumList::GetPList() {
  return plist;
}

Complex MomentumList::S0(int i, int j) {
  Complex li[2],lj[2];
  if (i<Size()) {
    if ((*(*this)[i])[0]+(*(*this)[i])[3]) {
      li[0]=csqrt((*(*this)[i])[0]+(*(*this)[i])[3]);
      li[1]=Complex((*(*this)[i])[1],(*(*this)[i])[2])/csqrt((*(*this)[i])[0]+(*(*this)[i])[3]);
    }
    else {
      li[0]=0;
      li[1]=csqrt((*(*this)[i])[0]-(*(*this)[i])[3]);
    }
  }
  else {
    li[0]=((*(*this)[i])[0]+(*(*this)[i])[3])*etha[0];
    li[0]+=Complex((*(*this)[i])[1],-(*(*this)[i])[2])*etha[1];
    li[1]=Complex((*(*this)[i])[1],(*(*this)[i])[2])*etha[0];
    li[1]+=((*(*this)[i])[0]-(*(*this)[i])[3])*etha[1];
  }

  if (j<Size()) {
    if ((*(*this)[j])[0]+(*(*this)[j])[3]) {
      lj[0]=csqrt((*(*this)[j])[0]+(*(*this)[j])[3]);
      lj[1]=Complex((*(*this)[j])[1],(*(*this)[j])[2])/csqrt((*(*this)[j])[0]+(*(*this)[j])[3]);
    }
    else {
      lj[0]=0;
      lj[1]=csqrt((*(*this)[j])[0]-(*(*this)[j])[3]);
    }
  }
  else {
    lj[0]=((*(*this)[j])[0]+(*(*this)[j])[3])*etha[0];
    lj[0]+=Complex((*(*this)[j])[1],-(*(*this)[j])[2])*etha[1];
    lj[1]=Complex((*(*this)[j])[1],(*(*this)[j])[2])*etha[0];
    lj[1]+=((*(*this)[j])[0]-(*(*this)[j])[3])*etha[1];
  }
  return li[0]*lj[1]-li[1]*lj[0];
}


Complex MomentumList::S1(int i, int j) { 
  Complex li[2],lj[2];
  if (i<Size()) {
    if ((*(*this)[i])[0]+(*(*this)[i])[3]) {
      li[0]=csqrt((*(*this)[i])[0]+(*(*this)[i])[3]);
      li[1]=Complex((*(*this)[i])[1],-(*(*this)[i])[2])/csqrt((*(*this)[i])[0]+(*(*this)[i])[3]);
    }
    else {
      li[0]=0;
      li[1]=csqrt((*(*this)[i])[0]-(*(*this)[i])[3]);
    }
  }
  else {
    li[0]=((*(*this)[i])[0]+(*(*this)[i])[3])*etha[0];
    li[0]+=Complex((*(*this)[i])[1],(*(*this)[i])[2])*etha[1];
    li[1]=Complex((*(*this)[i])[1],-(*(*this)[i])[2])*etha[0];
    li[1]+=((*(*this)[i])[0]-(*(*this)[i])[3])*etha[1];
  }

  if (j<Size()) {
    if ((*(*this)[j])[0]+(*(*this)[j])[3]) {
      lj[0]=csqrt((*(*this)[j])[0]+(*(*this)[j])[3]);
      lj[1]=Complex((*(*this)[j])[1],-(*(*this)[j])[2])/csqrt((*(*this)[j])[0]+(*(*this)[j])[3]);
    }
    else {
      lj[0]=0;
      lj[1]=csqrt((*(*this)[j])[0]-(*(*this)[j])[3]);
    }
  }
  else {
    lj[0]=((*(*this)[j])[0]+(*(*this)[j])[3])*etha[0];
    lj[0]+=Complex((*(*this)[j])[1],(*(*this)[j])[2])*etha[1];
    lj[1]=Complex((*(*this)[j])[1],-(*(*this)[j])[2])*etha[0];
    lj[1]+=((*(*this)[j])[0]-(*(*this)[j])[3])*etha[1];
  }
  return li[0]*lj[1]-li[1]*lj[0];
}


int MomentumList::GetMomNumber(Pfunc* pf) {
  Vec4D m;
  for (size_t i=0;i<4;i++) {
    m[i]=0;
    for (int j=1;j<pf->argnum;j++) m[i]+=(*(*this)[pf->arg[j]])[i];
  } 
  Mom* momentum = new Mom(m) ; 
  this->back()=momentum;
  return size()-1;
}

Vec4D MomentumList::Momentum(size_t mindex) {
  Vec4D m;
  for (size_t i=0;i<4;i++) m[i]=(*(*this)[mindex])[i];
  return m;
}





void MomentumList::Print() {
  for (size_t i=0;i<Size();i++) (*this)[i]->Print();
  msg_Info()<<endl<<"helicity list: ("<<hlist[0];
  for (size_t i=1;i<Size();i++)  msg_Info()<<","<<hlist[i]; 
  msg_Info()<<")"<<endl; 
  msg_Info()<<"particle list: ("<<plist[0];
  for (size_t i=1;i<Size();i++)  msg_Info()<<","<<plist[i]; 
  msg_Info()<<")"<<endl<<endl; 
  //for (size_t i=0;i<Size();i++) for (size_t j=i+1;j<Size();j++) { 
  //  msg_Info()<<"<"<<i<<","<<j<<"> = "<<S0(i,j); //print
  //  msg_Info()<<"     ["<<i<<","<<j<<"] = "<<S1(i,j)<<endl; //print
  //}
}



// class Fullamplitude_MHV_test

// constructor
Fullamplitude_MHV_test::Fullamplitude_MHV_test(MomentumList *momentuml,int *hl,int *pl): momentumlist(momentuml), hlist(hl), plist(pl) {} 

//destructor
//Fullamplitude_MHV_test::~Fullamplitude_MHV_test() {}


Complex Fullamplitude_MHV_test::Calculate(int pr) {
  print=pr;
  int qpos=1;
  amp=Complex(0.0,0.0);
  int part(momentumlist->Size());
  calc = new MHV_Calculator(part,momentumlist,hlist,plist,print);
  int* qlist(calc->GetQlist());
  perm = new int[part];
  int** perm_adr  = new int*[part-qlist[0]+1];
  if (!qlist[0]) {
    perm[part-1] = part-1;
    for (int i=0;i<part-1;i++) perm_adr[i]=&perm[i];
    Permutation_pureg(part-2,perm_adr);
    amp*=pow(3,part-2);
    amp*=8; amp/=pow(2,part);
  }
  else if (qlist[0]==2 || qlist[0]==4) {
    ma = new size_t[part-qlist[0]+1];
    ma[0] = part-qlist[0];
    for (int i=0;i<part;i++) {
      if (i==qlist[qpos]) {
	qpos++;
	perm[i]=i;
      }
      else { 
	ma[i+2-qpos]=i+1;
	perm_adr[i+1-qpos]=&perm[i];
      }
    }
    Permutation_quark(part-qlist[0]-1,perm_adr);
    amp/=pow(3,part-qlist[0]-1);
    amp*=8; amp/=pow(2,part-qlist[0]);
    delete [] ma;
  }
 
  delete [] perm_adr; 
  delete [] perm;
  delete calc;
  return amp; 
}

void Fullamplitude_MHV_test::Permutation_pureg(int p_number,int** p_adr) {
  if (p_number) {    
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number; 
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      Permutation_pureg(p_number-1,perm_adr);
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=0;
    Complex am=calc->Differential(perm);
    amp+=norm(am);
    
    // start 6 gluons TEST    
    if (momentumlist->Size()==6) {
      am=conj(am);
      int perm1[]={0,2,4,1,5,3};
      int permz[6];
      for (int l=0;l<6;l++) permz[l]=perm[perm1[l]];
      Complex am2=calc->Differential(permz);
      int perm2[]={0,4,2,5,1,3};
      for (int l=0;l<6;l++) permz[l]=perm[perm2[l]];
      am2 +=calc->Differential(permz);
      int perm3[]={0,2,5,3,1,4};
      for (int l=0;l<6;l++) permz[l]=perm[perm3[l]];
      am2+=calc->Differential(permz);
      am*=am2;
      am*=2;
      am/=pow(3,2);
      amp+=am;
    }
    // end 6 gluons TEST

    if (print) {                                                              //print
      msg_Info()<<"     perm: ("<<perm[0];                                    //print
      for (size_t i=1;i<momentumlist->Size();i++)  msg_Info()<<","<<perm[i];  //print 
      msg_Info()<<")"  ;                                                      //print
    }                                                                         //print
    if (print) msg_Info()<<"     "<<am<<endl;                                 //print 
    return;
  }
}
void Fullamplitude_MHV_test::Permutation_quark(int p_number,int** p_adr) {
   if (p_number>0) {
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=ma[p_number+1]-1;
      int** perm_adr = new int*[p_number]; 
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)]; 
      Permutation_quark(p_number-1,perm_adr); 
      delete [] perm_adr; 
    }
  }
   else {
    if (ma[0]) (*p_adr[0])=ma[1]-1;
    Complex am=calc->Differential(perm); 
    if (print) msg_Info()<<"     "<<am<<endl;                                 //print 
    
    // start 2 gluons + 2 quarks TEST    
    if ((momentumlist->Size()-ma[0])==2) {
      if (momentumlist->Size()==4) { 
	amp+= 9*norm(am);

	am=-conj(am);
	int perm1[]={0,1,2,3};
	int permz[4];
	for (int l=0;l<4;l++) permz[l]=perm[perm1[l]];
	Complex am2=calc->Differential(permz);
	int perm2[]={0,2,1,3};
	for (int l=0;l<4;l++) permz[l]=perm[perm2[l]];
	am2 +=calc->Differential(permz);
	am*=am2;
	amp+=am; 
      }
      // end 2 gluons + 2 quarks TEST

      // start 3 gluons + 2 quarks TEST 
      else if (momentumlist->Size()==5) {
	amp+= 81*norm(am);
	
	Complex am11=-conj(am);
	int perm1[]={0,1,2,3,4};
	int permz[5];
	for (int l=0;l<5;l++) permz[l]=perm[perm1[l]];
	Complex am12=calc->Differential(permz);
	am12*=2;
	int perm2[]={4,1,2,3,0};
	for (int l=0;l<5;l++) permz[l]=perm[perm2[l]];
	am12 +=calc->Differential(permz);
	int perm3[]={0,1,2,4,3};
	for (int l=0;l<5;l++) permz[l]=perm[perm3[l]];
	am12 +=calc->Differential(permz);
	int perm4[]={3,1,2,0,4};
	for (int l=0;l<5;l++) permz[l]=perm[perm4[l]];
	am12 -=calc->Differential(permz);
	am11*=am12;
	am11*=9;
	amp+=am11;
	
	am=conj(am);
	int perm21[]={0,1,2,3,4};
	for (int l=0;l<5;l++) permz[l]=perm[perm21[l]];
	Complex am22= calc->Differential(permz);
	int perm22[]={0,1,2,4,3};
	for (int l=0;l<5;l++) permz[l]=perm[perm22[l]];
	am22 +=calc->Differential(permz);
	int perm23[]={3,1,2,0,4};
	for (int l=0;l<5;l++) permz[l]=perm[perm23[l]];
	am22 +=calc->Differential(permz);
	int perm24[]={4,1,2,0,3};
	for (int l=0;l<5;l++) permz[l]=perm[perm24[l]];
	am22 +=calc->Differential(permz); 
	int perm25[]={3,1,2,4,0};
	for (int l=0;l<5;l++) permz[l]=perm[perm25[l]];
	am22 +=calc->Differential(permz);
	int perm26[]={4,1,2,3,0};
	for (int l=0;l<5;l++) permz[l]=perm[perm26[l]];
	am22 +=calc->Differential(permz);
	am*=am22;
	amp+=am;
      }
      // end 3 gluons + 2 quarks TEST
    }
    else amp+=norm(am);


    if (print) {                                                              //print
      msg_Info()<<"     perm: ("<<perm[0];                                    //print
      for (size_t i=1;i<momentumlist->Size();i++)  msg_Info()<<","<<perm[i];  //print 
      msg_Info()<<")"<<endl  ;                                                      //print
    }                                                                         //print
    return;
  }
}





// class Fullamplitude_MHV

// constructor
Fullamplitude_MHV::Fullamplitude_MHV(int np,int *pl,int pr): 
  n_part(np), ampsq(0.,0.), plist(pl), print(pr),colorstore(false) 
{
  perm = new int[n_part];
} 

//destructor
Fullamplitude_MHV::~Fullamplitude_MHV() 
{
  delete [] perm;
  if (colorstore) delete permstore;  
}




// Pure gluonic amplitude
void Fullamplitude_MHV::PermutationStoreColor_pureg(int p_number,int** p_adr) {
  if (p_number) {  
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number; 
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      PermutationStoreColor_pureg(p_number-1,perm_adr);
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=0;
    
    size_t *mact = new size_t[n_part+1];
    mact[0]=n_part;                
    size_t* maconj = new size_t[n_part+1];
    maconj[0]=n_part;
    for (int i=1;i<n_part+1;i++) {
      mact[i]=perm[i-1]+1;
      maconj[i]=n_part+1-i;
    } 
    Expression expression(n_part,0);
    expression[0] = Trace::New(mact,0,0);
    expression.push_back(Trace::New(maconj,0,0));
    expression.Evaluate();
    Complex col=expression.Result();
    size_t *perms = new size_t[n_part-1];
    for (int i=0;i<n_part-1;i++) perms[i]=perm[i];
    permstore->PutColor(perms,col);
  
    delete [] perms;
    delete [] maconj;
    delete [] mact;
    return;
  }
}


void Fullamplitude_MHV::PermutationStoreAmp_pureg(int p_number,int** p_adr) {
  if (p_number) {  
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number; 
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      PermutationStoreAmp_pureg(p_number-1,perm_adr);
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=0;
          
    Complex amp(calc->Differential(perm));
    size_t *perms = new size_t[n_part-1];
    for (int i=0;i<n_part-1;i++) perms[i]=perm[i];
    permstore->PutAmp(perms,amp);
  
    delete [] perms;
    return;
  }
}


void Fullamplitude_MHV::Permutation_pureg(int p_number,int** p_adr,bool conjflag) {
  if (p_number) {    
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number; 
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      Permutation_pureg(p_number-1,perm_adr,conjflag);
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=0;

    if (!conjflag) { 
      size_t *perms = new size_t[n_part-1];
      for (int i=0;i<n_part-1;i++) perms[i]=perm[i];
      ampnc = permstore->GetAmp(perms);

      int** perm_adr = new int*[n_part-1];
      for (int i=0;i<n_part-1;i++) perm_adr[i]=&permconj[i];
      Permutation_pureg(n_part-2,perm_adr,true);
      delete [] perm_adr;
    }
    else { 
      Complex amp(ampnc);
      size_t *perms = new size_t[n_part-1];
      for (int i=0;i<n_part-1;i++) perms[i]=permconj[i];
      amp*=conj(permstore->GetAmp(perms));

      size_t *permi = new size_t[n_part-1];
      for (int i=0;i<n_part-1;i++) permi[permconj[i]]=i;
      for (int i=0;i<n_part-1;i++) perms[i]=permi[perm[i]];
      amp*=permstore->GetColor(perms); 
     
      ampsq+=amp;
      delete [] permi;
      delete [] perms;
      return;
    } 
  }
}




// 2 quarks amplitude
void Fullamplitude_MHV::PermutationStoreColor_quark2(int p_number,std::vector<int>** p_adr) {
  if (p_number) {
    for (int l=0;l<p_number+1;l++) {
      (*p_adr[l])[0]=permgl[p_number];
      (*p_adr[l])[1]=p_number;
      std::vector<int>** perm_adr = new std::vector<int>*[p_number]; 
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)]; 
      PermutationStoreColor_quark2(p_number-1,perm_adr); 
      delete [] perm_adr; 
    }
  }
  else { 
    (*p_adr[0])[0]=permgl[0];
    (*p_adr[0])[1]=0;
   
    size_t *mact = new size_t[n_part-1];
    mact[0]=n_part-2;                
    size_t* maconj = new size_t[n_part-1];
    maconj[0]=n_part-2;
    size_t *perms = new size_t[n_part-2];
    for (int i=0;i<n_part-2;i++) {
      perms[i]=(permtb[i])[1];
      mact[i+1]=(permtb[i])[0]+1;
      maconj[n_part-2-i]=permgl[i]+1;
    }   

    Expression expression(n_part,2);
    expression[0] = Trace::New(mact,1,2);
    expression.push_back(Trace::New(maconj,2,1));
    expression.Evaluate();
    Complex col=expression.Result(); 
    permstore->PutColor(perms,col);
   
    delete [] perms;
    delete [] maconj;
    delete [] mact;
    return;
  }
}


void Fullamplitude_MHV::PermutationStoreAmp_quark2(int p_number,std::vector<int>** p_adr) {
  if (p_number) {
    for (int l=0;l<p_number+1;l++) {
      (*p_adr[l])[0]=permgl[p_number];
      (*p_adr[l])[1]=p_number;
      std::vector<int>** perm_adr = new std::vector<int>*[p_number]; 
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)]; 
      PermutationStoreAmp_quark2(p_number-1,perm_adr); 
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])[0]=permgl[0];
    (*p_adr[0])[1]=0;

    size_t *perms = new size_t[n_part-2];
    for (int i=0;i<n_part-2;i++) perms[i]=(permtb[i])[1];
    for (int i=0;i<n_part;i++) perm[i]=(permtb[i])[0];

    Complex amp(calc->Differential(perm));
    permstore->PutAmp(perms,amp);    

    delete [] perms;
    return;
  }
}


void Fullamplitude_MHV::Permutation_quark2(int p_number,int** p_adr,bool conjflag) {
  if (p_number) {    
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number; 
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      Permutation_quark2(p_number-1,perm_adr,conjflag);
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=0;

    if (!conjflag) { 
      size_t *perms = new size_t[n_part-2];
      for (int i=0;i<n_part-2;i++) perms[i]=perm[i];
      ampnc = permstore->GetAmp(perms)
;
      int** perm_adr = new int*[n_part-2];
      for (int i=0;i<n_part-2;i++) perm_adr[i]=&permconj[i];
      Permutation_quark2(n_part-3,perm_adr,true);
      delete [] perm_adr;
    }
    else {
      Complex amp(ampnc); 
      size_t *perms = new size_t[n_part-2];
      for (int i=0;i<n_part-2;i++) perms[i]=permconj[i];
      amp*=conj(permstore->GetAmp(perms));

      size_t *permi = new size_t[n_part-2];
      for (int i=0;i<n_part-2;i++) permi[permconj[i]]=i;
      for (int i=0;i<n_part-2;i++) perms[i]=permi[perm[i]];
      amp*=permstore->GetColor(perms); 
     
      ampsq+=amp;
      delete [] permi;
      delete [] perms;
      return;
    } 
  }
}






// 4 quarks amplitude
void Fullamplitude_MHV::Permutation_quark4(int p_number, int** p_adr,bool conjflag) {
  if (p_number>0) {
    for (int l=0;l<p_number+1;l++) {
      if (conjflag) *p_adr[l]=permconjgl[p_number];
      else *p_adr[l]=permgl[p_number];
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      Permutation_quark4(p_number-1,perm_adr,conjflag); 
      delete [] perm_adr; 
    }
  }
  else { 
    if (!conjflag) {
      (*p_adr[0])=permgl[0]; 
      
      int** perm_adr = new int*[n_part-4];
      int* plist(calc->GetPlist());
      int fq=0;
      if (plist[qlist[1]]>0) fq=1;
      ampnc = calc->Differential(perm);

      msg_Info()<<endl<<"     amp="<<ampnc<<endl<<endl;                       //print
   

      for (m_qposconj=1;m_qposconj<n_part-2;m_qposconj++) {
	permconj[m_qposconj-1]=qlist[3-fq];
	permconj[m_qposconj]=qlist[4-fq];

	int qpos=1;
	for (int i=0;i<n_part-2;i++) {
	  if (i==m_qposconj || i==(m_qposconj-1)) qpos++;
	  else perm_adr[i+1-qpos]=&permconj[i];
	}
	qpos=1;
	for (int i=0;i<n_part;i++) {
	  if (qlist[qpos]==i) qpos++;
	  else permconjgl[i+1-qpos]=i;	
	}
	Permutation_quark4(n_part-5,perm_adr,true); 
      }  
      
      int st=perm[m_qpos-1];
      perm[m_qpos-1]=perm[n_part-2];
      perm[n_part-2]=st;
      ampnc = calc->Differential(perm); msg_Info()<<endl<<"     amp="<<ampnc<<endl<<endl;                       //print
      ampnc/= (-NC);
      for (m_qposconj=1;m_qposconj<n_part-2;m_qposconj++) {
	permconj[m_qposconj-1]=qlist[3-fq];
	permconj[m_qposconj]=qlist[4-fq];

	int qpos=1;
	for (int i=0;i<n_part-2;i++) {
	  if (i==m_qposconj || i==(m_qposconj-1)) qpos++;
	  else perm_adr[i+1-qpos]=&permconj[i];
	}
	qpos=1;
	for (int i=0;i<n_part;i++) {
	  if (qlist[qpos]==i) qpos++;
	  else permconjgl[i+1-qpos]=i;	
	}
	Permutation_quark4(n_part-5,perm_adr,true); 
      }  
      st=perm[m_qpos-1];
      perm[m_qpos-1]=perm[n_part-2];
      perm[n_part-2]=st;
      
      delete [] perm_adr;
    }
    else {
      (*p_adr[0])=permconjgl[0];

      Complex amp(ampnc);
      amp*=conj(calc->Differential(permconj));
      size_t *mact1 = new size_t[m_qpos]; 
      size_t *mact2 = new size_t[n_part-2-m_qpos];
      mact1[0]=m_qpos-1;
      mact2[0]=n_part-3-m_qpos;	
      for (size_t i=1;i<mact1[0]+1;i++) mact1[i]=perm[i-1]+1;
      for (size_t i=1;i<mact2[0]+1;i++) mact2[i]=perm[m_qpos+i]+1;
      Expression expression(n_part,4);
      Expression expression2(n_part,4);
      expression[0] = Trace::New(mact1,1,2); 
      expression.push_back(Trace::New(mact2,3,4));
      expression2[0] = Trace::New(mact1,1,2); 
      expression2.push_back(Trace::New(mact2,3,4));
 
      delete [] mact1;
      delete [] mact2; 

      mact1 = new size_t[m_qposconj]; 
      mact2 = new size_t[n_part-2-m_qposconj];
      mact1[0]=m_qposconj-1;
      mact2[0]=n_part-3-m_qposconj;	
      for (size_t i=1;i<mact1[0]+1;i++) mact1[mact1[0]+1-i]=permconj[i-1]+1;
      for (size_t i=1;i<mact2[0]+1;i++) mact2[mact2[0]+1-i]=permconj[m_qposconj+i]+1;
      expression.push_back(Trace::New(mact1,2,1)); 
      expression.push_back(Trace::New(mact2,4,3));
      expression.Evaluate(); 
      amp*=expression.Result(); 

      delete [] mact1;
      delete [] mact2;
      
      ampsq+=amp;

      amp=-ampnc/NC;
      int st=permconj[m_qposconj-1];
      permconj[m_qposconj-1]=permconj[n_part-2];
      permconj[n_part-2]=st;
      amp*=conj(calc->Differential(permconj));
     
      mact1 = new size_t[m_qposconj]; 
      mact2 = new size_t[n_part-2-m_qposconj];
      mact1[0]=m_qposconj-1;
      mact2[0]=n_part-3-m_qposconj;	
      for (size_t i=1;i<mact1[0]+1;i++) mact1[mact1[0]+1-i]=permconj[i-1]+1;
      for (size_t i=1;i<mact2[0]+1;i++) mact2[mact2[0]+1-i]=permconj[m_qposconj+i]+1;
      expression2.push_back(Trace::New(mact1,2,1)); 
      expression2.push_back(Trace::New(mact2,4,3));
      expression2.Evaluate(); 
      amp*=expression2.Result(); 

      st=permconj[m_qposconj-1];
      permconj[m_qposconj-1]=permconj[n_part-2];
      permconj[n_part-2]=st;
      delete [] mact1;
      delete [] mact2;
     
      ampsq+=amp;
 
      return;
    }
  }
}




// Full amplitude 
Complex Fullamplitude_MHV::MSquare(MomentumList *p_BS,int *hlist)
{
  calc = new MHV_Calculator(n_part,p_BS,hlist,plist,print); 
  qlist=(calc->GetQlist());
 
  // pure gluons
  if (qlist[0]==0) {
    permconj = new int[n_part];
    int** perm_adr = new int*[n_part-1];

    perm[n_part-1] = n_part-1;
    permconj[n_part-1] = n_part-1;
    for (int i=0;i<n_part-1;i++) perm_adr[i]=&perm[i];
    if (!colorstore){
      permstore = new PermStore(n_part-1);
      PermutationStoreColor_pureg(n_part-2,perm_adr);
      colorstore=true;
    }
    PermutationStoreAmp_pureg(n_part-2,perm_adr);
    Permutation_pureg(n_part-2,perm_adr,false);

    delete [] perm_adr; 
    delete [] permconj;
  }

  // 2 quarks
  else if (qlist[0]==2) {  
    permgl = new int[n_part-2];
    permtb = new std::vector<int>[n_part];
    std::vector<int>** permtb_adr = new std::vector<int>*[n_part-2];

    if (qlist[1]==0 && qlist[2]==n_part-1) {
      (permtb[n_part-2]).push_back(qlist[2]);
      (permtb[n_part-2]).push_back(-1);
      (permtb[n_part-1]).push_back(qlist[1]);
      (permtb[n_part-1]).push_back(-1);
    }
    else {
      (permtb[n_part-2]).push_back(qlist[1]);
      (permtb[n_part-2]).push_back(-1);
      (permtb[n_part-1]).push_back(qlist[2]);
      (permtb[n_part-1]).push_back(-1);
    }
    int qpos=1;
    for (int i=0;i<n_part;i++) {
      if (qlist[qpos]==i) qpos++;
      else {
	permtb_adr[i+1-qpos]=&permtb[i+1-qpos];
	permgl[i+1-qpos]=i;	
	(permtb[i+1-qpos]).push_back(0);
	(permtb[i+1-qpos]).push_back(0);
      }
    }
    if (!colorstore){
      permstore = new PermStore(n_part-2);
      PermutationStoreColor_quark2(n_part-3,permtb_adr);
      colorstore=true;
    }
    PermutationStoreAmp_quark2(n_part-3,permtb_adr);

    delete [] permtb_adr;
    delete [] permtb;
    delete [] permgl;

    permconj = new int[n_part-2];
    int** perm_adr = new int*[n_part-2];
    for (int i=0;i<n_part-2;i++) perm_adr[i]=&perm[i];
    Permutation_quark2(n_part-3,perm_adr,false); 

    delete [] perm_adr; 
    delete [] permconj;
  }
  
  // 4 quarks
  else if (qlist[0]==4) {
    int** perm_adr = new int*[n_part-4];
    permconj = new int[n_part];
    permgl = new int[n_part-4];
    permconjgl = new int[n_part-4];

    int* plist(calc->GetPlist());
    int fq=0;
    if (plist[qlist[1]]>0) {
      fq=1;
      perm[n_part-2]=qlist[4];
      perm[n_part-1]=qlist[1];
      permconj[n_part-2]=qlist[4];
      permconj[n_part-1]=qlist[1];
    }
    else {
      perm[n_part-2]=qlist[1];
      perm[n_part-1]=qlist[2];
      permconj[n_part-2]=qlist[1];
      permconj[n_part-1]=qlist[2];
    }
    
    for (m_qpos=1;m_qpos<n_part-2;m_qpos++) {
      perm[m_qpos-1]=qlist[3-fq];
      perm[m_qpos]=qlist[4-fq];

      int qpos=1;
      for (int i=0;i<n_part-2;i++) {
	if (i==m_qpos || i==(m_qpos-1)) qpos++;
	else perm_adr[i+1-qpos]=&perm[i];
      }
      qpos=1;
      for (int i=0;i<n_part;i++) {
	if (qlist[qpos]==i) qpos++;
	else permgl[i+1-qpos]=i;
      }
      Permutation_quark4(n_part-5,perm_adr,false);
    }

    delete [] permconjgl;
    delete [] permgl;
    delete [] permconj;
    delete [] perm_adr; 
  
  }

  delete calc;
  return ampsq;
}





