#include "AMEGIC++/Amplitude/Zfunctions/Mom.H"

#ifndef Basic_Sfuncs_In_MHV
#include "ATOOLS/Phys/Color.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <iostream.h>
#include "ATOOLS/Math/MathTools.H"


using namespace ATOOLS;
using namespace AMEGIC;


// class Mom

// constructors

Mom::Mom(Vec4D& m,int h,int p):hel(h),part(p),etha(0) 
{
    for (size_t i=0;i<4;i++) (*this)[i]=m[i];
    Compute_l();
} 


Mom::Mom(Vec4D& m, double *et):hel(0),part(0),etha(et)  
{
    for (size_t i=0;i<4;i++) (*this)[i]=m[i]; 
    Compute_l();
} 



// public functions

void Mom::Print() 
{
    msg_Info()<<(*this)<<"    p^2=";
    if (ATOOLS::IsZero(Abs2())) {msg_Info()<<0;}
    else 
    msg_Info()<<Abs2(); 
    msg_Info()<<"    helicity="<<hel<<"    particle ("<<part<<")"<<endl;
}



// private functions

void  Mom::Compute_l() 
{
    if (etha) {	
	l0[0]=((*this)[0]+(*this)[3])*etha[0];
	l0[0]+=Complex((*this)[1],-(*this)[2])*etha[1];
	l0[1]=Complex((*this)[1],(*this)[2])*etha[0];
	l0[1]+=((*this)[0]-(*this)[3])*etha[1];
       	l1[0]=((*this)[0]+(*this)[3])*etha[0];
	l1[0]+=Complex((*this)[1],(*this)[2])*etha[1];
	l1[1]=Complex((*this)[1],-(*this)[2])*etha[0];
	l1[1]+=((*this)[0]-(*this)[3])*etha[1];
    }
    else {
      if (std::abs((*this)[0]+(*this)[3])>1.e+2*Accu()?1:0) {
	l0[0]=csqrt((*this)[0]+(*this)[3]);
	l0[1]=Complex((*this)[1],(*this)[2])/csqrt((*this)[0]+(*this)[3]);
	l1[0]=csqrt((*this)[0]+(*this)[3]);
	l1[1]=Complex((*this)[1],-(*this)[2])/csqrt((*this)[0]+(*this)[3]);
      }
      else {
	l0[0]=0;
	l0[1]=csqrt((*this)[0]-(*this)[3]); 	    
	l1[0]=0;
	l1[1]=csqrt((*this)[0]-(*this)[3]);
      } 
    }
}








// class MomentumList

// constructors

MomentumList::MomentumList(const char*  file):m_propsize(0) 
{
    etha[0]=1; etha[1]=1;
    if (!Get(file))  THROW(fatal_error,"Error in initializing class MomentumList");
    MakeHlist(); 
    MakePlist();
    if (!etha[0] && !etha[1]) {
	etha[0]=0; 
	etha[1]=1;
    }
    if (m_size>5) {
	Vec4D m;
	Mom* momentum; 
	m_propsize=(1<<(m_size-1))-1;
	for (size_t y=m_size;y<m_propsize;y++) {
	    momentum = new Mom(m,etha) ; 
	    this->push_back(momentum);
	}
	BuildMomList();
    }
    s0= new Complex*[size()]; 	
    s1= new Complex*[size()]; 
    for (size_t y=1;y<size();y++) {
	s0[y]= new Complex[y];
	s1[y]= new Complex[y];
    }
    Compute_s(); 
    cout<<"etha = ("<<etha[0]<<","<<etha[1]<<")"<<endl; //print
} 


MomentumList::MomentumList(size_t in, size_t out):hlist(0),plist(0), m_in(in), m_out(out), m_propsize(0) 
{
    etha[0]=1; etha[1]=1;
    m_size=in+out;
    Vec4D m;
    Mom* momentum; 
    for (size_t y=0;y<m_size;y++) {
	momentum = new Mom(m) ; 
	this->push_back(momentum);
    }  
    if (m_size>5) {
	m_propsize=(1<<(m_size-1))-1;
	for (size_t y=m_size;y<m_propsize;y++) {
	    momentum = new Mom(m,etha) ; 
	    this->push_back(momentum);
	}  
    }
    s0= new Complex*[size()]; 	
    s1= new Complex*[size()]; 
    for (size_t y=1;y<size();y++) {
	s0[y]= new Complex[y];
	s1[y]= new Complex[y];
    }
    Compute_s();
} 



// destructor

MomentumList::~MomentumList() 
{  
    for (size_t i=0;i<size();i++) delete (*this)[i];
    if (hlist) delete [] hlist;
    if (plist) delete [] plist; 
    for (size_t y=1;y<size();y++) {  
	delete [] s0[y];
	delete [] s1[y];
    }
    delete [] s0;
    delete [] s1;
}



// public functions

void MomentumList::PutMomenta(const ATOOLS::Vec4D* mom) 
{
    for (size_t y=0;y<m_in;y++)   {
	for (size_t z=0;z<4;z++) (*(*this)[y])[z]=(mom[y])[z];
	(*this)[y]->Compute_l();
    }
    for (size_t y=m_in;y<m_in+m_out;y++)   {
	for (size_t z=0;z<4;z++) (*(*this)[y])[z]=-(mom[y])[z];
	(*this)[y]->Compute_l();
    }
    if (m_size>5) BuildMomList();
    Compute_s();
}


Vec4D MomentumList::Momentum(int mindex) 
{
    Vec4D m;
    for (size_t i=0;i<4;i++) m[i]=(*(*this)[mindex])[i];
    return m;
}


int MomentumList::GetMomNumber(Pfunc* pf) 
{
    size_t res(0); 
    int max(-1);
    for (int y=1;y<pf->argnum;y++) {
	res+=(1<<pf->arg[y]);
	if (pf->arg[y]>max) max=pf->arg[y];
    }
    if (max==(int)m_size-1) {
	res= m_propsize*2+1-res;
	size_t tst=res;
	for (max=-1;tst>0;max++) tst= tst>>1;
    } 
    res+= m_size-2-max;
    return res;
}


int * MomentumList::GetHList() 
{
    if (hlist) return hlist;
    else THROW(fatal_error,"List of helicities is not properly constructed.");
}


int * MomentumList::GetPList() 
{
    if (plist) return plist;
    else THROW(fatal_error,"List of particles is not properly constructed.");
}


void MomentumList::Print() 
{
    for (size_t i=0;i<Size();i++) (*this)[i]->Print();
    if (hlist) { 
	msg_Info()<<endl<<"helicity list: ("<<hlist[0];
	for (size_t i=1;i<Size();i++)  msg_Info()<<","<<hlist[i]; 
	msg_Info()<<")"<<endl; 
    }
    if (plist) {    
	msg_Info()<<"particle list: ("<<plist[0];
	for (size_t i=1;i<Size();i++)  msg_Info()<<","<<plist[i]; 
	msg_Info()<<")"<<endl<<endl; 
    }
}


// private functions

bool MomentumList::Get(const char* file) 
{
    std::string  line;
    My_In_File dat(file);
    if (!dat.Open())
	THROW(fatal_error,"Error in opening a data file");
    while (getline(*dat,line)) {
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
		for (int i=1;i<4;i++) m[i]=m[i]*m[0]/sp;
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
    *dat.close();
    m_size=size();
    if (size()>0) return 1;
    return 0;
}


void MomentumList::MakeHlist() 
{
    hlist = new int[Size()];
    for (size_t i=0;i<Size();i++) hlist[i]=(*this)[i]->hel;
}


void MomentumList::MakePlist() 
{
    plist = new int[Size()];
    int aqpos1=-1, aqtype1=0;
    int aqpos2=-1, aqtype2=0;
    for (size_t i=0;i<Size();i++) {
	plist[i]=(*this)[i]->part;
	if (plist[i]<0 && plist[i]>-10) {
	    if (aqpos1==-1) {aqpos1=i; aqtype1=plist[i];}
	    else if (aqpos2==-1) {aqpos2=i; aqtype2=plist[i];}
	    else THROW(fatal_error,"To many quark-antiquark pairs");
	}
    }
    if (aqpos1>-1 && plist[(aqpos1+1)%Size()]!=-aqtype1) THROW(fatal_error,"Wrong quark-antiquark structure"); 
    if (aqpos2>-1 && plist[(aqpos2+1)%Size()]!=-aqtype2) THROW(fatal_error,"Wrong quark-antiquark structure");
}


void MomentumList::BuildMomList() 
{
    size_t cn(m_size);
    for (size_t i=3;i<m_propsize;i++) {
	size_t kk = 0;
	for (size_t j=0;j<m_size-1;j++) if (i&(1<<j)) kk++;
	if (kk>1) {
	    for (size_t y=0;y<4;y++) (*(*this)[cn])[y]=0;
	    for (size_t j=0;j<m_size-1;j++) {
		if (i&(1<<j)) {
		    for (size_t y=0;y<4;y++) (*(*this)[cn])[y]+=(*(*this)[j])[y];
		}
	    }
	    (*this)[cn]->Compute_l();
	    cn++;
	}
    }
}


void MomentumList::Compute_s() 
{
    for (size_t i=1;i<size();i++) {  
	for (size_t j=0;j<i;j++) {
	    s0[i][j]=(*this)[i]->l0[0] * (*this)[j]->l0[1] - (*this)[i]->l0[1] * (*this)[j]->l0[0]; 
	    s1[i][j]=(*this)[i]->l1[0] * (*this)[j]->l1[1] - (*this)[i]->l1[1] * (*this)[j]->l1[0];
	}
    }
}






// class Fullamplitude_MHV_test


// constructor
Fullamplitude_MHV_test::Fullamplitude_MHV_test(MomentumList *momentuml,int *hl,int *pl): momentumlist(momentuml), hlist(hl), plist(pl) {} 

// destructor
//Fullamplitude_MHV_test::~Fullamplitude_MHV_test() {}


Complex Fullamplitude_MHV_test::Calculate() {
  int qpos=1;
  amp=Complex(0.0,0.0);
  int part(momentumlist->Size());
  calc = new MHVCalculator(part,momentumlist,plist);
  const int* qlist(calc->GetQlist());
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
    Complex am=calc->Differential(perm,hlist);
    amp+=norm(am);
    
    // start 6 gluons TEST    
    if (momentumlist->Size()==6) {
      am=conj(am);
      int perm1[]={0,2,4,1,5,3};
      int permz[6];
      for (int l=0;l<6;l++) permz[l]=perm[perm1[l]];
      Complex am2=calc->Differential(permz,hlist);
      int perm2[]={0,4,2,5,1,3};
      for (int l=0;l<6;l++) permz[l]=perm[perm2[l]];
      am2 +=calc->Differential(permz,hlist);
      int perm3[]={0,2,5,3,1,4};
      for (int l=0;l<6;l++) permz[l]=perm[perm3[l]];
      am2+=calc->Differential(permz,hlist);
      am*=am2;
      am*=2;
      am/=pow(3,2);
      amp+=am;
    }
    // end 6 gluons TEST

    if (1) {                                                              //print
      msg_Info()<<"     perm: ("<<perm[0];                                    //print
      for (size_t i=1;i<momentumlist->Size();i++)  msg_Info()<<","<<perm[i];  //print 
      msg_Info()<<")"  ;                                                      //print
    }                                                                         //print
    if (1) msg_Info()<<"     "<<am<<endl;                                 //print 

    THROW(fatal_error,"Test");        
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
      

    int permtest[]={0,1,2,3,4,5,0,1,2,3,4,5}; 
    for (int y=0;y<6;y++) {
      Complex a=calc->Differential(&permtest[y],hlist); 
      if (1) msg_Info()<<"     "<<a<<endl;                                 //print
    } 
     THROW(fatal_error,"Test"); 
     Complex am=calc->Differential(permtest,hlist);

    // start 2 gluons + 2 quarks TEST    
    if ((momentumlist->Size()-ma[0])==2) {
      if (momentumlist->Size()==4) { 
	amp+= 9*norm(am);

	am=-conj(am);
	int perm1[]={0,1,2,3};
	int permz[4];
	for (int l=0;l<4;l++) permz[l]=perm[perm1[l]];
	Complex am2=calc->Differential(permz,hlist);
	int perm2[]={0,2,1,3};
	for (int l=0;l<4;l++) permz[l]=perm[perm2[l]];
	am2 +=calc->Differential(permz,hlist);
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
	Complex am12=calc->Differential(permz,hlist);
	am12*=2;
	int perm2[]={4,1,2,3,0};
	for (int l=0;l<5;l++) permz[l]=perm[perm2[l]];
	am12 +=calc->Differential(permz,hlist);
	int perm3[]={0,1,2,4,3};
	for (int l=0;l<5;l++) permz[l]=perm[perm3[l]];
	am12 +=calc->Differential(permz,hlist);
	int perm4[]={3,1,2,0,4};
	for (int l=0;l<5;l++) permz[l]=perm[perm4[l]];
	am12 -=calc->Differential(permz,hlist);
	am11*=am12;
	am11*=9;
	amp+=am11;
	
	am=conj(am);
	int perm21[]={0,1,2,3,4};
	for (int l=0;l<5;l++) permz[l]=perm[perm21[l]];
	Complex am22= calc->Differential(permz,hlist);
	int perm22[]={0,1,2,4,3};
	for (int l=0;l<5;l++) permz[l]=perm[perm22[l]];
	am22 +=calc->Differential(permz,hlist);
	int perm23[]={3,1,2,0,4};
	for (int l=0;l<5;l++) permz[l]=perm[perm23[l]];
	am22 +=calc->Differential(permz,hlist);
	int perm24[]={4,1,2,0,3};
	for (int l=0;l<5;l++) permz[l]=perm[perm24[l]];
	am22 +=calc->Differential(permz,hlist); 
	int perm25[]={3,1,2,4,0};
	for (int l=0;l<5;l++) permz[l]=perm[perm25[l]];
	am22 +=calc->Differential(permz,hlist);
	int perm26[]={4,1,2,3,0};
	for (int l=0;l<5;l++) permz[l]=perm[perm26[l]];
	am22 +=calc->Differential(permz,hlist);
	am*=am22;
	amp+=am;
      }
      // end 3 gluons + 2 quarks TEST
    }
    else amp+=norm(am);


    if (1) {                                                              //print
      msg_Info()<<"     perm: ("<<perm[0];                                    //print
      for (size_t i=1;i<momentumlist->Size();i++)  msg_Info()<<","<<perm[i];  //print 
      msg_Info()<<")"<<endl  ;                                                      //print
    }                                                                  //print
    return;
  }
}









// class Fullamplitude_MHV_old

// constructor

Fullamplitude_MHV_old::Fullamplitude_MHV_old(int np,int *pl): 
    calc(0), n_part(np), ampsq(0.,0.), p_norm(1),colorstore(false) 
{
    perm = new int[n_part];
    plist = new int[n_part];
    for (int y=0;y<np;y++) plist[y]=pl[y];
    m_cpl=1;//pow(4.*M_PI*MODEL::s_model->ScalarFunction(std::string("alpha_S"),sqr(rpa->gen.Ecms())),(double)n_part-2.);
} 



//destructor

Fullamplitude_MHV_old::~Fullamplitude_MHV_old() 
{
  delete [] perm;
  delete [] plist;
  if (calc) delete calc;
  if (colorstore) delete permstore;  
}




// Full amplitude 

double Fullamplitude_MHV_old::MSquare(MomentumList *p_BS,int *hlist)
{
  if (!calc) calc = new MHVCalculator(n_part,p_BS,plist); 
  qlist = calc->GetQlist();
  p_hlist=hlist;
  ampsq=0;
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
      p_norm/=4;
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
      p_norm/=16;
      int** perm_adr = new int*[n_part-4];
      permconj = new int[n_part];
      permgl = new int[n_part-4];
      permconjgl = new int[n_part-4];

      const int* plist(calc->GetPlist());
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
  return m_cpl*real(ampsq);
}





Complex Fullamplitude_MHV_old::MSquareHel(MomentumList *p_BS)
{
    int* hlist = new int[n_part]; 
    int tt = (1<<n_part)-1;
    Complex res(0.,0.);
    for (int i=2;i<tt-1;i++) {
	int hps(0);
	for (int j=0;j<n_part;j++) {
	    if (i&(1<<j)) {
		hlist[j]=1; 
		hps+=1;
	    }
	    else hlist[j]=-1;
	}	
	if (hps!=1 && hps!=(n_part-1)) {Complex rest=MSquare(p_BS,hlist); res+=rest;
	msg_Info()<<endl<<"h=("; for (int y=0;y<n_part;y++)  	msg_Info()<<hlist[y]<<",";	msg_Info()<<") = "<<rest  <<endl;
	}
    }
    delete [] hlist;
    return res;
}



// Pure gluonic amplitude
void Fullamplitude_MHV_old::PermutationStoreColor_pureg(int p_number,int** p_adr) {
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


void Fullamplitude_MHV_old::PermutationStoreAmp_pureg(int p_number,int** p_adr) {
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
          
    Complex amp(calc->Differential(perm,p_hlist));
    size_t *perms = new size_t[n_part-1];
    for (int i=0;i<n_part-1;i++) perms[i]=perm[i];
    permstore->PutAmp(perms,amp);
  
    delete [] perms;
    return;
  }
}


void Fullamplitude_MHV_old::Permutation_pureg(int p_number,int** p_adr,bool conjflag) {
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
void Fullamplitude_MHV_old::PermutationStoreColor_quark2(int p_number,std::vector<int>** p_adr) {
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


void Fullamplitude_MHV_old::PermutationStoreAmp_quark2(int p_number,std::vector<int>** p_adr) {
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
  
    Complex amp(calc->Differential(perm,p_hlist));
    permstore->PutAmp(perms,amp);    

    delete [] perms; 
    return;
  }
}


void Fullamplitude_MHV_old::Permutation_quark2(int p_number,int** p_adr,bool conjflag) {
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
void Fullamplitude_MHV_old::Permutation_quark4(int p_number, int** p_adr,bool conjflag) {
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
      const int* plist(calc->GetPlist());
      int fq=0;
      if (plist[qlist[1]]>0) fq=1;
      ampnc = calc->Differential(perm,p_hlist);   

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
	mact1 = new size_t[m_qpos]; 
	mact2 = new size_t[n_part-2-m_qpos];
	mact1[0]=m_qpos-1;
	mact2[0]=n_part-3-m_qpos;	
	for (size_t i=1;i<mact1[0]+1;i++) mact1[i]=perm[i-1]+1;
	for (size_t i=1;i<mact2[0]+1;i++) mact2[i]=perm[m_qpos+i]+1;
	qindex[0]=2;
	qindex[1]=4;

	Permutation_quark4(n_part-5,perm_adr,true); 

	delete [] mact1;
	delete [] mact2;
      }  
      
      int st=perm[m_qpos-1];
      perm[m_qpos-1]=perm[n_part-2];
      perm[n_part-2]=st;

      ampnc = calc->Differential(perm,p_hlist);
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
	mact1 = new size_t[m_qpos]; 
	mact2 = new size_t[n_part-2-m_qpos];
	mact1[0]=m_qpos-1;
	mact2[0]=n_part-3-m_qpos;	
	for (size_t i=1;i<mact1[0]+1;i++) mact1[i]=perm[i-1]+1;
	for (size_t i=1;i<mact2[0]+1;i++) mact2[i]=perm[m_qpos+i]+1;
	qindex[0]=4;
	qindex[1]=2;

	Permutation_quark4(n_part-5,perm_adr,true); 

	delete [] mact1;
	delete [] mact2;
      }  
      st=perm[m_qpos-1];
      perm[m_qpos-1]=perm[n_part-2];
      perm[n_part-2]=st;
      
      delete [] perm_adr;
    }
    else {
      (*p_adr[0])=permconjgl[0];

      Complex amp(ampnc);
      amp*=conj(calc->Differential(permconj,p_hlist));      
     
      Expression expression(n_part,4);
      Expression expression2(n_part,4);
      expression[0] = Trace::New(mact1,1,qindex[0]); 
      expression.push_back(Trace::New(mact2,3,qindex[1]));
      expression2[0] = Trace::New(mact1,1,qindex[0]); 
      expression2.push_back(Trace::New(mact2,3,qindex[1]));

      size_t *mactconj1 = new size_t[m_qposconj]; 
      size_t *mactconj2 = new size_t[n_part-2-m_qposconj];
      mactconj1[0]=m_qposconj-1;
      mactconj2[0]=n_part-3-m_qposconj;	
      for (size_t i=1;i<mactconj1[0]+1;i++) mactconj1[mactconj1[0]+1-i]=permconj[i-1]+1;
      for (size_t i=1;i<mactconj2[0]+1;i++) mactconj2[mactconj2[0]+1-i]=permconj[m_qposconj+i]+1;
      expression.push_back(Trace::New(mactconj1,2,1)); 
      expression.push_back(Trace::New(mactconj2,4,3));
      expression.Evaluate(); 
      amp*=expression.Result();

      delete [] mactconj1;
      delete [] mactconj2;
      
      ampsq+=amp;

      amp=-ampnc/NC;
      int st=permconj[m_qposconj-1];
      permconj[m_qposconj-1]=permconj[n_part-2];
      permconj[n_part-2]=st;
      amp*=conj(calc->Differential(permconj,p_hlist));      
     
      mactconj1 = new size_t[m_qposconj]; 
      mactconj2 = new size_t[n_part-2-m_qposconj];
      mactconj1[0]=m_qposconj-1;
      mactconj2[0]=n_part-3-m_qposconj;	
      for (size_t i=1;i<mactconj1[0]+1;i++) mactconj1[mactconj1[0]+1-i]=permconj[i-1]+1;
      for (size_t i=1;i<mactconj2[0]+1;i++) mactconj2[mactconj2[0]+1-i]=permconj[m_qposconj+i]+1;
      expression2.push_back(Trace::New(mactconj1,4,1)); 
      expression2.push_back(Trace::New(mactconj2,2,3));
      expression2.Evaluate(); 
      amp*=expression2.Result(); 

      st=permconj[m_qposconj-1];
      permconj[m_qposconj-1]=permconj[n_part-2];
      permconj[n_part-2]=st;
      delete [] mactconj1;
      delete [] mactconj2;
     
      ampsq+=amp;
 
      return;
    }
  }
}


#endif





