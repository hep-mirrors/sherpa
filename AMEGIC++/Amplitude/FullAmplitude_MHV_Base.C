#include "FullAmplitude_MHV_Base.H"
#include "Color.H"
#include "Exception.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"
#include "MathTools.H"
#include "Interaction_Model_Base.H"
#include <iostream>

using namespace ATOOLS;
using namespace AMEGIC;
using namespace MODEL;
using namespace std;





// class FullAmplitude_MHV_Base


// constructor

FullAmplitude_MHV_Base::FullAmplitude_MHV_Base(Model_Base *model,int np,int *pl,MomentumList* BS): 
  p_model(model->GetInteractionModel()),
  p_permstore(0), p_permutation(0), p_calc(0), p_BS(BS), m_colorstore(0), m_ampstore(0), p_norm(1), colorflag(false),
  n_part(np), m_plist(0), m_perm(0), m_permgl(0)
{ 
  m_plist= new int[np];
  for (int y=0;y<np;y++) {
    m_plist[y]=pl[y];
    Flavour *fl = new Flavour((kf_code)abs(pl[y]),pl[y]<0);
    m_flist.push_back(fl);
  }
  m_perm= new int[np]; 
  p_calc = new MHVCalculator(n_part,p_BS,m_plist);
  m_cpl=pow(4.*M_PI*p_model->ScalarFunction(std::string("alpha_S"),rpa.gen.CplScale()),(double)np-2.);
  //msg_Info()<<"Amplitude ("; 
  //for (int y=0;y<np-1;y++) msg_Info()<<pl[y]<<",";
  //msg_Info()<<pl[np-1]<<") has been initialized"<<endl;
} 



// destructor

FullAmplitude_MHV_Base::~FullAmplitude_MHV_Base() 
{ 
  if (m_colorstore) { 
    for (int y=0;y<maxn;y++) delete [] m_colorstore[y];
    delete [] m_colorstore;
  }  
  if (p_permutation)      delete p_permutation;
  if (m_ampstore)         delete [] m_ampstore;
  if (m_perm)             delete [] m_perm; 
  if (m_permgl)           delete [] m_permgl; 
  if (m_plist)            delete [] m_plist;  
  if (p_calc)             delete p_calc;
}


// public functions

double FullAmplitude_MHV_Base::MSquare(int *hlist)   
{  
  m_hlist=hlist;
  if (!colorflag) { 
    InitAmplitude();
    colorflag=true;
  }	
  double res(0.);
  if (AmpStore()) res=Result();     
  return m_cpl*res;
}

 

double FullAmplitude_MHV_Base::MSquareHel() 
{
  int* hlist = new int[n_part]; 
  int tt = (1<<(n_part-1))-1;
  double res(0);
  for (int i=2;i<tt;i++) {
    int hps(0);
    for (int j=0;j<n_part;j++) {
      if (i&(1<<j)) {
	hlist[j]=1; 
	hps+=1;
      }
      else hlist[j]=-1;
    }
    if (hps!=1 && hps!=(n_part-1)) {
      double rest=MSquare(hlist); 
      res+=rest;
      msg_Info()<<endl<<i<<"("<<tt-1<<") h=(";
      for (int y=0;y<n_part;y++) msg_Info()<<hlist[y]<<","; 
      msg_Info()<<") = "<<rest<<endl;
    }
  }
  delete [] hlist;
  return 2*res;
}




// private functions

void FullAmplitude_MHV_Base::InitAmplitude()
{
  THROW(fatal_error,"Virtual function called.");
}


bool FullAmplitude_MHV_Base::AmpStore() 
{
  THROW(fatal_error,"Virtual function called.");
  return false;
}


double FullAmplitude_MHV_Base::Result() 
{
  Complex ampsq(0.,0.);
  for (int y=0;y<maxn;y++) {
    for (int z=0;z<maxn;z++) {     
      Complex amp(m_ampstore[y]);
      amp*=conj(m_ampstore[z]);	
      amp*=m_colorstore[y][z];
      ampsq+=amp;
    }
  }
  return real(ampsq);
}





// class FullAmplitude_MHV_PureG


// constructor

FullAmplitude_MHV_PureG::FullAmplitude_MHV_PureG(Model_Base *model,int np,int *pl,MomentumList* BS):
  FullAmplitude_MHV_Base(model,np,pl,BS)
{ 
  p_norm=pow((double)2.,(int)n_part);
  p_permutation = new Permutation(n_part-2);
  m_perm[n_part-1] = n_part-1;
  m_perm[n_part-2] = n_part-2; 
  maxn= p_permutation->MaxNumber();
  m_colorstore = new Complex*[maxn];
  for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[maxn];
  m_ampstore =  new Complex[maxn];
} 



// destructor

FullAmplitude_MHV_PureG::~FullAmplitude_MHV_PureG() { }



// private functions

void FullAmplitude_MHV_PureG::InitAmplitude()
{
  int** perm_adr = new int*[n_part-2];
  for (int i=0;i<n_part-2;i++) perm_adr[i]=&m_perm[i];
  p_permstore = new PermStore(n_part-2);
  PermutationStoreColor(n_part-3,perm_adr);
  ColorStore();
		
  delete p_permstore;
  delete [] perm_adr;	
}

	 	    	  
void FullAmplitude_MHV_PureG::PermutationStoreColor(int p_number,int** p_adr) 
{
  if (p_number) {  
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number; 
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      PermutationStoreColor(p_number-1,perm_adr);
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
      mact[i]=m_perm[i-1]+1;
      maconj[i]=n_part+1-i;
    } 
    Expression expression(2*n_part,0);
    expression[0] = Adjoint::New(n_part,m_perm[0]+1,m_perm[0]+1+n_part);
    for (int i=1;i<n_part-3;i++) expression.push_back(Adjoint::New(m_perm[i-1]+1+n_part,m_perm[i]+1,m_perm[i]+1+n_part));
    expression.push_back(Adjoint::New(m_perm[n_part-4]+1+n_part,m_perm[n_part-3]+1,n_part-1));
    expression.push_back(Adjoint::New(n_part,1,m_perm[0]+1+n_part));
    for (int i=1;i<n_part-3;i++) expression.push_back(Adjoint::New(m_perm[i-1]+1+n_part,i+1,m_perm[i]+1+n_part));
    expression.push_back(Adjoint::New(m_perm[n_part-4]+1+n_part,n_part-2,n_part-1));
    expression.Evaluate();
    Complex col=expression.Result();
    col/=4;
    size_t *perms = new size_t[n_part-2];
    for (int i=0;i<n_part-2;i++) perms[i]=m_perm[i];
    p_permstore->PutColor(perms,col);
    
    delete [] perms;
    return;
  }
}



void FullAmplitude_MHV_PureG::ColorStore() 
{
  int* permt;
  size_t *perms = new size_t[n_part-2];
  size_t *permi = new size_t[n_part-2];
  for (int y=0;y<maxn;y++) {
    permt=p_permutation->Get(y);	    
    for (int i=0;i<n_part-2;i++) permi[permt[i]]=i;	
    for (int z=0;z<maxn;z++) {	    
      permt=p_permutation->Get(z);
      for (int i=0;i<n_part-2;i++) m_perm[i]=permt[n_part-3-i];	
      for (int i=0;i<n_part-2;i++) perms[i]=permi[permt[i]];
      Complex col=p_permstore->GetColor(perms); 
      m_colorstore[z][y]=col;
    }
  }  
  delete [] permi;
  delete [] perms;
  return;
}


bool FullAmplitude_MHV_PureG::AmpStore() 
{
  int *permamp;;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    for (int z=0;z<n_part-2;z++) m_perm[z]=permamp[z];
    Complex amp(p_calc->Differential(m_perm,m_hlist));
    m_ampstore[y]=amp;
  }
  return true;
}





// class FullAmplitude_MHV_Q2


// constructor

FullAmplitude_MHV_Q2::FullAmplitude_MHV_Q2(Model_Base *model,int np,int *pl,MomentumList* BS): 
  FullAmplitude_MHV_Base(model,np,pl,BS)
{ 
  p_norm=pow((double)2.,(int)n_part-2);
  p_permutation = new Permutation(n_part-2);
  maxn= p_permutation->MaxNumber();
  m_colorstore = new Complex*[maxn];
  for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[maxn];
  m_ampstore =  new Complex[maxn];
  m_permgl = new int[n_part-2];
} 



// destructor

FullAmplitude_MHV_Q2::~FullAmplitude_MHV_Q2() { }



// private functions

void FullAmplitude_MHV_Q2::InitAmplitude()
{ 
  const int *qlist=(p_calc->GetQlist());
  for (int i=1;i<3;i++) if (!m_flist[qlist[i]]->IsQuark() || m_flist[qlist[i]]->IsMassive()) {
    THROW(fatal_error,"FullAmplitude_MHV_Q2::InitAmplitude: Amplitude is not implemented");  
  }
	 
  int** perm_adr = new int*[n_part-2];
  
  if (qlist[3]>0) {
    m_perm[n_part-2]=qlist[2];
    m_perm[n_part-1]=qlist[1];
  }
  else {
    m_perm[n_part-2]=qlist[1];
    m_perm[n_part-1]=qlist[2];
  }
  for (int i=0;i<n_part-2;i++) perm_adr[i]=&m_permgl[i];
  
  p_permstore = new PermStore(n_part-2);
  PermutationStoreColor(n_part-3,perm_adr);
  
  int qpos=1;
  for (int i=0;i<n_part;i++) {
    if (qlist[qpos]==i && qpos<3) qpos++; 
    else m_permgl[i+1-qpos]=i;
  }			 
	 
  ColorStore();
	    
  delete p_permstore;
  delete [] perm_adr;
}   
  

void FullAmplitude_MHV_Q2::PermutationStoreColor(int p_number,int** p_adr) 
{
  if (p_number) {  
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number;
      int** perm_adr = new int*[p_number]; 
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      PermutationStoreColor(p_number-1,perm_adr);
      delete [] perm_adr; 
    }
  }
  else {
    *p_adr[0]=0;
    
    size_t *mact = new size_t[n_part-1];
    mact[0]=n_part-2;                
    size_t* maconj = new size_t[n_part-1];
    maconj[0]=n_part-2;
    size_t *perms = new size_t[n_part-2];
    for (int i=0;i<n_part-2;i++) {
      perms[i]=m_permgl[i];
      mact[i+1]=m_permgl[i]+1;
      maconj[n_part-2-i]=i+1;
    }   

    Expression expression(n_part,2);
    expression[0] = Trace::New(mact,1,2);
    expression.push_back(Trace::New(maconj,2,1));
    expression.Evaluate();
    Complex col=expression.Result(); 
    p_permstore->PutColor(perms,col);
    
    delete [] perms;
    delete [] maconj;
    delete [] mact;
    return;
  }
}


void FullAmplitude_MHV_Q2::ColorStore() 
{
  int* permt;
  size_t *perms = new size_t[n_part-2];
  size_t *permi = new size_t[n_part-2];
  for (int y=0;y<maxn;y++) {
    permt=p_permutation->Get(y);
    for (int i=0;i<n_part-2;i++) m_perm[i]=permt[i];
    for (int z=0;z<maxn;z++) {	    
      permt=p_permutation->Get(z);
      for (int i=0;i<n_part-2;i++) permi[permt[i]]=i;
      for (int i=0;i<n_part-2;i++) perms[i]=permi[m_perm[i]];
      Complex col=p_permstore->GetColor(perms); 
      m_colorstore[y][z]=col;
    }
  }  
  delete [] permi;
  delete [] perms;
  return;
}


bool FullAmplitude_MHV_Q2::AmpStore() 
{
  int *permamp;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    for (int z=0;z<n_part-2;z++) m_perm[z]=m_permgl[permamp[z]];
    Complex amp(p_calc->Differential(m_perm,m_hlist));
    m_ampstore[y]=amp;	
  }
  return true;
}






// class FullAmplitude_MHV_Q4


// constructor

FullAmplitude_MHV_Q4::FullAmplitude_MHV_Q4(Model_Base *model,int np,int *pl,MomentumList* BS): 
  FullAmplitude_MHV_Base(model,np,pl,BS), p_calc_partner(0)
{ 
  p_norm=pow((double)2.,(int)n_part-4);
  p_permutation = new Permutation(n_part-3);
  maxn= p_permutation->MaxNumber();
  m_colorstore = new Complex*[maxn];
  for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[2*maxn];
  m_ampstore =  new Complex[2*maxn];
  m_permgl = new int[n_part-2];
} 
 



// destructor

FullAmplitude_MHV_Q4::~FullAmplitude_MHV_Q4() 
{   
  if (p_calc_partner)     delete p_calc_partner;
}



// private functions

void FullAmplitude_MHV_Q4::InitAmplitude()
{
  const int *qlist=(p_calc->GetQlist());
  	  
  int partn(0);
  if (qlist[5]>0) { 
    m_perm[n_part-1]=qlist[1];
    for (int y=2;y<5;y++) {
      if (qlist[y+4]== -qlist[5]) {
	if (!partn)  m_perm[n_part-2]=qlist[y];
	partn++;
      }
      if (qlist[y+4]>0) m_permgl[n_part-3]=qlist[y];
      if (qlist[y+4]<0 && ((qlist[y+4]== -qlist[5] && partn>1) || qlist[y+4]!= -qlist[5])) m_permgl[n_part-4]=qlist[y];
    }     
  }
  else {
    m_perm[n_part-2]=qlist[1];
    for (int y=2;y<5;y++) {
      if (qlist[y+4]== -qlist[5]) {
	if (!partn) m_perm[n_part-1]=qlist[y];
	partn++;
      }
      if (qlist[y+4]>0 && ((qlist[y+4]== -qlist[5] && partn>1) || qlist[y+4]!= -qlist[5])) m_permgl[n_part-3]=qlist[y];
      if (qlist[y+4]<0) m_permgl[n_part-4]=qlist[y];
    }      
  }
  
  // 2 identical quark pairs 
  if (partn==2) {     
    int *plist_partner= new int[n_part];
    int *plist2= new int[n_part];
    for (int y=0;y<n_part;y++) {
      plist_partner[y]=m_plist[y];
      plist2[y]=m_plist[y];
    }	  
    if (qlist[5]>0) {
      plist_partner[qlist[1]]=qlist[5]%6+1; 
      plist2[qlist[1]]=qlist[5]%6+1;
    }
    else  {
      plist_partner[qlist[1]]=-(-qlist[5]%6+1);
      plist2[qlist[1]]=-(-qlist[5]%6+1);
    }
    bool nofpartn(true);
    for (int y=2;y<5;y++) {
      if (qlist[y+4]== -qlist[5]) {
	if (nofpartn)  {
	  plist2[qlist[y]]=-plist2[qlist[1]];
	  nofpartn=false;
	}
	else plist_partner[qlist[y]]=-plist_partner[qlist[1]];
      }
    }
    
    delete p_calc;
    p_calc = new MHVCalculator(n_part,p_BS,plist2);
    p_calc_partner = new MHVCalculator(n_part,p_BS,plist_partner);
    qlist=(p_calc->GetQlist());
    delete [] plist2;
    delete [] plist_partner;			
  }
  
  int qpos=1;
  for (int i=0;i<n_part;i++) {
    if (qlist[qpos]==i && qpos<5) qpos++; 
    else m_permgl[i+1-qpos]=i;
  }			 
  
  ColorStore();
}


void FullAmplitude_MHV_Q4::ColorStore() 
{
  int* permt;
  for (int y=0;y<maxn;y++) {
    permt=p_permutation->Get(y);
    int qpos(-1);
    for (int i=0;qpos<0;i++) if (permt[i]==n_part-4) qpos=i;
    size_t *mact1 = new size_t[qpos+1]; 
    size_t *mact2 = new size_t[n_part-3-qpos];
    mact1[0]=qpos;
    mact2[0]=n_part-4-qpos;
    for (size_t i=1;i<mact1[0]+1;i++) mact1[i]=permt[i-1]+1;
    for (size_t i=1;i<mact2[0]+1;i++) mact2[i]=permt[qpos+i]+1;
    
    for (int z=0;z<maxn;z++) {  
      permt=p_permutation->Get(z);
      qpos=-1;
      for (int i=0;qpos<0;i++) if (permt[i]==n_part-4) qpos=i;
      size_t *mactconj1 = new size_t[qpos+1]; 
      size_t *mactconj2 = new size_t[n_part-3-qpos];
      mactconj1[0]=qpos;
      mactconj2[0]=n_part-4-qpos;
      for (size_t i=1;i<mactconj1[0]+1;i++) mactconj1[mactconj1[0]+1-i]=permt[i-1]+1;
      for (size_t i=1;i<mactconj2[0]+1;i++) mactconj2[mactconj2[0]+1-i]=permt[qpos+i]+1;
      
      Expression expression1(n_part,4);
      expression1[0] = Trace::New(mact1,1,2); 
      expression1.push_back(Trace::New(mact2,3,4));
      expression1.push_back(Trace::New(mactconj1,2,1)); 
      expression1.push_back(Trace::New(mactconj2,4,3));
      expression1.Evaluate(); 
      Complex col=expression1.Result(); 	   
      m_colorstore[y][z]=col; 
      
      Expression expression2(n_part,4);
      expression2[0] = Trace::New(mact1,1,2); 
      expression2.push_back(Trace::New(mact2,3,4));
      expression2.push_back(Trace::New(mactconj1,4,1)); 
      expression2.push_back(Trace::New(mactconj2,2,3));
      expression2.Evaluate(); 
      col=expression2.Result(); 	   
      m_colorstore[y][z+maxn]=col;      
  
      delete [] mactconj1;
      delete [] mactconj2;
    }
    delete [] mact1;
    delete [] mact2; 
  }  
  return;
}


bool FullAmplitude_MHV_Q4::AmpStore() 
{
  int *permamp;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    
    int qpos(-1);
    for (int i=0;qpos<0;i++) if (permamp[i]==n_part-4) qpos=i;
    for (int i=0;i<qpos;i++) m_perm[i]=m_permgl[permamp[i]];
    m_perm[qpos]=m_permgl[n_part-4];	
    m_perm[qpos+1]=m_permgl[n_part-3];
    for (int i=qpos+1;i<n_part-3;i++) m_perm[i+1]=m_permgl[permamp[i]];	 
    Complex amp(p_calc->Differential(m_perm,m_hlist));
    if (p_calc_partner) amp += p_calc_partner->Differential(m_perm,m_hlist)/NC;	
    m_ampstore[y]=amp; 
    
    m_perm[qpos]=m_perm[n_part-2];
    m_perm[n_part-2]=m_permgl[n_part-4];	
    amp= -p_calc->Differential(m_perm,m_hlist)/NC;   	
    if (p_calc_partner) amp -= p_calc_partner->Differential(m_perm,m_hlist);	
    m_ampstore[y+maxn]= amp;
    m_perm[n_part-2]=m_perm[qpos];
  }
  return true;
}


double FullAmplitude_MHV_Q4::Result() 
{ 
  Complex ampsq(0.,0.);
  for (int y=0;y<maxn;y++) {
    for (int z=0;z<maxn;z++) {   
      Complex amp(m_ampstore[y]);
      amp*=conj(m_ampstore[z]);	
      amp*=m_colorstore[y][z];
      ampsq+=amp;
      
      amp=m_ampstore[y+maxn];
      amp*=conj(m_ampstore[z]);	
      amp*=m_colorstore[y][z+maxn];
      ampsq+=amp;

      amp=m_ampstore[y];
      amp*=conj(m_ampstore[z+maxn]);	
      amp*=m_colorstore[y][z+maxn];
      ampsq+=amp;

      amp=m_ampstore[y+maxn];
      amp*=conj(m_ampstore[z+maxn]);	
      amp*=m_colorstore[y][z];
      ampsq+=amp;
    }
  }
  return real(ampsq);
}







// class FullAmplitude_MHV_Q2L2


// constructor

FullAmplitude_MHV_Q2L2::FullAmplitude_MHV_Q2L2(Model_Base *model,int np,int *pl,MomentumList* BS): 
  FullAmplitude_MHV_Base(model,np,pl,BS), m_qlist(0), m_llist(0)
{ 
  m_cpl=pow(4.*M_PI*p_model->ScalarFunction(std::string("alpha_S"),rpa.gen.CplScale()),(double)n_part-4.);
  m_cpl*=4*pow(4.*M_PI*p_model->ScalarFunction(std::string("alpha_QED"),rpa.gen.CplScale()),(double)2.);
  p_norm=pow((double)2.,(int)n_part-4);
  p_permutation = new Permutation(n_part-4);
  maxn= p_permutation->MaxNumber();
  m_colorstore = new Complex*[maxn];
  for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[maxn];
  m_ampstore =  new Complex[maxn];
  m_permgl = new int[n_part-4];
  m_qlist = new int[5];
  m_llist = new int[5];
} 



// destructor

FullAmplitude_MHV_Q2L2::~FullAmplitude_MHV_Q2L2() 
{   
  if (m_qlist)          delete [] m_qlist;    
  if (m_llist)          delete [] m_llist; 
}



// private functions

void FullAmplitude_MHV_Q2L2::InitAmplitude()
{
    const int *qlist=(p_calc->GetQlist());

    int l1(-1),l2(-1),q1(-1),q2(-1);
    for (int i=1;i<5;i++) {
      if (!(qlist[i+4]/6)) {
	if (qlist[i+4]>0) q2=qlist[i];
	else q1=qlist[i];
      }
      else if (qlist[i+4]>0) l2=qlist[i];
      else l1=qlist[i];
    }
    m_qlist[0]=2;
    m_qlist[1]=q1;
    m_qlist[2]=q2;
    m_qlist[3]=m_plist[q1];
    m_qlist[4]=m_plist[q2];
    m_llist[0]=2;
    m_llist[1]=l1;
    m_llist[2]=l2;
    m_llist[3]=m_plist[l1];
    m_llist[4]=m_plist[l2];

    m_perm[n_part-1]=q2; 
    m_perm[n_part-2]=l1; 
    m_perm[n_part-3]=l2; 
    m_perm[n_part-4]=q1;

    int *plist2= new int[n_part];
    for (int y=0;y<n_part;y++) plist2[y]=m_plist[y];
    plist2[q1]=-m_plist[q2];
    plist2[l1]=-m_plist[l2];	    
    delete p_calc;
    p_calc = new MHVCalculator(n_part,p_BS,plist2);
    delete [] plist2;			
    
    int qpos=1;
    for (int i=0;i<n_part;i++) {
      if (qlist[qpos]==i && qpos<5) qpos++; 
      else m_permgl[i+1-qpos]=i;
    }			 
    
    ColorStore();
}


void FullAmplitude_MHV_Q2L2::ColorStore() 
{
  int* permt;
  for (int y=0;y<maxn;y++) {
    permt=p_permutation->Get(y);	
    size_t *mact = new size_t[n_part-3]; 
    mact[0]=n_part-4;
    for (size_t i=1;i<mact[0]+1;i++) mact[i]=permt[i-1]+1;
    
    for (int z=0;z<maxn;z++) {  
      permt=p_permutation->Get(z);
      size_t *mactconj = new size_t[n_part-3]; 
      mactconj[0]=n_part-4;
      for (size_t i=1;i<mactconj[0]+1;i++) mactconj[mactconj[0]+1-i]=permt[i-1]+1;
      
      Expression expression(n_part,4);
      expression[0] = Trace::New(mact,1,2); 
      expression.push_back(Trace::New(mactconj,2,1)); 
      expression.Evaluate(); 
      Complex col=expression.Result(); 	   
      m_colorstore[y][z]=col;    
      
      delete [] mactconj;
    }
    delete [] mact;
  }  
  return;
}


bool FullAmplitude_MHV_Q2L2::AmpStore() 
{
  double q_em(0.);
  Complex q_w(0.,0.);	
  double y_quark(2*m_flist[m_qlist[2]]->IsoWeak()); 
  double y_lepton(2*m_flist[m_llist[2]]->IsoWeak());
  double q_quark(m_flist[m_qlist[2]]->Charge());
  double q_lepton(m_flist[m_llist[2]]->Charge());
  double sintw(sqrt(p_model->ScalarConstant(std::string("sin2_thetaW"))));
  double costw(sqrt(1-p_model->ScalarConstant(std::string("sin2_thetaW"))));
  
  Pfunc pf(3);
  pf.arg[1]=m_llist[1];	
  pf.arg[2]=m_llist[2];
  int pn = p_BS->GetMomNumber(&pf);	
  
  if (m_llist[4]==-m_llist[3] && m_qlist[4]==-m_qlist[3]) {
    
    // photon exchange 
    if (m_flist[m_llist[2]]->IsDowntype()) q_em+=q_quark;
    
    // Z exchange
    q_w=Complex((p_BS->Momentum(pn)).Abs2(),0.);	  
    q_w/= q_w-pow(Flavour(kf_Z).Mass(),2)+Complex(0.,1.)*Flavour(kf_Z).Width()*Flavour(kf_Z).Mass();
    
    if (m_hlist[m_llist[2]]<0) q_w*=(y_lepton/(2*sintw*costw)-q_lepton*sintw/costw);
    else q_w*=(-q_lepton*sintw/costw);
    if (m_hlist[m_qlist[2]]<0) q_w*=(y_quark/(2*sintw*costw)-q_quark*sintw/costw);
    else q_w*=(-q_quark*sintw/costw);
    }

  // W exchange
  else if (m_flist[m_llist[1]]->LeptonFamily()==m_flist[m_llist[2]]->LeptonFamily() && m_hlist[m_llist[2]]<0 && m_hlist[m_qlist[2]]<0) {
    Complex ckm=p_model->ComplexMatrixElement(std::string("CKM"),m_flist[m_qlist[2]]->QuarkFamily()-1,m_flist[m_qlist[1]]->QuarkFamily()-1);
    q_w=Complex((p_BS->Momentum(pn)).Abs2(),0.);
    q_w/= q_w-pow(Flavour(kf_Wplus).Mass(),2)+Complex(0.,1.)*Flavour(kf_Wplus).Width()*Flavour(kf_Wplus).Mass();
    q_w*= ckm/(2*pow(sintw,2));
  }

  Complex q_eff(q_em-q_w);
  if (q_eff==Complex(0.,0.)) return false;
  
  int *permamp;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    for (int i=0;i<n_part-4;i++) m_perm[i]=m_permgl[permamp[i]];
    m_ampstore[y]=q_eff*p_calc->Differential(m_perm,m_hlist); 
  }
  
  return true;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////

AMEGIC::FullAmplitude_MHV_Base* AMEGIC::FullAmplitude_MHV_Handler(Model_Base *model,int part,int* plist,MomentumList* BS) 
{
  FullAmplitude_MHV_Base* fullamp(0);
  MHVCalculator calc(part,BS,plist);
  const int *qlist=(calc.GetQlist());
  
  if (qlist[0]==0) fullamp = new FullAmplitude_MHV_PureG(model,part,plist,BS);  // pure gluons
  else if (qlist[0]==2) {
    for (int i=1;i<3;i++) {
      if (!Flavour((kf_code)abs(qlist[i+2]),qlist[i+2]<0).IsQuark() || Flavour((kf_code)abs(qlist[i+2]),qlist[i+2]<0).IsMassive()) {
	THROW(fatal_error,"Fullamplitude_MHV_Handler: Amplitude is not implemented");
      }  
    }
    fullamp = new FullAmplitude_MHV_Q2(model,part,plist,BS);                    // 2 massless quarks
  }
  else if (qlist[0]==4) {
    int nq(0), nl(0);
    for (int i=1;i<5;i++) {
      if (Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsQuark() && !Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsMassive()) nq++;
      else if (Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsLepton() && !Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsMassive()) nl++;
      else THROW(fatal_error,"Fullamplitude_MHV_Handler: Amplitude is not implemented"); 
    }
    if (nq==4)  fullamp = new FullAmplitude_MHV_Q4(model,part,plist,BS);       // 4 massless quarks
    else if (nq==2 && nl==2) fullamp = new FullAmplitude_MHV_Q2L2(model,part,plist,BS);  // 2 massless quarks + 2 leptons
  }
  
  if (fullamp) return fullamp;
  else THROW(fatal_error,"Fullamplitude_MHV_Handler: Amplitude is not implemented");
}








