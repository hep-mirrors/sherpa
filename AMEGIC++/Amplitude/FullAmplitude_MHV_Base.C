#include "AMEGIC++/Amplitude/FullAmplitude_MHV_Base.H"
#include "ATOOLS/Phys/Color.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include <iostream>

using namespace ATOOLS;
using namespace AMEGIC;
using namespace MODEL;
using namespace std;


static map<string,FullAmplitude_MHV_Base*> s_ampmap;


// class FullAmplitude_MHV_Base


// constructor

FullAmplitude_MHV_Base::FullAmplitude_MHV_Base(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl): 
  p_model(model->GetInteractionModel()),m_cpls(cpls),
  p_permstore(0), p_permutation(0), p_calc(0), m_colorstore(0), m_ampstore(0), m_ampstore2(0), p_norm(1), colorflag(false),
  n_part(np), m_plist(0), m_perm(0), m_permgl(0), m_emit(0), m_spect(0), m_dptgluon(false), m_A(1.), m_conv(1.), m_phase2(1.,0.)
{ 
  m_plist= new int[np];
  for (int y=0;y<np;y++) {
    m_plist[y]=pl[y];
    Flavour *fl = new Flavour((kf_code)abs(pl[y]),pl[y]<0);
    m_flist.push_back(fl);
  }
  m_perm= new int[np]; 
  p_calc = new MHVCalculator(n_part,m_plist);
  m_cpl=pow(4.*M_PI*p_model->ScalarFunction(std::string("alpha_S"),rpa.gen.CplScale()),(double)np-2.);
  m_oqcd = (double)n_part-2;
  m_oqed = (double)0;
} 



// destructor

FullAmplitude_MHV_Base::~FullAmplitude_MHV_Base() 
{ 
  if (m_colorstore) { 
    for (int y=0;y<maxn;y++) delete [] m_colorstore[y];
    delete [] m_colorstore;
  }  
  if (m_emit==128) {
    for (ColorStoreMap::iterator it=m_colormap.begin();it!=m_colormap.end();it++) {
      for (int y=0;y<maxn;y++) delete [] (it->second)[y];
      delete [] (it->second);      
    }
  }
  if (p_permutation)      delete p_permutation;
  if (m_ampstore)         delete [] m_ampstore;
  if (m_ampstore2)        delete [] m_ampstore2;
  if (m_perm)             delete [] m_perm; 
  if (m_permgl)           delete [] m_permgl; 
  if (m_plist)            delete [] m_plist;  
  if (p_calc)             delete p_calc;
}


// public functions

double FullAmplitude_MHV_Base::MSquare(int *hlist,MomentumList* BS)   
{  
  m_hlist=hlist;
  if (!colorflag) { 
    InitAmplitude();
    colorflag=true;
  }	
  double res(0.);
  p_calc->SetMomentumList(BS);
  if (m_dptgluon) {
    if (AmpStoreDPT(BS)) res=ResultDPT();     
  }
  else {
    if (AmpStore(BS)) res=Result(m_colorstore);     
  }
  
  if (m_cpls->find("Alpha_QCD")!=m_cpls->end()) p_aqcd=(*m_cpls)["Alpha_QCD"];
  if (m_cpls->find("Alpha_QED")!=m_cpls->end()) p_aqed=(*m_cpls)["Alpha_QED"];
  
  double cplfac(1.0);
  if (p_aqcd && m_oqcd) {
    cplfac *= pow(p_aqcd->Factor(),(double)m_oqcd);
  }
  if (p_aqed && m_oqed) {
    cplfac *= pow(p_aqed->Factor(),(double)m_oqed);
  }
  
  return m_cpl*cplfac*res;
}

void FullAmplitude_MHV_Base::CalculateAmps(int* hlist,MomentumList* BS)
{
  m_hlist=hlist;
  if (!colorflag) { 
    InitAmplitude();
    colorflag=true;
  }	
  p_calc->SetMomentumList(BS);
  AmpStore(BS);
}

double FullAmplitude_MHV_Base::MSquare(int i,int j)   
{  
  if (i+j==0) return m_cpl*Result(m_colorstore);
  return m_cpl*Result(m_colormap[i*100+j]);
}



double FullAmplitude_MHV_Base::MSquareHel(MomentumList* BS) 
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
      double rest=MSquare(hlist,BS); 
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


bool FullAmplitude_MHV_Base::AmpStore(MomentumList* BS) 
{
  THROW(fatal_error,"Virtual function called.");
  return false;
}

bool FullAmplitude_MHV_Base::AmpStoreDPT(MomentumList* BS) 
{
  THROW(fatal_error,"Virtual function called.");
  return false;
}


double FullAmplitude_MHV_Base::Result(Complex** colorstore) 
{
  Complex ampsq(0.,0.);
  for (int y=0;y<maxn;y++) {
    for (int z=0;z<maxn;z++) {     
      Complex amp(m_ampstore[y]);
      amp*=conj(m_ampstore[z]);	
      amp*=colorstore[y][z];
      ampsq+=amp;
    }
  }
  return real(ampsq);
}


double FullAmplitude_MHV_Base::ResultDPT() 
{
  Complex ampsqp(0.,0.);
  Complex ampsqm(0.,0.);
  Complex ampsqpm(0.,0.);

  for (int y=0;y<maxn;y++) {
    for (int z=0;z<maxn;z++) {     
      Complex amp(m_ampstore[y]);
      amp*=conj(m_ampstore[z]);	
      amp*=m_colorstore[y][z];
      ampsqm+=amp;
      amp=m_ampstore2[y];
      amp*=conj(m_ampstore2[z]);	
      amp*=m_colorstore[y][z];
      ampsqp+=amp;
      amp=m_ampstore[y];
      amp*=conj(m_ampstore2[z]);	
      amp*=m_colorstore[y][z];//*m_phase2;
      ampsqpm+=amp;
    }
  }

  return real(0.5*(1.+m_A)*(ampsqm+ampsqp)+m_conv*(1.-m_A)*ampsqpm*m_phase2);
}





// class FullAmplitude_MHV_PureG


// constructor

FullAmplitude_MHV_PureG::FullAmplitude_MHV_PureG(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl):
  FullAmplitude_MHV_Base(model,cpls,np,pl)
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


FullAmplitude_MHV_PureG::FullAmplitude_MHV_PureG(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl,int emit,int spect):
  FullAmplitude_MHV_Base(model,cpls,np,pl)
{ 
  m_emit=emit+1; m_spect=spect+1;
  p_norm=pow((double)2.,(int)n_part);
  p_permutation = new Permutation(n_part-2);
  m_perm[n_part-1] = n_part-1;
  m_perm[n_part-2] = n_part-2; 
  maxn= p_permutation->MaxNumber();

  m_colorstore = new Complex*[maxn];
  for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[maxn];

  if (emit!=127) m_dptgluon = true;
  else {
    for (int i=0;i<n_part-1;i++) {
      for (int j=i+1;j<n_part;j++) {
	Complex** col = new Complex*[maxn];
	for (int y=0;y<maxn;y++) col[y]= new Complex[maxn];
	m_colormap[i*100+j] = col;
      }
    }
  }
  m_ampstore  =  new Complex[maxn];
  m_ampstore2 =  new Complex[maxn];
} 



// destructor

FullAmplitude_MHV_PureG::~FullAmplitude_MHV_PureG() { }



// private functions

void FullAmplitude_MHV_PureG::InitAmplitude()
{
  if (m_emit!=m_spect) ColorStoreDPT(m_emit,m_spect,m_colorstore);
  if (m_emit==128) {
    for (int i=0;i<n_part-1;i++) {
      for (int j=i+1;j<n_part;j++) {
	ColorStoreDPT(i+1,j+1,m_colormap[i*100+j]);
      }
    }
  }
  if (m_emit==m_spect) {
    int** perm_adr = new int*[n_part-2];
    for (int i=0;i<n_part-2;i++) perm_adr[i]=&m_perm[i];
    p_permstore = new PermStore(n_part-2);
    PermutationStoreColor(n_part-3,perm_adr);
    ColorStore();
    delete p_permstore;
    delete [] perm_adr;	
  }
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
    
    Expression expression(2*n_part,0);
    expression[0] = Adjoint::New(n_part,m_perm[0]+1,m_perm[0]+1+n_part);
    for (int i=1;i<n_part-3;i++) expression.push_back(Adjoint::New(m_perm[i-1]+1+n_part,m_perm[i]+1,m_perm[i]+1+n_part));
    expression.push_back(Adjoint::New(m_perm[n_part-4]+1+n_part,m_perm[n_part-3]+1,n_part-1));
    expression.push_back(Adjoint::New(n_part,1,m_perm[0]+1+2*n_part));
    for (int i=1;i<n_part-3;i++) expression.push_back(Adjoint::New(m_perm[i-1]+1+2*n_part,i+1,m_perm[i]+1+2*n_part));
    expression.push_back(Adjoint::New(m_perm[n_part-4]+1+2*n_part,n_part-2,n_part-1));
//     expression.Print();
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

void FullAmplitude_MHV_PureG::ColorStoreDPT(int emit, int spect,Complex **colorstore) 
{
  int* permt;
  size_t *perms = new size_t[n_part];
  size_t *permc = new size_t[n_part-2];
  for (int y=0;y<maxn;y++) {
    permt=p_permutation->Get(y);	    
    for (int i=0;i<n_part-2;i++) perms[i]=permt[i]+1;
    perms[n_part-2]=n_part-1;
    perms[n_part-1]=n_part;
    for (int i=0;i<n_part;i++) if (perms[i]==(size_t)emit||perms[i]==(size_t)spect) perms[i]+=3*n_part;
    for (int z=y;z<maxn;z++) {	    
      permt=p_permutation->Get(z);
      for (int i=0;i<n_part-2;i++) permc[i]=permt[i]+1;

      Expression expression(4*n_part+2,0);
      expression[0] = Adjoint::New(perms[n_part-1],perms[0],perms[0]+n_part);
      for (int i=1;i<n_part-3;i++) expression.push_back(Adjoint::New(perms[i-1]+n_part,perms[i],perms[i]+n_part));
      expression.push_back(Adjoint::New(perms[n_part-4]+n_part,perms[n_part-3],perms[n_part-2]));
      expression.push_back(Adjoint::New(emit+3*n_part,4*n_part+1,emit));

      expression.push_back(Adjoint::New(spect+3*n_part,4*n_part+1,spect));
      expression.push_back(Adjoint::New(n_part,permc[0],permc[0]+2*n_part));
      for (int i=1;i<n_part-3;i++) expression.push_back(Adjoint::New(permc[i-1]+2*n_part,permc[i],permc[i]+2*n_part));
      expression.push_back(Adjoint::New(permc[n_part-4]+2*n_part,permc[n_part-3],n_part-1));
//       expression.Print();
      expression.Evaluate();
      Complex col=expression.Result();
      col/=-4.;

      colorstore[y][z]=col;
      colorstore[z][y]=conj(col);
    }
  }  
  delete [] perms;
  delete [] permc;
  return;
}


bool FullAmplitude_MHV_PureG::AmpStore(MomentumList* BS) 
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


bool FullAmplitude_MHV_PureG::AmpStoreDPT(MomentumList* BS) 
{
  int *permamp;
  if (m_hlist[m_emit-1]!=90) THROW(fatal_error,"FullAmplitude_MHV_PureG::AmpStoreDPT: unexpected helicity label");  
  m_hlist[m_emit-1]=-1;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    for (int z=0;z<n_part-2;z++) m_perm[z]=permamp[z];
    Complex amp(p_calc->Differential(m_perm,m_hlist));
    m_ampstore[y]=amp;
  }
  m_hlist[m_emit-1]=1;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    for (int z=0;z<n_part-2;z++) m_perm[z]=permamp[z];
    Complex amp(p_calc->Differential(m_perm,m_hlist));
    m_ampstore2[y]=amp;
  }
  m_hlist[m_emit-1]=90;
  return true;
}




// class FullAmplitude_MHV_Q2


// constructor

FullAmplitude_MHV_Q2::FullAmplitude_MHV_Q2(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl): 
  FullAmplitude_MHV_Base(model,cpls,np,pl)
{ 
  p_norm=pow((double)2.,(int)n_part-2);
  p_permutation = new Permutation(n_part-2);
  maxn= p_permutation->MaxNumber();
  m_colorstore = new Complex*[maxn];
  for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[maxn];
  m_ampstore =  new Complex[maxn];
  m_permgl = new int[n_part-2];
} 

FullAmplitude_MHV_Q2::FullAmplitude_MHV_Q2(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl,int emit,int spect): 
  FullAmplitude_MHV_Base(model,cpls,np,pl)
{ 
  if (m_flist[0]->IsGluon()&&m_flist[1]->IsQuark()) m_conv=-1.;
  m_emit=emit+1; m_spect=spect+1;

  p_norm=pow((double)2.,(int)n_part-2);
  p_permutation = new Permutation(n_part-2);
  maxn= p_permutation->MaxNumber();

  if (emit!=127) {
    m_efl=Flavour((kf_code)abs(pl[emit]));
    m_sfl=Flavour((kf_code)abs(pl[spect]));
    if (m_efl.IsGluon()) m_dptgluon = true;
    m_colorstore = new Complex*[maxn];
    for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[maxn];
  }
  else {
    m_colorstore = new Complex*[maxn];
    for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[maxn];
    for (int i=0;i<n_part-1;i++) {
      for (int j=i+1;j<n_part;j++) {
	Complex** col = new Complex*[maxn];
	for (int y=0;y<maxn;y++) col[y]= new Complex[maxn];
	m_colormap[i*100+j] = col;
      }
    }
  }
  m_ampstore =  new Complex[maxn];
  if (m_dptgluon) m_ampstore2 =  new Complex[maxn];
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
	 
  int e=m_emit-1;
  int s=m_spect-1;
  if (m_emit==m_spect&&m_emit!=128)  e=s=-1;

  if (qlist[3]>0) { //order is: g...g qbar q
    m_perm[n_part-2]=qlist[2];
    m_perm[n_part-1]=qlist[1];
  }
  else {
    m_perm[n_part-2]=qlist[1];
    m_perm[n_part-1]=qlist[2];
  }
  int cnt=0;
  for (int i=0;i<n_part;i++) {
    if (i!=qlist[1]&&i!=qlist[2]) {
      m_perm[cnt]=i;
      cnt++;
    }
  }

  if (e!=s) ColorStoreDPT(e,s,m_colorstore);
  if (e==127) {
    for (int i=0;i<n_part-1;i++) {
      for (int j=i+1;j<n_part;j++) {
	ColorStoreDPT(i,j,m_colormap[i*100+j]);
      }
    }
  }
  if (e==s) {
    int** perm_adr = new int*[n_part-2];
    for (int i=0;i<n_part-2;i++) perm_adr[i]=&m_permgl[i];
    p_permstore = new PermStore(n_part-2);
    PermutationStoreColor(n_part-3,perm_adr);
    ColorStore();
    delete p_permstore;
    delete [] perm_adr;
  }
  int qpos=1;
  for (int i=0;i<n_part;i++) {
    if (qlist[qpos]==i && qpos<3) qpos++; 
    else m_permgl[i+1-qpos]=i;
  }			 
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

void FullAmplitude_MHV_Q2::ColorStoreDPT(int emit,int spect,Complex **colorstore) 
{
  int e=emit,s=spect;
  for (int i=0;i<n_part;i++) {
    if (m_perm[i]==emit) e=i+1;
    if (m_perm[i]==spect) s=i+1;
  }

  int *permt;
  size_t *mact = new size_t[n_part-1];
  size_t *maconj = new size_t[n_part-1];
  size_t qb=1,q=2;
  Complex fac(1.,0.);
  if (e==n_part-1) qb=3;
  else {
    if (e==n_part) q=4;
    else fac*=Complex(0.,1.);
  }
  if (s==n_part-1) qb=3; 
  else {
    if (s==n_part)q=4;
    else fac*=Complex(0.,1.);
  }
  if (q==4) fac*=-1.;

  mact[0]=maconj[0]=n_part-2;
  for (int y=0;y<maxn;y++) {
    permt=p_permutation->Get(y);
    for (int i=0;i<n_part-2;i++) {
      mact[i+1]=permt[i]+1;
      if (mact[i+1]==(size_t)e||mact[i+1]==(size_t)s) mact[i+1]+=n_part;
    }

    for (int z=y;z<maxn;z++) {	    
      permt=p_permutation->Get(z);
      for (int i=0;i<n_part-2;i++) maconj[i+1]=permt[n_part-3-i]+1;
      Expression expression(2*n_part+1,4);
      expression[0] = Trace::New(mact,q,qb);
      expression.push_back(Trace::New(maconj,1,2));

      if (e<n_part-1) expression.push_back(Adjoint::New(e+n_part,2*n_part,e));
      if (e==n_part-1) expression.push_back(Fundamental::New(2*n_part,qb,1));
      if (e==n_part)   expression.push_back(Fundamental::New(2*n_part,2,q));

      if (s<n_part-1) expression.push_back(Adjoint::New(s+n_part,2*n_part,s));
      if (s==n_part-1) expression.push_back(Fundamental::New(2*n_part,qb,1));
      if (s==n_part)   expression.push_back(Fundamental::New(2*n_part,2,q));

      expression.Evaluate();
//       expression.Print();
      Complex col=expression.Result(); 
      colorstore[y][z]=col*fac;
      colorstore[z][y]=conj(colorstore[y][z]);
    }
  }  
  delete [] mact;
  delete [] maconj;
  return;
}


bool FullAmplitude_MHV_Q2::AmpStore(MomentumList* BS) 
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

bool FullAmplitude_MHV_Q2::AmpStoreDPT(MomentumList* BS) 
{
  int *permamp;
  if (m_hlist[m_emit-1]!=90) THROW(fatal_error,"FullAmplitude_MHV_Q2::AmpStoreDPT: unexpected helicity label");  
  m_hlist[m_emit-1]=-1;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    for (int z=0;z<n_part-2;z++) m_perm[z]=m_permgl[permamp[z]];
    Complex amp(p_calc->Differential(m_perm,m_hlist));
    m_ampstore[y]=amp;
  }
  m_hlist[m_emit-1]=1;
  for (int y=0;y<maxn;y++) {
    permamp=p_permutation->Get(y);
    for (int z=0;z<n_part-2;z++) m_perm[z]=m_permgl[permamp[z]];
    Complex amp(p_calc->Differential(m_perm,m_hlist));
    m_ampstore2[y]=amp;
  }
  m_hlist[m_emit-1]=90;

  return true;
}






// class FullAmplitude_MHV_Q4


// constructor

FullAmplitude_MHV_Q4::FullAmplitude_MHV_Q4(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl): 
  FullAmplitude_MHV_Base(model,cpls,np,pl), p_calc_partner(0)
{ 
  p_norm=pow((double)2.,(int)n_part-4);
  p_permutation = new Permutation(n_part-3);
  maxn= p_permutation->MaxNumber();
  m_colorstore = new Complex*[maxn];
  for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[2*maxn];
  m_ampstore =  new Complex[2*maxn];
  m_permgl = new int[n_part-2];
} 
 
FullAmplitude_MHV_Q4::FullAmplitude_MHV_Q4(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl,int emit,int spect): 
  FullAmplitude_MHV_Base(model,cpls,np,pl), p_calc_partner(0)
{ 
  if (m_flist[0]->IsGluon()&&m_flist[1]->IsQuark()) m_conv=-1.;
  m_emit=emit+1; m_spect=spect+1;

  p_norm=pow((double)2.,(int)n_part-4);
  p_permutation = new Permutation(n_part-3);
  maxn= p_permutation->MaxNumber();

  if (emit!=127) {
    m_efl=Flavour((kf_code)abs(pl[emit]));
    m_sfl=Flavour((kf_code)abs(pl[spect]));
    if (m_efl.IsGluon()) m_dptgluon = true;
    m_colorstore = new Complex*[2*maxn];
    for (int y=0;y<2*maxn;y++) m_colorstore[y]= new Complex[2*maxn];
  }
  else {
    m_colorstore = new Complex*[maxn];
    for (int y=0;y<maxn;y++) m_colorstore[y]= new Complex[2*maxn];
    for (int i=0;i<n_part-1;i++) {
      for (int j=i+1;j<n_part;j++) {
	Complex** col = new Complex*[2*maxn];
	for (int y=0;y<2*maxn;y++) col[y]= new Complex[2*maxn];
	m_colormap[i*100+j] = col;
      }
    }
  }

  m_ampstore =  new Complex[2*maxn];
  if (m_dptgluon) m_ampstore2 =  new Complex[2*maxn];
  m_permgl = new int[n_part-2];
} 
 



// destructor

FullAmplitude_MHV_Q4::~FullAmplitude_MHV_Q4() 
{   
  if (p_calc_partner)     delete p_calc_partner;
  if (m_colorstore) {
    int k=1;
    if (m_emit!=m_spect) k=2;
    for (int y=0;y<k*maxn;y++) delete [] m_colorstore[y];
    delete [] m_colorstore;
    m_colorstore=0;
  }  
  if (m_emit==128) {
    for (ColorStoreMap::iterator it=m_colormap.begin();it!=m_colormap.end();it++) {
      for (int y=maxn;y<2*maxn;y++) delete [] (it->second)[y];
    }
  }
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
    p_calc = new MHVCalculator(n_part,plist2);
    p_calc_partner = new MHVCalculator(n_part,plist_partner);
    qlist=(p_calc->GetQlist());
    delete [] plist2;
    delete [] plist_partner;			
  }
  
  int qpos=1;
  for (int i=0;i<n_part;i++) {
    if (qlist[qpos]==i && qpos<5) qpos++; 
    else m_permgl[i+1-qpos]=i;
  }			 
 
  if (m_emit!=m_spect) ColorStore(m_emit,m_spect,m_colorstore);
  if (m_emit==128) {
    for (int i=0;i<n_part-1;i++) {
      for (int j=i+1;j<n_part;j++) {
	ColorStore(i+1,j+1,m_colormap[i*100+j]);
      }
    }
  }
  if (m_emit==m_spect) ColorStore();
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

void FullAmplitude_MHV_Q4::ColorStore(int e,int s,Complex** colorstore) 
{
  int* permt;
  for (int i=0;i<n_part-4;i++) {
    if (e==m_permgl[i]+1) e=-i-1;
    if (s==m_permgl[i]+1) s=-i-1;
  }
  bool ef=true;
  bool sf=true; 
  for (int i=1;i<3;i++) {
    if (ef && e==m_perm[n_part-i]+1) {e=i;ef=false;}
    if (sf && s==m_perm[n_part-i]+1) {s=i;sf=false;}
  }
  for (int i=3;i<5;i++) {
    if (ef && e==m_permgl[n_part-i]+1) {e=i;ef=false;}
    if (sf && s==m_permgl[n_part-i]+1) {s=i;sf=false;}
  }

  int ql[4];
  for (int i=0;i<4;i++) ql[i]=i+1;
  if (e>0) ql[e-1]+=n_part;
  if (s>0) ql[s-1]+=n_part;
  Complex fac(1.,0.);
  if (e<0) fac*=Complex(0.,1.);
  if (s<0) fac*=Complex(0.,1.);
  if (e%2==0) fac*=-1.;
  if (s%2==0) fac*=-1.;
  for (int y=0;y<2*maxn;y++) {
    permt=p_permutation->Get(y);
    int qpos(-1);
    for (int i=0;qpos<0;i++) if (permt[i]==n_part-4) qpos=i;
    size_t *mact1 = new size_t[qpos+1]; 
    size_t *mact2 = new size_t[n_part-3-qpos];
    mact1[0]=qpos;
    mact2[0]=n_part-4-qpos;
    for (size_t i=1;i<mact1[0]+1;i++) mact1[i]=permt[i-1]+1;
    for (size_t i=1;i<mact2[0]+1;i++) mact2[i]=permt[qpos+i]+1;
    if (ef) {
      for (size_t i=1;i<mact1[0]+1;i++) if (mact1[i]==size_t(-e)) mact1[i]+=n_part;
      for (size_t i=1;i<mact2[0]+1;i++) if (mact2[i]==size_t(-e)) mact2[i]+=n_part;
    }
    if (sf) {
      for (size_t i=1;i<mact1[0]+1;i++) if (mact1[i]==size_t(-s)) mact1[i]+=n_part;
      for (size_t i=1;i<mact2[0]+1;i++) if (mact2[i]==size_t(-s)) mact2[i]+=n_part;
    }

    for (int z=0;z<2*maxn;z++) {  
      permt=p_permutation->Get(z);
      qpos=-1;
      for (int i=0;qpos<0;i++) if (permt[i]==n_part-4) qpos=i;
      size_t *mactconj1 = new size_t[qpos+1]; 
      size_t *mactconj2 = new size_t[n_part-3-qpos];
      mactconj1[0]=qpos;
      mactconj2[0]=n_part-4-qpos;
      for (size_t i=1;i<mactconj1[0]+1;i++) mactconj1[mactconj1[0]+1-i]=permt[i-1]+1;
      for (size_t i=1;i<mactconj2[0]+1;i++) mactconj2[mactconj2[0]+1-i]=permt[qpos+i]+1;
     
      Expression expression1(2*n_part+1,n_part+5);
      if (y>=maxn) {
	expression1[0] = Trace::New(mact1,ql[1],ql[0]); 
	expression1.push_back(Trace::New(mact2,ql[3],ql[2]));
      }
      else {
	expression1[0] = Trace::New(mact1,ql[3],ql[0]); 
	expression1.push_back(Trace::New(mact2,ql[1],ql[2]));
      }
      if (z>=maxn) {
	expression1.push_back(Trace::New(mactconj1,1,2)); 
	expression1.push_back(Trace::New(mactconj2,3,4));
      }
      else {
	expression1.push_back(Trace::New(mactconj1,1,4)); 
	expression1.push_back(Trace::New(mactconj2,3,2));
      }
      if (ef) expression1.push_back(Adjoint::New(n_part-e,2*n_part+1,-e));
      else {
	if (e%2==1) expression1.push_back(Fundamental::New(2*n_part+1,e+n_part,e));
	else expression1.push_back(Fundamental::New(2*n_part+1,e,e+n_part));
      }
      if (sf) expression1.push_back(Adjoint::New(n_part-s,2*n_part+1,-s));
      else {
	if (s%2==1) expression1.push_back(Fundamental::New(2*n_part+1,s+n_part,s));
	else expression1.push_back(Fundamental::New(2*n_part+1,s,s+n_part));
      }
//       expression1.Print();
      expression1.Evaluate(); 
      Complex col=expression1.Result(); 	   
      colorstore[y][z]=col*fac; 
  
      delete [] mactconj1;
      delete [] mactconj2;
    }
    delete [] mact1;
    delete [] mact2; 
  }  

  return;
}


bool FullAmplitude_MHV_Q4::AmpStore(MomentumList* BS) 
{
  int *permamp;
  if (p_calc_partner) p_calc_partner->SetMomentumList(BS);
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

bool FullAmplitude_MHV_Q4::AmpStoreDPT(MomentumList* BS) 
{
  int *permamp;
  if (p_calc_partner) p_calc_partner->SetMomentumList(BS);
  if (m_hlist[m_emit-1]!=90) THROW(fatal_error,"FullAmplitude_MHV_Q4::AmpStoreDPT: unexpected helicity label");  
  m_hlist[m_emit-1]=-1;
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
  m_hlist[m_emit-1]=1;
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
//     else amp *=-1.;
    m_ampstore2[y]=amp; 
    
    m_perm[qpos]=m_perm[n_part-2];
    m_perm[n_part-2]=m_permgl[n_part-4];	
    amp=-p_calc->Differential(m_perm,m_hlist)/NC;   	
    if (p_calc_partner) amp -= p_calc_partner->Differential(m_perm,m_hlist);	
//     else amp*=-1.;
    m_ampstore2[y+maxn]= amp;
    m_perm[n_part-2]=m_perm[qpos];
  }
  m_hlist[m_emit-1]=90;

  return true;
}


double FullAmplitude_MHV_Q4::Result(Complex** colorstore) 
{ 
  if (m_emit!=m_spect) return ResultDPT();
  Complex ampsq(0.,0.);
  if (colorstore!=m_colorstore) {
    for (int y=0;y<2*maxn;y++) {
      for (int z=0;z<2*maxn;z++) {   
	Complex amp(m_ampstore[y]);
	amp*=conj(m_ampstore[z]);	
	amp*=colorstore[y][z];
	ampsq+=amp;
      }
    }
    return real(ampsq);
  }
  for (int y=0;y<maxn;y++) {
    for (int z=0;z<maxn;z++) {   
      Complex amp(m_ampstore[y]);
      amp*=conj(m_ampstore[z]);	
      amp*=colorstore[y][z];
      ampsq+=amp;
      
      amp=m_ampstore[y+maxn];
      amp*=conj(m_ampstore[z]);	
      amp*=colorstore[y][z+maxn];
      ampsq+=amp;

      amp=m_ampstore[y];
      amp*=conj(m_ampstore[z+maxn]);	
      amp*=colorstore[y][z+maxn];
      ampsq+=amp;

      amp=m_ampstore[y+maxn];
      amp*=conj(m_ampstore[z+maxn]);	
      amp*=colorstore[y][z];
      ampsq+=amp;
    }
  }
  return real(ampsq);
}

double FullAmplitude_MHV_Q4::ResultDPT() 
{ 
  Complex ampsqm(0.,0.);
  Complex ampsqp(0.,0.);
  Complex ampsqpm(0.,0.);
  for (int y=0;y<2*maxn;y++) {
    for (int z=0;z<2*maxn;z++) {   
      Complex amp(m_ampstore[y]);
      amp*=conj(m_ampstore[z]);	
      amp*=m_colorstore[y][z];
      ampsqm+=amp;

      if (m_dptgluon) {
	amp=m_ampstore2[y];
	amp*=conj(m_ampstore2[z]);	
	amp*=m_colorstore[y][z];
	ampsqp+=amp;

	amp=m_ampstore[y];
	amp*=conj(m_ampstore2[z]);	
	amp*=m_colorstore[y][z];
	ampsqpm+=amp;
      }
    }
  }
  return real(0.5*(1.+m_A)*(ampsqm+ampsqp)+m_conv*(1.-m_A)*ampsqpm*m_phase2);
}







// class FullAmplitude_MHV_Q2L2


// constructor

FullAmplitude_MHV_Q2L2::FullAmplitude_MHV_Q2L2(Model_Base *model,MODEL::Coupling_Map *const cpls,int np,int *pl): 
  FullAmplitude_MHV_Base(model,cpls,np,pl), m_qlist(0), m_llist(0)
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
  m_oqcd = (double)n_part-4;
  m_oqed = (double)2;
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
    p_calc = new MHVCalculator(n_part,plist2);
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


bool FullAmplitude_MHV_Q2L2::AmpStore(MomentumList* BS) 
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
  int pn = BS->GetMomNumber(&pf);	
  
  if (m_llist[4]==-m_llist[3] && m_qlist[4]==-m_qlist[3]) {
    
    // photon exchange 
    if (m_flist[m_llist[2]]->IsDowntype()) q_em+=q_quark;
    
    // Z exchange
    q_w=Complex((BS->Momentum(pn)).Abs2(),0.);	  
    q_w/= q_w-pow(Flavour(kf_Z).Mass(),2)+Complex(0.,1.)*Flavour(kf_Z).Width()*Flavour(kf_Z).Mass();
    
    if (m_hlist[m_llist[2]]<0) q_w*=(y_lepton/(2*sintw*costw)-q_lepton*sintw/costw);
    else q_w*=(-q_lepton*sintw/costw);
    if (m_hlist[m_qlist[2]]<0) q_w*=(y_quark/(2*sintw*costw)-q_quark*sintw/costw);
    else q_w*=(-q_quark*sintw/costw);
    }

  // W exchange
  else if (m_flist[m_llist[1]]->LeptonFamily()==m_flist[m_llist[2]]->LeptonFamily() && m_hlist[m_llist[2]]<0 && m_hlist[m_qlist[2]]<0) {
    Complex ckm=p_model->ComplexMatrixElement(std::string("CKM"),m_flist[m_qlist[2]]->QuarkFamily()-1,m_flist[m_qlist[1]]->QuarkFamily()-1);
    q_w=Complex((BS->Momentum(pn)).Abs2(),0.);
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

AMEGIC::FullAmplitude_MHV_Base* AMEGIC::FullAmplitude_MHV_Handler(Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist,MomentumList* BS,bool& newamp,int emit,int spec) 
{
  newamp=false;
#ifdef Basic_Sfuncs_In_MHV
  if (part>5) {
    Pfunc_List pl;
    int tt = (1<<(part))-1;
    for (int i=3;i<tt;i++) {
      int kk = 0;
      for (int j=0;j<part;j++) if (i&(1<<j)) kk++;
      if (kk>1) {
	int n=1;
	Pfunc* pf = new Pfunc(kk+1);
	for (int j=0;j<part;j++) if (i&(1<<j)) {
	  pf->arg[n]=j;
	  n++;
	}
	pl.push_back(pf);
      }
    }
    BS->BuildMomlist(pl);
    for (Pfunc_Iterator it=pl.begin();it!=pl.end();it++) delete (*it);
  }
#endif
  string ampID("");
  bool valid = true;
  int fl1=0,fl2=0;
  for (int i=0;i<part;i++) {
    if (abs(plist[i])==21) ampID+="G";
    else {
      if (abs(plist[i])>=1 && abs(plist[i])<=6) {
	if (fl1==0||fl1==abs(plist[i])) {
	  if (plist[i]<0) { ampID+="qb"; fl1=abs(plist[i]);}
	  else            { ampID+="q";  fl1=abs(plist[i]);}
	}
	else {
	  if (fl2==0) fl2=abs(plist[i]);
	  if (fl2!=abs(plist[i])) valid=false;
	  if (plist[i]<0) { ampID+="pb"; }
	  else            { ampID+="p";  }	  
	}
      }
      else ampID+="_"+ToString(plist[i])+"_";
    }
  }
  if (emit!=spec) ampID+=ToString(emit)+ToString(spec);
  if (emit==127) ampID+="vr";
  if (valid) if (s_ampmap.find(ampID)!=s_ampmap.end()) {
    return s_ampmap[ampID];
  }

  FullAmplitude_MHV_Base* fullamp(0);
  MHVCalculator calc(part,plist);
  const int *qlist=(calc.GetQlist());
  
  if (qlist[0]==0) {
    if (emit!=spec||emit==127) fullamp = new FullAmplitude_MHV_PureG(model,cpls,part,plist,emit,spec);
    else fullamp = new FullAmplitude_MHV_PureG(model,cpls,part,plist);                                // pure gluons
  }
  else if (qlist[0]==2) {
    for (int i=1;i<3;i++) {
      if (!Flavour((kf_code)abs(qlist[i+2]),qlist[i+2]<0).IsQuark() || Flavour((kf_code)abs(qlist[i+2]),qlist[i+2]<0).IsMassive()) {
	THROW(fatal_error,"Fullamplitude_MHV_Handler: Amplitude is not implemented");
      }  
    }
    if (emit!=spec||emit==127) fullamp = new FullAmplitude_MHV_Q2(model,cpls,part,plist,emit,spec);
    else fullamp = new FullAmplitude_MHV_Q2(model,cpls,part,plist);                                  // 2 massless quarks
  }
  else if (qlist[0]==4) {
    int nq(0), nl(0);
    for (int i=1;i<5;i++) {
      if (Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsQuark() && !Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsMassive()) nq++;
      else if (Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsLepton() && !Flavour((kf_code)abs(qlist[i+4]),qlist[i+4]<0).IsMassive()) nl++;
      else THROW(fatal_error,"Fullamplitude_MHV_Handler: Amplitude is not implemented"); 
    }
    if (nq==4)  {
      if (emit!=spec||emit==127) fullamp = new FullAmplitude_MHV_Q4(model,cpls,part,plist,emit,spec);
      else fullamp = new FullAmplitude_MHV_Q4(model,cpls,part,plist);                                 // 4 massless quarks
    }
    else if (nq==2 && nl==2 && emit==spec) fullamp = new FullAmplitude_MHV_Q2L2(model,cpls,part,plist);  // 2 massless quarks + 2 leptons
  }

  if (!valid) return fullamp; 
  if (fullamp) {
    newamp = true;
    // s_ampmap[ampID]=fullamp;
    return fullamp;
  }
  else THROW(fatal_error,"Fullamplitude_MHV_Handler: Amplitude is not implemented");
}








