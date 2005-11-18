#include "Color.H"
#include "Exception.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"
#include <iostream>

using namespace ATOOLS;
using namespace AMEGIC;

CNumber::CNumber_Vector CNumber::s_cnumbers;
Delta::Delta_Vector Delta::s_deltas;
Fundamental::Fundamental_Vector Fundamental::s_fundamentals;
Adjoint::Adjoint_Vector Adjoint::s_adjoints;
Trace::Trace_Vector Trace::s_traces;
FullAmplitude_MHV::Amplitude_Vector FullAmplitude_MHV::s_amplitudes;

Expression::Expression_Vector Expression::s_expressions;

class Destructor {
public:

  // destructor
  inline ~Destructor() { Expression::DeleteAll(); }

};// end of class Destructor

static Destructor s_destructor;

Color_Term::~Color_Term() 
{
}

bool Color_Term::Evaluate(Expression *const expression)
{
  THROW(fatal_error,"Virtual function called.");
  return false;
}

void Color_Term::Print() const
{
  THROW(fatal_error,"Virtual function called.");
}

Color_Term *Color_Term::GetCopy() const
{
  THROW(fatal_error,"Virtual function called.");
  return NULL;
}

void Color_Term::Delete()
{
  THROW(fatal_error,"Virtual function called.");
}

bool CNumber::Evaluate(Expression *const expression)
{
  bool evaluated(false);
  for (Expression::Color_Term_Vector::iterator 
	 cit(expression->begin());cit!=expression->end() &&
	 (*cit)->Type()==ctt::number;++cit) {
    if (*cit!=this) {
      m_n*=((CNumber*)*cit)->m_n;
      ((CNumber*)*cit)->Delete();
      cit=--expression->erase(cit);
      evaluated=true;
    }
  }
  return evaluated;
}

void CNumber::Print() const
{
  msg_Info()<<"("<<this<<"): { "<<m_n<<" }";
}

Color_Term *CNumber::GetCopy() const
{
  return New(m_n);
}

CNumber *CNumber::New(const Complex &n)
{
  if (s_cnumbers.empty()) return new CNumber(n);
  CNumber *cnumber(s_cnumbers.back());
  s_cnumbers.pop_back();
  cnumber->m_n=n;
  return cnumber;
}

void CNumber::Delete() 
{
  s_cnumbers.push_back(this);
}

void CNumber::DeleteAll()
{
  while (!s_cnumbers.empty()) {
    delete s_cnumbers.back();
    s_cnumbers.pop_back();
  }
}

bool Delta::Evaluate(Expression *const expression)
{
  if (m_i==m_j) {
    Delete();
    (*expression)[expression->CIndex()] = CNumber::New(Complex(NC,0.0));
    return true;
  }
  bool evaluated(false);
  for (Expression::Color_Term_Vector::iterator 
	 tit(expression->begin());tit!=expression->end();++tit) {
    if ((*tit)->Type()==ctt::delta) {
      if (*tit!=this) {
	Delta *delta((Delta*)*tit);
	if (m_j==delta->m_i) {
	  m_j=delta->m_j;
	  delta->Delete();
	  tit=--expression->erase(tit);
	  evaluated=true;
	}
	else if (m_i==delta->m_j) {
	  m_i=delta->m_i;
	  delta->Delete();
	  tit=--expression->erase(tit);
	  evaluated=true;
	}
      }
    }
    else if ((*tit)->Type()==ctt::fundamental) {
      Fundamental *fundamental((Fundamental*)*tit);
      if (m_j==fundamental->m_i) {
	fundamental->m_i=m_i;
	(*expression)[expression->CIndex()]=fundamental;
	expression->erase(tit);
	Delete();
	return true;
      }	
      if (m_i==fundamental->m_j) {
	fundamental->m_j=m_j;
	(*expression)[expression->CIndex()]=fundamental;
	expression->erase(tit);
	Delete();
	return true;
      }	
    }
  }
  return evaluated;
}

void Delta::Print() const
{
  msg_Info()<<"("<<this<<"): { d_("<<m_i<<","<<m_j<<") }";
}

Color_Term *Delta::GetCopy() const
{
  return New(m_i,m_j);
}

Delta *Delta::New(const size_t &i,const size_t &j)
{
  if (s_deltas.empty()) return new Delta(i,j);
  Delta *delta(s_deltas.back());
  s_deltas.pop_back();
  delta->m_i=i;
  delta->m_j=j;
  return delta;
}

void Delta::Delete() 
{
  s_deltas.push_back(this);
}

void Delta::DeleteAll()
{
  while (!s_deltas.empty()) {
    delete s_deltas.back();
    s_deltas.pop_back();
  }
}

bool Fundamental::Evaluate(Expression *const expression)
{
  size_t size(expression->size()), j(expression->CIndex());
  for (size_t i(0);i<size;++i) {
    if ((*expression)[i]->Type()==ctt::fundamental) {
      if ((*expression)[i]!=this) {
	Fundamental *fundamental((Fundamental*)(*expression)[i]);
	if (m_a==fundamental->m_a) {
	  if (m_j==fundamental->m_i) {
	    if (m_i==fundamental->m_j)
	      (*expression)[j] = CNumber::New(NC);
	    else
	      (*expression)[j] = Delta::New(m_i,fundamental->m_j);
	    if (m_fromf || fundamental->m_fromf) 
	      (*expression)[i] = CNumber::New(Complex(0.5*NC,0.0));
	    else
	      (*expression)[i] = CNumber::New(Complex(CF,0.0));
	    fundamental->Delete();
	    Delete();
	    return true;
	  }
	  if (m_i==fundamental->m_j) {
	    (*expression)[j] = Delta::New(fundamental->m_i,m_j);
	    if (m_fromf || fundamental->m_fromf) 
	      (*expression)[i] = CNumber::New(Complex(0.5*NC,0.0));
	    else
	      (*expression)[i] = CNumber::New(Complex(CF,0.0));
	    fundamental->Delete();
	    Delete();
	    return true;
	  }
	  if (!(m_fromf || fundamental->m_fromf)) {
	    Expression *copy(expression->GetCopy());
	    expression->Add(copy);
	    (*copy)[j]->Delete();
	    (*copy)[i]->Delete();
	    (*copy)[j] = Delta::New(m_i,m_j);
	    (*copy)[i] = Delta::New(fundamental->m_i,fundamental->m_j);
	    copy->push_back(CNumber::New(Complex(-0.5/NC,0.0)));
	  }
	  (*expression)[j] = Delta::New(m_i,fundamental->m_j);
	  (*expression)[i] = Delta::New(fundamental->m_i,m_j);
	  expression->push_back(CNumber::New(Complex(0.5,0.0)));
	  fundamental->Delete();
	  Delete();
	  return true;
	}
      }
    }
  }
  return false;
}

void Fundamental::Print() const
{
  msg_Info()<<"("<<this<<"): { t_("<<m_a<<","<<m_i<<","<<m_j<<") }";
}

Color_Term *Fundamental::GetCopy() const
{
  return New(m_a,m_i,m_j,m_fromf);
}

Fundamental *Fundamental::New(const size_t &a,const size_t &i,const size_t &j,
			      const bool &fromf)
{
  if (s_fundamentals.empty()) return new Fundamental(a,i,j,fromf);
  Fundamental *fundamental(s_fundamentals.back());
  s_fundamentals.pop_back();
  fundamental->m_a=a;
  fundamental->m_i=i;
  fundamental->m_j=j;
  fundamental->m_fromf=fromf;
  return fundamental;
}

void Fundamental::Delete() 
{
  s_fundamentals.push_back(this);
}

void Fundamental::DeleteAll()
{
  while (!s_fundamentals.empty()) {
    delete s_fundamentals.back();
    s_fundamentals.pop_back();
  }
}

bool Adjoint::Evaluate(Expression *const expression)
{
  size_t size(expression->size()), j(expression->CIndex());
  for (size_t i(0);i<size;++i) {
    if ((*expression)[i]->Type()==ctt::fundamental) {
      Fundamental *fundamental((Fundamental*)(*expression)[i]);
      if (m_b==fundamental->m_a) {
	std::swap<size_t>(m_b,m_c);
	std::swap<size_t>(m_a,m_b);
      }
      else if (m_a==fundamental->m_a) {
	std::swap<size_t>(m_a,m_b);
	std::swap<size_t>(m_b,m_c);
      }
      if (m_c==fundamental->m_a) {
	size_t im(expression->FIndex());
	Expression *copy(expression->GetCopy());
	expression->Add(copy);
	(*copy)[j]->Delete();
	(*copy)[j] = Fundamental::New(m_a,im,fundamental->m_j,true);
 	((Fundamental*)(*copy)[i])->m_fromf=true;
	((Fundamental*)(*copy)[i])->m_a=m_b;
	((Fundamental*)(*copy)[i])->m_j=im;
	copy->push_back(CNumber::New(Complex(0.0,1.0)));
	(*expression)[j] = Fundamental::New(m_b,im,fundamental->m_j,true);
 	fundamental->m_fromf=true;
	fundamental->m_a=m_a;
	fundamental->m_j=im;
	expression->push_back(CNumber::New(Complex(0.0,-1.0)));
	Delete();
	return true;
      }
    }
  }
  size_t ii(expression->FIndex());
  size_t ij(expression->FIndex());
  size_t ik(expression->FIndex());
  Expression *copy(expression->GetCopy());
  expression->Add(copy);
  (*copy)[j]->Delete();
  (*copy)[j] = Fundamental::New(m_a,ii,ij,true);
  copy->push_back(Fundamental::New(m_c,ij,ik,true));
  copy->push_back(Fundamental::New(m_b,ik,ii,true));
  copy->push_back(CNumber::New(Complex(0.0,2.0)));
  (*expression)[j] = Fundamental::New(m_a,ii,ij,true);
  expression->push_back(Fundamental::New(m_b,ij,ik,true));
  expression->push_back(Fundamental::New(m_c,ik,ii,true));
  expression->push_back(CNumber::New(Complex(0.0,-2.0)));
  Delete();
  return true;
}

void Adjoint::Print() const
{
  msg_Info()<<"("<<this<<"): { f_("<<m_a<<","<<m_b<<","<<m_c<<") }";
}

Color_Term *Adjoint::GetCopy() const
{
  return New(m_a,m_b,m_c);
}

Adjoint *Adjoint::New(const size_t &a,const size_t &b,const size_t &c)
{
  if (s_adjoints.empty()) return new Adjoint(a,b,c);
  Adjoint *adjoint(s_adjoints.back());
  s_adjoints.pop_back();
  adjoint->m_a=a;
  adjoint->m_b=b;
  adjoint->m_c=c;
  return adjoint;
}

void Adjoint::Delete() 
{
  s_adjoints.push_back(this);
}

void Adjoint::DeleteAll()
{
  while (!s_adjoints.empty()) {
    delete s_adjoints.back();
    s_adjoints.pop_back();
  }
}


Trace::Trace(size_t *a,const size_t &i,const size_t &j):
      Color_Term(ctt::trace), m_i(i), m_j(j) {
  ma = new size_t[a[0]+1];
  for (size_t i=0;i<a[0]+1;i++) ma[i]=a[i];
}

Trace::~Trace() { 
  delete [] ma;    
}

bool Trace::Evaluate(Expression *const expression)
{ 
  size_t j(expression->CIndex());
  size_t ii;
  if (this->m_i==0 && this->m_j==0) ii=expression->FIndex();
  else ii=this->m_i;
  size_t ij(expression->FIndex());
  size_t ik;
  (*expression)[j] = Fundamental::New(ma[1],ii,ij,false);
  size_t n;
  for (n=2;n<ma[0];n++) {
    ik=ij;
    ij=expression->FIndex();
    expression->push_back(Fundamental::New(ma[n],ik,ij,false));
  }
  if (this->m_i!=0 || this->m_j!=0) ii=this->m_j;
  expression->push_back(Fundamental::New(ma[n],ij,ii,false));
  Delete();
  return true;
}

void Trace::Print() const
{
  msg.Out()<<"("<<this<<"): { tr_("<<ma[1];
  for (size_t i=2;i<ma[0]+1;i++) {
    msg.Out()<<","<<ma[i];
  }
  msg.Out()<<")";
  if (m_i!=0 || m_j!=0) msg.Out()<<"_("<<m_i<<","<<m_j<<")";
  msg.Out() << " }";
}

Color_Term *Trace::GetCopy() const
{
  return New(ma,m_i,m_j);
}

Trace *Trace::New(size_t *a,const size_t &i,const size_t &j)
{
  if (s_traces.empty()) return new Trace(a,i,j);
  Trace *trace(s_traces.back());
  s_traces.pop_back();
  for (size_t l=0;l<a[0]+1;l++) trace->ma[l]=a[l];
  trace->m_i=i;
  trace->m_j=j;
  return trace;
}

void Trace::Delete() 
{
  s_traces.push_back(this);
}

void Trace::DeleteAll()
{
  while (!s_traces.empty()) {
    delete s_traces.back();
    s_traces.pop_back();
  }
}




//class FullAmplitude_MHV

//constructor

FullAmplitude_MHV::FullAmplitude_MHV(size_t *a,size_t *f,MHV__Calculator *ca,int con,int pr):
  Color_Term(ctt::amplitude), calc(ca), conjugate(con), print(pr) {
  part=a[0]+f[0];
  ma = new size_t [a[0]+1];
  mf = new size_t [f[0]+1];
  for (size_t i=0;i<a[0]+1;i++) ma[i]=a[i]; 
  for (size_t i=0;i<f[0]+1;i++) mf[i]=f[i];
}

//destructor
FullAmplitude_MHV::~FullAmplitude_MHV() {
  delete [] ma;
  delete [] mf;
}


void FullAmplitude_MHV::Permutation_pureg(int p_number,int** p_adr,Expression *const expression) {
  if (p_number) {    
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=p_number; 
      int** perm_adr = new int*[p_number];
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)];
      Permutation_pureg(p_number-1,perm_adr,expression);
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=0;

    size_t j(expression->CIndex());  
    Complex amp=calc->Differential(perm); 
    size_t *mact = new size_t[part+1];
    mact[0]=part;                     
    if (conjugate) {
      amp=conj(amp);
      for (int i=1;i<part+1;i++) mact[i]=ma[perm[part-i]+1];
    }
    else for (int i=1;i<part+1;i++) mact[i]=ma[perm[i-1]+1];

    if ((*expression)[j]==NULL) {
      (*expression)[j] = Trace::New(mact,0,0);
       expression->push_back(CNumber::New(amp));
    } 
    else {
      Expression *copy(expression->GetCopy());
      expression->Add(copy);
      (*copy)[j]->Delete();
      (*copy)[j] = Trace::New(mact,0,0);
      (*copy).back() = CNumber::New(amp);
    }

                                                              
    msg_Info()<<"     perm: ("<<perm[0];                        //print
    for (int i=1;i<part;i++)  msg_Info()<<","<<perm[i];         //print 
    msg_Info()<<")"  ;                                          //print
    //msg_Info()<<"  "<<amp<<"  ";                              //print
    msg_Info()<<std::endl;                                           //print 


    delete [] mact;
    return;
  }
}

void FullAmplitude_MHV::Permutation_quark2(int p_number,int** p_adr,Expression *const expression) {
  if (p_number) {
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=ma[p_number+1]-1;
      int** perm_adr = new int*[p_number]; 
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)]; 
      Permutation_quark2(p_number-1,perm_adr,expression); 
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=ma[1]-1;

    size_t j(expression->CIndex());  
    Complex amp=calc->Differential(perm);
    size_t *mact = new size_t[part-1];
    size_t mf1,mf2;
    mact[0]=ma[0];                      
    if (conjugate) {
      amp=conj(amp);
      for (size_t i=mf[2];i<mf[2]+ma[0];i++) mact[ma[0]-i+mf[2]]=perm[i%part];
      mf1=mf[2]; mf2=mf[1];
    }
    else {
      for (size_t i=mf[2];i<mf[2]+ma[0];i++) mact[i-mf[2]+1]=perm[i%part];
      mf1=mf[1]; mf2=mf[2];
    }
    if ((*expression)[j]==NULL) {
      (*expression)[j] = Trace::New(mact,mf1,mf2);
       expression->push_back(CNumber::New(amp));
    } 
    else {
      Expression *copy(expression->GetCopy());
      expression->Add(copy);
      (*copy)[j]->Delete();
      (*copy)[j] = Trace::New(mact,mf1,mf2);
      (*copy).back() = CNumber::New(amp);
    }

    msg_Info()<<"     perm: ("<<perm[0];                        //print
    for (int i=1;i<part;i++)  msg_Info()<<","<<perm[i];         //print
    msg_Info()<<")"  ;                                          //print
    //msg_Info()<<"  "<<amp<<"  ";                              //print
    msg_Info()<<std::endl;                                           //print

    delete [] mact;
    return;
  }
}

void FullAmplitude_MHV::Permutation_quark4(int p_number,int** p_adr,Expression *const expression,int* qlist) {
  if (p_number) {
    for (int l=0;l<p_number+1;l++) {
      *p_adr[l]=ma[p_number+1]-1;
      int** perm_adr = new int*[p_number]; 
      for (int i=0;i<p_number;i++) perm_adr[i]=p_adr[(l+i+1)%(p_number+1)]; 
      Permutation_quark4(p_number-1,perm_adr,expression,qlist); 
      delete [] perm_adr; 
    }
  }
  else {
    (*p_adr[0])=ma[1]-1;  msg.Out()<<std::endl; for(int z=0;z<6;z++) msg.Out()<<perm[z]<<",";msg.Out()<<std::endl;

    size_t j(expression->CIndex());  
    Complex amp=calc->Differential(perm);
    size_t *mact1 = new size_t[(qlist[3]-qlist[1]+part-2)%part]; 
    size_t *mact2 = new size_t[(qlist[1]-qlist[3]+part-2)%part];
    size_t mf1,mf2,mf3,mf4;
    mact1[0]=(qlist[3]-qlist[1]+part-2)%part;
    mact2[0]=(qlist[1]-qlist[3]+part-2)%part;           
    if (conjugate) {
      amp=conj(amp);
      for (size_t i=1;i<mact1[0]+1;i++) mact1[mact1[0]-i+1]=perm[(qlist[1]+i)%part];
      for (size_t i=1;i<mact2[0]+1;i++) mact2[mact2[0]-i+1]=perm[(qlist[3]+i)%part];
      mf1=mf[2]; mf2=mf[1]; 
      mf3=mf[4]; mf4=mf[3];
    }
    else {
      for (size_t i=1;i<mact1[0]+1;i++) mact1[i]=perm[(qlist[1]+i)%part];
      for (size_t i=1;i<mact2[0]+1;i++) mact2[i]=perm[(qlist[3]+i)%part];
      mf1=mf[1]; mf2=mf[2]; 
      mf3=mf[3]; mf4=mf[4];
    }
    if ((*expression)[j]==NULL) {
      if (mact1[0]>1) (*expression)[j] = Trace::New(mact1,mf1,mf2);
      else if (mact1[0]==1) (*expression)[j] = Fundamental::New(mact1[1],mf1,mf2);
      else (*expression)[j] = Delta::New(mf1,mf2);
      if (mact2[0]) expression->push_back(Trace::New(mact2,mf3,mf4));
      else if (mact2[0]==1) (*expression)[j] = Fundamental::New(mact2[1],mf3,mf4);
      else expression->push_back(Delta::New(mf3,mf4)); 
      expression->push_back(CNumber::New(amp));
    } 
    else {
      Expression *copy(expression->GetCopy());
      expression->Add(copy);
      (*copy)[j]->Delete();
      if (mact1[0]>1) (*copy)[j] = Trace::New(mact1,mf1,mf2);
      else if (mact1[0]==1) (*copy)[j] = Fundamental::New(mact1[1],mf1,mf2);
      else (*copy)[j] = Delta::New(mf1,mf2);
      if (mact2[0]>1) (*copy).back() = Trace::New(mact2,mf3,mf4);
      else if (mact2[0]==1) (*copy)[j] = Fundamental::New(mact2[1],mf3,mf4);
      else (*copy).back() = Delta::New(mf3,mf4);
      (*copy).back() = CNumber::New(amp);
    }

    msg_Info()<<"     perm: ("<<perm[0];                        //print
    for (int i=1;i<part;i++)  msg_Info()<<","<<perm[i];         //print
    msg_Info()<<")"  ;                                          //print
    //msg_Info()<<"  "<<amp<<"  ";                              //print
    msg_Info()<<std::endl;                                           //print

    delete [] mact1;
    delete [] mact2;
    return;
  }
}


bool FullAmplitude_MHV::Evaluate(Expression *const expression)
{
  perm = new int[part];
  int *qlist(calc->GetQlist());
  int** perm_adr = new int*[ma[0]]; 
  (*expression)[expression->CIndex()]=NULL;
  if (!mf[0]) {
    perm[part-1] = part-1;
    for (int i=0;i<part-1;i++) perm_adr[i]=&perm[i];
    Permutation_pureg(part-2,perm_adr,expression);
  }
  else if (mf[0]==2) { 
    int qpos=1;
    for (int i=0;i<part;i++) {
      if (i==qlist[qpos]) {
	qpos++;
	perm[i]=i;
      }
      else perm_adr[i+1-qpos]=&perm[i]; 
    }
    Permutation_quark2(part-3,perm_adr,expression);     
  }
  else if (mf[0]==4) {
    int qpos;
    int* plist(calc->GetPlist());
    int* qlistbis = new int[5];
    qlistbis[0]=qlist[0];
    int fq=0;
    if (plist[qlist[1]]<0) fq=1;
    qlistbis[1]=qlist[1+fq];
    for (int j=2;j<part-1;j++) {
      qlistbis[3]=(qlistbis[1]+j)%part;
      qlistbis[2]=(qlistbis[3]-1+part)%part;
      qlistbis[4]=(qlistbis[1]-1+part)%part;
      qpos=1;
      for (int i=0;i<part;i++) {
      if (i==qlistbis[qpos]) {
	qpos++;
	perm[i]=qlist[(qpos+fq-2)%4+1];
      }
      else {
	perm_adr[i+1-qpos]=&perm[i];
      }
      }
      Permutation_quark4(part-5,perm_adr,expression,qlistbis);  
    }
    delete [] qlistbis;
  }

  delete [] perm_adr; 
  delete [] perm;
  Delete();
  return true;
}

void FullAmplitude_MHV::Print() const {}

Color_Term *FullAmplitude_MHV::GetCopy() const
{
  return New(ma,mf,calc,conjugate,print);
}

FullAmplitude_MHV *FullAmplitude_MHV::New(size_t *a,size_t *f,MHV__Calculator *ca,int con,int pr)
{
  if (s_amplitudes.empty()) return new FullAmplitude_MHV(a,f,ca,con,pr);
  FullAmplitude_MHV *amplitude(s_amplitudes.back());
  s_amplitudes.pop_back();
  for (size_t i=0;i<a[0]+1;i++) amplitude->ma[i]=a[i];
  for (size_t i=0;i<f[0]+1;i++) amplitude->mf[i]=f[i];
  amplitude->calc=ca;
  amplitude->conjugate=con;
  amplitude->print=pr;
  return amplitude;
}


void FullAmplitude_MHV::Delete() 
{
  s_amplitudes.push_back(this);
}

void FullAmplitude_MHV::DeleteAll()
{
  while (!s_amplitudes.empty()) {
    delete s_amplitudes.back();
    s_amplitudes.pop_back();
  }
}






Expression::Expression(const std::string &expression): 
  Node<Color_Term*>(NULL,true),
  m_result(0.0,0.0),
  m_findex(0), m_aindex(0), m_evaluated(0), calc(NULL)
{
  if (expression.find('+')!=std::string::npos) 
    THROW(not_implemented,"No read in routine for sums yet.");
  std::string expr(expression);
  for (size_t i(0), mpos(expr.find('*'));
       mpos!=std::string::npos || expr.length()>0;mpos=expr.find('*')) {
    if (i>0) push_back(NULL);
    ++i;
    std::string factor;
    if (mpos==std::string::npos) {
      factor=expr;
      expr="";
    }
    else {
      factor=expr.substr(0,mpos);
      expr=expr.substr(mpos+1);
    }
    if  (factor.length()==0) THROW(fatal_error,"Missing factor");
    if (factor.find("f_[")==0 && factor[factor.length()-1]==']') {
      size_t c1pos(factor.find(','));
      if (c1pos==std::string::npos)
	THROW(fatal_error,"Invalid number of indices for t.");
      size_t c2pos(factor.find(',',c1pos+1));
      if (c2pos==std::string::npos || 
	  factor.find(',',c2pos+1)!=std::string::npos)
	THROW(fatal_error,"Invalid number of indices for t.");
      size_t a(ToType<int>(factor.substr(3,c1pos-3)));
      size_t b(ToType<int>(factor.substr(c1pos+1,c2pos-c1pos-1)));
      size_t c(ToType<int>(factor.substr(c2pos+1,factor.length()-c2pos-2)));
      back() = Adjoint::New(a,b,c);
      m_aindex=Max(m_aindex,Max(a,Max(b,c)));
    }
    else if (factor.find("t_[")==0 && factor[factor.length()-1]==']') {
      size_t c1pos(factor.find(','));
      if (c1pos==std::string::npos)
	THROW(fatal_error,"Invalid number of indices for t.");
      size_t c2pos(factor.find(',',c1pos+1));
      if (c2pos==std::string::npos || 
	  factor.find(',',c2pos+1)!=std::string::npos)
	THROW(fatal_error,"Invalid number of indices for t.");
      size_t a(ToType<int>(factor.substr(3,c1pos-3)));
      size_t i(ToType<int>(factor.substr(c1pos+1,c2pos-c1pos-1)));
      size_t j(ToType<int>(factor.substr(c2pos+1,factor.length()-c2pos-2)));
      back() = Fundamental::New(a,i,j);
      m_findex=Max(m_findex,Max(i,j));
      m_aindex=Max(m_aindex,a);
    }
    else if (factor.find("tr_[")==0 && factor[factor.length()-1]==']') {
      size_t c1pos(factor.find(','));
      size_t c2pos;
      size_t fpos(factor.find(']'));
      size_t i=0, j=0;
      if (c1pos==std::string::npos || c1pos>fpos)
	THROW(fatal_error,"Invalid number of color indices for tr.");
      size_t l;
      for (l=1;(c1pos!=std::string::npos && c1pos<fpos);l++) {
	c1pos=factor.find(',',c1pos+1);
      }
      size_t *a = new size_t[l+1];
      a[0]=l;
      c1pos=factor.find(',');
      a[1]=ToType<int>(factor.substr(4,c1pos-4));
      for(l=2;l<a[0];l++) {
	c2pos=factor.find(',',c1pos+1);
	a[l]=ToType<int>(factor.substr(c1pos+1,c2pos-c1pos-1));
	c1pos=c2pos;
      }
      a[l]=ToType<int>(factor.substr(c1pos+1,fpos-c1pos-1));
      if (factor.find("]_[")!=std::string::npos) {
	fpos=factor.find(',',fpos+3);
	if (fpos==std::string::npos || 
	    factor.find(',',fpos+1)!=std::string::npos)
	  THROW(fatal_error,"Invalid number of fundamental indices in tr.");
	i = ToType<int>(factor.substr(factor.find("]_[")+3,fpos-3-factor.find("]_[")));
	j = ToType<int>(factor.substr(fpos+1,factor.length()-fpos-2));
      }
      back() = Trace::New(a,i,j);
      m_findex=Max(m_findex,Max(i,j));
      for (l=1;l<a[0]+1;l++) m_aindex=Max(m_aindex,a[l]);
    }
    else if (factor.find("d_[")==0 && factor[factor.length()-1]==']') {
      size_t cpos(factor.find(','));
      if (cpos==std::string::npos || 
	  factor.find(',',cpos+1)!=std::string::npos)
	THROW(fatal_error,"Invalid number of indices for \\delta.");
      std::string i(factor.substr(3,cpos-3));
      std::string j(factor.substr(cpos+1,factor.length()-cpos-2));
      back() = Delta::New(ToType<int>(i),ToType<int>(j));
    }
    else if (factor=="i_") {
      back() = CNumber::New(Complex(0.0,1.0));
    }
    else {
      back() = CNumber::New(Complex(ToType<double>(factor),0.0));
    }
  }
}


Expression::Expression(int part,Basic_Sfuncs *p_BS,int *hlist,int *plist,int print):
  Node<Color_Term*>(NULL,true),
  m_result(0.0,0.0), m_evaluated(0)
{ 
  calc = new MHV__Calculator(part,p_BS,hlist,plist,print);
  int* qlist(calc->GetQlist());
  m_findex=part+1; 
  m_aindex=part+1;

  size_t *ma = new size_t[part+1-qlist[0]];
  size_t *mf = new size_t[qlist[0]+1];
  ma[0]=part-qlist[0];
  mf[0]=qlist[0];
  int qpos=1;
  for (int i=1;i<part+1;i++)  {
    if (qlist[qpos]==i-1) {
      mf[qpos]=i;
      qpos++;
    }
    else ma[i+1-qpos]=i;
  } 
  if (mf[1]==1 && mf[mf[0]]==(size_t)part) {
    for (size_t i=1;i<mf[0];i++) mf[i]=mf[i+1];
    mf[mf[0]]=1;
  }
  back() = FullAmplitude_MHV::New(ma,mf,calc,0,print);
  push_back(FullAmplitude_MHV::New(ma,mf,calc,1,print));
  delete [] ma;
  delete [] mf;
 }

Expression::~Expression()
{
  for (Color_Term_Vector::iterator tit(begin());tit!=end();++tit) 
    (*tit)->Delete(); 
  if (calc!=NULL) delete calc;
}



size_t Expression::Size()
{
  size_t terms(1);
  for (size_t i(0);i<(*this)().size();++i) {
    Expression *expression((Expression*)(*this)()[i]);
    terms+=expression->Size();
  }
  return terms;
}

static double oldtime(rpa.gen.Timer().UserTime());

void Expression::PrintStatus(const bool endline,const bool print)
{
  double usertime(rpa.gen.Timer().UserTime());
  if (usertime-oldtime>0.5 || print) {
    oldtime=usertime;
    Expression *root(this);
    while (--*root) root=(Expression*)--*root;
    msg.Out()<<"Terms evaluated: "<<root->Evaluated()<<"     \n"
	     <<"Terms left     : "<<root->Size()<<"     \n"
	     <<"Time  elapsed  : "<<usertime<<" s     ";
    if (endline) msg.Out()<<std::endl;
    else msg.Out()<<mm(2,mm::up)<<bm::cr<<std::flush;
  }
}

size_t Expression::Evaluated()
{
  size_t terms(m_evaluated);
  for (size_t i(0);i<(*this)().size();++i) {
    Expression *expression((Expression*)(*this)()[i]);
    terms+=expression->Evaluated();
  }
  return terms;
}

class Order_Type {
public:

  bool operator()(const Color_Term *a,const Color_Term *b)
  {
    return a->Type()<b->Type();
  }

};// end of class Order_Type

bool Expression::Evaluate()
{
  if (--*this==NULL) rpa.gen.Timer().Start();
  m_result=Complex(1.0,0.0);
  if (size()<1 || (*this)[0]==NULL) return false;
  Complex result2(0.0,0.0);
  bool treat(false);
  do {
    treat=false;
    std::sort(begin(),end(),Order_Type());
    m_cindex=0;
    for (Color_Term_Vector::iterator tit(begin());tit!=end();++tit) {
      size_t oldsize((*this)().size());
      if ((*tit)->Evaluate(this)) {
	if ((*this)().size()!=oldsize) {
	  while ((*this)().size()>0) {
	    Expression *expression((Expression*)(*this)().back());
	    if (!expression->Evaluate()) {
	      if (--*this==NULL) expression->Print();
	      m_result=Complex(sqrt(-1.0),sqrt(-1.0));
	      return false;
	    }
	    result2+=expression->Result();
	    m_evaluated+=expression->Evaluated();
	    expression->Delete();
	    (*this)().pop_back();
	  }
	}
	treat=true;
	break;
      }
      ++m_cindex;
    }
    if (msg.LevelIsTracking()) PrintStatus(false,false);
  } while (treat==true);
  for (Color_Term_Vector::iterator tit(begin());tit!=end();++tit) {
    if ((*tit)->Type()!=ctt::number) {
      msg.Error()<<"Expression::Evaluate(): Result is nan."<<std::endl;
      m_result=Complex(sqrt(-1.0),sqrt(-1.0));
      return false;
    }
    else {
      CNumber *number((CNumber*)*tit);
      m_result*=(*number)();
    }
  }
  m_result+=result2;
  m_evaluated+=1;
  if (--*this==NULL) {
    rpa.gen.Timer().Stop();
    if (msg.LevelIsTracking()) PrintStatus();
  }
  return true;
}

void Expression::Print()
{
  msg_Info()<<"("<<this<<"): {\n";
  {
    msg_Indent();
    for (Color_Term_Vector::iterator tit(begin());tit!=end();++tit) {
      (*tit)->Print();
      msg_Info()<<"\n";
    }
  }
  msg_Info()<<"}\n";
  {
    msg_Indent();
    if (operator->()!=NULL) {
      for (size_t i(0);i<(*this)().size();++i) {
	Expression *expression((Expression*)(*this)()[i]);
	expression->Print();
      }
    }
  }
}

Expression *Expression::GetCopy() const
{
  Expression *expression(New(size()));
  size_t esize(size());
  for (size_t i(0);i<esize;++i) 
    (*expression)[i] = (*this)[i]->GetCopy();
  expression->m_findex=m_findex;
  expression->m_aindex=m_aindex;
  return expression;
}

Expression *Expression::New(const size_t &terms)
{
  if (s_expressions.empty()) {
    Expression *expression(new Expression());
    expression->resize(terms,NULL);
    return expression;
  }
  Expression *expression((Expression*)s_expressions.back());
  s_expressions.pop_back();
  expression->resize(terms,NULL);
  return expression;
}

void Expression::Delete() 
{
  for (Color_Term_Vector::iterator tit(begin());
       tit!=end();++tit) (*tit)->Delete();
  resize(0);
  m_evaluated=0;
  s_expressions.push_back(this);
}

void Expression::DeleteAll()
{
  while (!s_expressions.empty()) {
    delete s_expressions.back();
    s_expressions.pop_back();
  }
  CNumber::DeleteAll();
  Delta::DeleteAll();
  Fundamental::DeleteAll();
  Adjoint::DeleteAll();
  Trace::DeleteAll();
  FullAmplitude_MHV::DeleteAll();
}

