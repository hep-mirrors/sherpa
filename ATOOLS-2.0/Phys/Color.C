#include "Color.H"

#include "Exception.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"

using namespace ATOOLS;

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

bool CNumber::Evaluate(Expression *const expression)
{
  bool evaluated(false);
  for (Expression::Color_Term_Vector::iterator 
	 cit(expression->begin());cit!=expression->end();++cit) {
    if (*cit!=this && (*cit)->Type()==ctt::number) {
      m_n*=((CNumber*)*cit)->m_n;
      delete *cit;
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
  return new CNumber(m_n);
}

bool Delta::Evaluate(Expression *const expression)
{
  Expression::Color_Term_Vector::iterator end(expression->end());
  if (m_i==m_j) {
    for (Expression::Color_Term_Vector::iterator 
	   tit(expression->begin()); tit!=end;++tit) {
      if (*tit==this) {
	delete this;
	*tit = new CNumber(Complex(NC,0.0));
	return true;
      }
    }
  }
  for (Expression::Color_Term_Vector::iterator 
	 tit(expression->begin()); tit!=end;++tit) {
    if (*tit!=this && (*tit)->Type()==ctt::delta) {
      Delta *delta((Delta*)*tit);
      if (m_j==delta->m_i) {
	m_j=delta->m_j;
	delete delta;
	expression->erase(tit);
	return true;
      }
      else if (m_i==delta->m_j) {
	m_i=delta->m_i;
	delete delta;
	expression->erase(tit);
	return true;
      }
    }
  }
  return false;
}

void Delta::Print() const
{
  msg_Info()<<"("<<this<<"): { d_("<<m_i<<","<<m_j<<") }";
}

Color_Term *Delta::GetCopy() const
{
  return new Delta(m_i,m_j);
}

bool Fundamental::Evaluate(Expression *const expression)
{
  size_t size(expression->size());
  for (size_t j(0);j<size;++j) {
    if ((*expression)[j]==this) {
      for (size_t i(0);i<size;++i) {
	if ((*expression)[i]!=this && 
	    (*expression)[i]->Type()==ctt::fundamental) {
	  Fundamental *fundamental((Fundamental*)(*expression)[i]);
	  if (m_a==fundamental->m_a) {
	    if (!FromF() && !fundamental->FromF()) {
	      Expression *copy(expression->GetCopy());
	      expression->Add(copy);
	      delete (*copy)[j];
	      delete (*copy)[i];
	      (*copy)[j] = new Delta(m_i,m_j);
	      (*copy)[i] = new Delta(fundamental->m_i,fundamental->m_j);
	      copy->push_back(new CNumber(Complex(-0.5/NC,0.0)));
	    }
	    (*expression)[j] = new Delta(m_i,fundamental->m_j);
	    (*expression)[i] = new Delta(fundamental->m_i,m_j);
	    expression->push_back(new CNumber(Complex(0.5,0.0)));
	    delete fundamental;
	    delete this;
	    return true;
	  }
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
  return new Fundamental(m_a,m_i,m_j);
}

bool Adjoint::Evaluate(Expression *const expression)
{
  size_t size(expression->size());
  for (size_t j(0);j<size;++j) {
    if ((*expression)[j]==this) {
      size_t ii(expression->FIndex());
      size_t ij(expression->FIndex());
      size_t ik(expression->FIndex());
      Expression *copy(expression->GetCopy());
      expression->Add(copy);
      delete (*copy)[j];
      (*copy)[j] = new Fundamental(m_a,ii,ij,true);
      copy->push_back(new Fundamental(m_c,ij,ik,true));
      copy->push_back(new Fundamental(m_b,ik,ii,true));
      copy->push_back(new CNumber(Complex(0.0,2.0)));
      (*expression)[j] = new Fundamental(m_a,ii,ij,true);
      expression->push_back(new Fundamental(m_b,ij,ik,true));
      expression->push_back(new Fundamental(m_c,ik,ii,true));
      expression->push_back(new CNumber(Complex(0.0,-2.0)));
      delete this;
      return true;
    }
  }
  return false;
}

void Adjoint::Print() const
{
  msg_Info()<<"("<<this<<"): { f_("<<m_a<<","<<m_b<<","<<m_c<<") }";
}

Color_Term *Adjoint::GetCopy() const
{
  return new Adjoint(m_a,m_b,m_c);
}

Expression::Expression(const std::string &expression): 
  Node<Color_Term*>(NULL,true),
  m_result(0.0,0.0),
  m_findex(0), m_aindex(0), m_evaluated(0)
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
      back() = new Adjoint(a,b,c);
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
      back() = new Fundamental(a,i,j);
      m_findex=Max(m_findex,Max(i,j));
      m_aindex=Max(m_aindex,a);
    }
    else if (factor.find("d_[")==0 && factor[factor.length()-1]==']') {
      size_t cpos(factor.find(','));
      if (cpos==std::string::npos || 
	  factor.find(',',cpos+1)!=std::string::npos)
	THROW(fatal_error,"Invalid number of indices for \\delta.");
      std::string i(factor.substr(3,cpos-3));
      std::string j(factor.substr(cpos+1,factor.length()-cpos-2));
      back() = new Delta(ToType<int>(i),ToType<int>(j));
    }
    else if (factor=="i_") {
      back() = new CNumber(Complex(0.0,1.0));
    }
    else {
      back() = new CNumber(Complex(ToType<double>(factor),0.0));
    }
  }
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

void Expression::PrintStatus(const bool endline,const bool print)
{
  static double oldtime(rpa.gen.Timer().UserTime());
  if (rpa.gen.Timer().UserTime()-oldtime>0.5 || print) {
    oldtime=rpa.gen.Timer().UserTime();
    Expression *root(this);
    while (--*root) root=(Expression*)--*root;
    msg.Out()<<"Terms evaluated: "<<root->Evaluated()<<"     \n"
	     <<"Terms left     : "<<root->Size()<<"     \n"
	     <<"Time  elapsed  : "<<rpa.gen.Timer().UserTime()<<" s     ";
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
    for (Color_Term_Vector::iterator tit(begin());tit!=end();++tit) {
      size_t oldsize((*this)().size());
      if ((*tit)->Evaluate(this)) {
	if ((*this)().size()!=oldsize) {
	  while ((*this)().size()>0) {
	    Expression *expression((Expression*)(*this)().back());
	    if (!expression->Evaluate()) return false;
	    result2+=expression->Result();
	    m_evaluated+=expression->Evaluated();
	    delete expression;
	    (*this)().pop_back();
	  }
	}
	treat=true;
	break;
      }
    }
    if (msg.LevelIsTracking()) PrintStatus(false,false);
  } while (treat==true);
  for (Color_Term_Vector::iterator tit(begin());tit!=end();++tit) {
    if ((*tit)->Type()!=ctt::number) {
      msg.Error()<<"Expression::Evaluate(): Result is nan."<<std::endl;
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
  Expression *expression(new Expression());
  expression->resize(size(),NULL);
  size_t size(size());
  for (size_t i(0);i<size;++i) 
    (*expression)[i] = (*this)[i]->GetCopy();
  expression->m_findex=m_findex;
  expression->m_aindex=m_aindex;
  return expression;
}
