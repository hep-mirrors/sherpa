#include "Color.H"

#include "Exception.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"

using namespace COLOR;
using namespace ATOOLS;

Term::Term(const ctt::type &type,Expression *const owner): 
  m_type(type),
  p_owner(owner) {}

Term::~Term() 
{
}

bool Term::Evaluate(Expression *const expression)
{
  THROW(fatal_error,"Virtual function called.");
  return false;
}

void Term::Print() const
{
  THROW(fatal_error,"Virtual function called.");
}

Term *Term::GetCopy(Expression *const expression) const
{
  THROW(fatal_error,"Virtual function called.");
  return NULL;
}

Number::Number(Expression *const owner,const Complex n):
  Term(ctt::number,owner),
  m_n(n) {}

bool Number::Evaluate(Expression *const expression)
{
  return false;
}

void Number::Print() const
{
  msg_Info()<<"("<<this<<"): { "<<m_n<<" }";
}

Term *Number::GetCopy(Expression *const expression) const
{
  return new Number(expression,m_n);
}

Delta::Delta(Expression *const owner,
	     const size_t i,const size_t j):
  Term(ctt::delta,owner),
  m_i(i), m_j(j) {}

bool Delta::Evaluate(Expression *const expression)
{
  if (m_i==m_j) {
    for (Expression::Term_Vector::iterator tit(expression->begin());
	 tit!=expression->end();++tit) {
      if (*tit==this) {
	delete this;
	*tit = new Number(expression,Complex(NC,0.0));
	return true;
      }
    }
  }
  for (Expression::Term_Vector::iterator tit(expression->begin());
       tit!=expression->end();++tit) {
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

Term *Delta::GetCopy(Expression *const expression) const
{
  return new Delta(expression,m_i,m_j);
}

Fundamental::Fundamental(Expression *const owner,
			 const size_t a,const size_t i,const size_t j):
  Term(ctt::fundamental,owner),
  m_a(a), m_i(i), m_j(j) {}

bool Fundamental::Evaluate(Expression *const expression)
{
  for (size_t j(0);j<expression->size();++j) {
    if ((*expression)[j]==this) {
      for (size_t i(0);i<expression->size();++i) {
	if ((*expression)[i]!=this && 
	    (*expression)[i]->Type()==ctt::fundamental) {
	  Fundamental *fundamental((Fundamental*)(*expression)[i]);
	  if (m_a==fundamental->m_a) {
	    Expression *copy(expression->GetCopy());
	    expression->Add(copy);
	    delete (*copy)[j];
	    delete (*copy)[i];
	    (*copy)[j] = new Delta(copy,m_i,m_j);
	    (*copy)[i] = new Delta(copy,fundamental->m_i,fundamental->m_j);
	    copy->push_back(new Number(copy,Complex(-0.5/NC,0.0)));
	    (*expression)[j] = new Delta(expression,m_i,fundamental->m_j);
	    (*expression)[i] = new Delta(expression,fundamental->m_i,m_j);
	    expression->push_back(new Number(expression,Complex(0.5,0.0)));
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

Term *Fundamental::GetCopy(Expression *const expression) const
{
  return new Fundamental(expression,m_a,m_i,m_j);
}

Adjoint::Adjoint(Expression *const owner,
		 const size_t a,const size_t b,const size_t c):
  Term(ctt::adjoint,owner),
  m_a(a), m_b(b), m_c(c) {}

double Adjoint::SwapIndices(Adjoint *const adjoint)
{
  if (m_a==adjoint->m_b) {
    std::swap<size_t>(adjoint->m_b,adjoint->m_a);
    return -1.0;
  }
  if (m_a==adjoint->m_c) {
    std::swap<size_t>(adjoint->m_c,adjoint->m_a);
    return -1.0;
  }
  if (m_b==adjoint->m_a) {
    std::swap<size_t>(m_b,m_a);
    return -1.0;
  }
  if (m_b==adjoint->m_b) {
    std::swap<size_t>(m_b,m_a);
    std::swap<size_t>(adjoint->m_b,adjoint->m_a);
    return 1.0;
  }
  if (m_b==adjoint->m_c) {
    std::swap<size_t>(m_b,m_a);
    std::swap<size_t>(adjoint->m_c,adjoint->m_a);
    return 1.0;
  }
  if (m_c==adjoint->m_a) {
    std::swap<size_t>(m_c,m_a);
    return -1.0;
  }
  if (m_c==adjoint->m_b) {
    std::swap<size_t>(m_c,m_a);
    std::swap<size_t>(adjoint->m_b,adjoint->m_a);
    return 1.0;
  }
  if (m_c==adjoint->m_c) {
    std::swap<size_t>(m_c,m_a);
    std::swap<size_t>(adjoint->m_c,adjoint->m_a);
    return 1.0;
  }
  return 1.0;
}

void Adjoint::ReplacePair(Adjoint *const adjoint,
			  const size_t &i,const size_t &j,
			  const double &factor,
			  Expression *const expression)
{
  size_t ii(expression->FIndex());
  size_t ij(expression->FIndex());
  size_t ik(expression->FIndex());
  size_t il(expression->FIndex());

  Expression *copy(expression->GetCopy());
  expression->Add(copy);
  delete (*copy)[j];
  delete (*copy)[i];
  (*copy)[j] = new Fundamental(copy,m_b,ii,ij);
  (*copy)[i] = new Fundamental(copy,m_c,ij,ik);
  copy->push_back(new Fundamental(copy,adjoint->m_c,ik,il));
  copy->push_back(new Fundamental(copy,adjoint->m_b,il,ii));
  copy->push_back(new Number(copy,Complex(factor*2.0,0.0)));

  copy = expression->GetCopy();
  expression->Add(copy);
  delete (*copy)[j];
  delete (*copy)[i];
  (*copy)[j] = new Fundamental(copy,m_c,ii,ij);
  (*copy)[i] = new Fundamental(copy,m_b,ij,ik);
  copy->push_back(new Fundamental(copy,adjoint->m_b,ik,il));
  copy->push_back(new Fundamental(copy,adjoint->m_c,il,ii));
  copy->push_back(new Number(copy,Complex(factor*2.0,0.0)));

  copy = expression->GetCopy();
  expression->Add(copy);
  delete (*copy)[j];
  delete (*copy)[i];
  (*copy)[j] = new Fundamental(copy,m_c,ii,ij);
  (*copy)[i] = new Fundamental(copy,m_b,ij,ik);
  copy->push_back(new Fundamental(copy,adjoint->m_c,ik,il));
  copy->push_back(new Fundamental(copy,adjoint->m_b,il,ii));
  copy->push_back(new Number(copy,Complex(-factor*2.0,0.0)));

  (*expression)[j] = new Fundamental(expression,m_b,ii,ij);
  (*expression)[i] = new Fundamental(expression,m_c,ij,ik);
  expression->
    push_back(new Fundamental(expression,adjoint->m_b,ik,il));
  expression->
    push_back(new Fundamental(expression,adjoint->m_c,il,ii));
  expression->
    push_back(new Number(expression,Complex(-factor*2.0,0.0)));
}

bool Adjoint::Evaluate(Expression *const expression)
{
  for (size_t j(0);j<expression->size();++j) {
    if ((*expression)[j]==this) {
      for (size_t i(0);i<expression->size();++i) {
	if ((*expression)[i]!=this && 
	    (*expression)[i]->Type()==ctt::adjoint) {
	  Adjoint *adjoint((Adjoint*)(*expression)[i]);
	  double factor(SwapIndices(adjoint));
	  if (m_a==adjoint->m_a) {
	    ReplacePair(adjoint,i,j,factor,expression);
	    delete adjoint;
	    delete this;
	    return true;
	  }
	}
      }
    }
  }
  for (size_t j(0);j<expression->size();++j) {
    if ((*expression)[j]==this) {
      size_t ii(expression->FIndex());
      size_t ij(expression->FIndex());
      size_t ik(expression->FIndex());
      Expression *copy(expression->GetCopy());
      expression->Add(copy);
      delete (*copy)[j];
      (*copy)[j] = new Fundamental(copy,m_a,ii,ij);
      copy->push_back(new Fundamental(copy,m_c,ij,ik));
      copy->push_back(new Fundamental(copy,m_b,ik,ii));
      copy->push_back(new Number(copy,Complex(0.0,2.0)));
      (*expression)[j] = new Fundamental(expression,m_a,ii,ij);
      expression->push_back(new Fundamental(expression,m_b,ij,ik));
      expression->push_back(new Fundamental(expression,m_c,ik,ii));
      expression->push_back(new Number(expression,Complex(0.0,-2.0)));
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

Term *Adjoint::GetCopy(Expression *const expression) const
{
  return new Adjoint(expression,m_a,m_b,m_c);
}

Expression::Expression(): 
  Node<Term*>(NULL,true),
  m_result(0.0,0.0),
  m_findex(0), m_aindex(0), m_evaluated(0) {}

Expression::Expression(const std::string &expression): 
  Node<Term*>(NULL,true),
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
      if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos)
	THROW(fatal_error,"Invalid number of indices for t.");
      size_t a(ToType<int>(factor.substr(3,c1pos-3)));
      size_t b(ToType<int>(factor.substr(c1pos+1,c2pos-c1pos-1)));
      size_t c(ToType<int>(factor.substr(c2pos+1,factor.length()-c2pos-2)));
      back() = new Adjoint(this,a,b,c);
      m_aindex=Max(m_aindex,Max(a,Max(b,c)));
    }
    else if (factor.find("t_[")==0 && factor[factor.length()-1]==']') {
      size_t c1pos(factor.find(','));
      if (c1pos==std::string::npos)
	THROW(fatal_error,"Invalid number of indices for t.");
      size_t c2pos(factor.find(',',c1pos+1));
      if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos)
	THROW(fatal_error,"Invalid number of indices for t.");
      size_t a(ToType<int>(factor.substr(3,c1pos-3)));
      size_t i(ToType<int>(factor.substr(c1pos+1,c2pos-c1pos-1)));
      size_t j(ToType<int>(factor.substr(c2pos+1,factor.length()-c2pos-2)));
      back() = new Fundamental(this,a,i,j);
      m_findex=Max(m_findex,Max(i,j));
      m_aindex=Max(m_aindex,a);
    }
    else if (factor.find("d_[")==0 && factor[factor.length()-1]==']') {
      size_t cpos(factor.find(','));
      if (cpos==std::string::npos || factor.find(',',cpos+1)!=std::string::npos)
	THROW(fatal_error,"Invalid number of indices for \\delta.");
      std::string i(factor.substr(3,cpos-3));
      std::string j(factor.substr(cpos+1,factor.length()-cpos-2));
      back() = new Delta(this,ToType<int>(i),ToType<int>(j));
    }
    else {
      back() = new Number(this,ToType<double>(factor));
    }
  }
}

Expression::~Expression()
{
  for (Term_Vector::iterator tit(begin());tit!=end();++tit) 
    delete *tit;
}

size_t Expression::Add(Expression *const expression)
{
  (*this)().push_back(expression);
  (*expression)<<this;
  return (*this)().size();
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
    msg.Out()<<"Terms evaluated: "<<root->Evaluated()<<"\n"
	     <<"Terms left     : "<<root->Size()<<"\n"
	     <<"Time  elapsed  : "<<rpa.gen.Timer().UserTime()<<" s";
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

bool Expression::Evaluate()
{
  if (--*this==NULL) rpa.gen.Timer().Start();
  m_result=Complex(1.0,0.0);
  if (size()<1 || (*this)[0]==NULL) return false;
  Complex result2(0.0,0.0);
  bool treat(false);
  do {
    treat=false;
    for (Term_Vector::iterator tit(begin());tit!=end();++tit) {
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
  for (Term_Vector::iterator tit(begin());tit!=end();++tit) {
    if ((*tit)->Type()!=ctt::number) {
      msg.Error()<<"Expression::Evaluate(): Result is nan."<<std::endl;
      return false;
    }
    else {
      Number *number((Number*)*tit);
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
    for (Term_Vector::iterator tit(begin());tit!=end();++tit) {
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
  for (size_t i(0);i<size();++i) {
    (*expression)[i] = (*this)[i]->GetCopy(expression);
  }
  expression->m_findex=m_findex;
  expression->m_aindex=m_aindex;
  return expression;
}
