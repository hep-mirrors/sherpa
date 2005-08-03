#include "Color.H"

#include "Exception.H"

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

Number::Number(Expression *const owner,const double n):
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
	*tit = new Number(expression,NC);
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
	    copy->push_back(new Number(copy,-0.5/NC));
	    (*expression)[j] = new Delta(expression,m_i,fundamental->m_j);
	    (*expression)[i] = new Delta(expression,fundamental->m_i,m_j);
	    expression->push_back(new Number(expression,0.5));
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

bool Adjoint::Evaluate(Expression *const expression)
{
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
      copy->push_back(new Number(copy,-1.0));
      (*expression)[j] = new Fundamental(expression,m_a,ii,ij);
      expression->push_back(new Fundamental(expression,m_b,ij,ik));
      expression->push_back(new Fundamental(expression,m_c,ik,ii));
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
  m_result(0.0),
  m_findex(0), m_aindex(0), m_evaluated(0) {}

Expression::Expression(const std::string &expression): 
  Node<Term*>(NULL,true),
  m_result(0.0),
  m_findex(0), m_aindex(0), m_evaluated(0)
{
  // test: fictious T matrix string
  // ->
  /*
  back() = new Fundamental(this,0,0,1);
  push_back(new Fundamental(this,1,1,2));
  push_back(new Fundamental(this,0,2,3));
  push_back(new Fundamental(this,2,3,4));
  push_back(new Fundamental(this,1,4,5));
  push_back(new Fundamental(this,2,5,0));
  */
  // test: t-channel 2->2 gluons
  // brute force -> 0m0.036s
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,4));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,4));
  */
  // test: t-channel 2->3 gluons
  // brute force -> 0m0.184s
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,102));
  push_back(new Adjoint(this,102,4,5));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,112));
  push_back(new Adjoint(this,112,4,5));
  */
  // test: t-channel 2->4 gluons
  // brute force -> 0m5.419s
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,102));
  push_back(new Adjoint(this,102,4,103));
  push_back(new Adjoint(this,103,5,6));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,112));
  push_back(new Adjoint(this,112,4,113));
  push_back(new Adjoint(this,113,5,6));
  // test: t-channel 2->5 gluons  
  // brute force -> 3m19.615s
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,102));
  push_back(new Adjoint(this,102,4,103));
  push_back(new Adjoint(this,103,5,104));
  push_back(new Adjoint(this,104,6,7));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,112));
  push_back(new Adjoint(this,112,4,113));
  push_back(new Adjoint(this,113,5,114));
  push_back(new Adjoint(this,114,6,7));
  */
  // test: t-channel 2->6 gluons 
  // brute force -> 128m31.599s
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,102));
  push_back(new Adjoint(this,102,4,103));
  push_back(new Adjoint(this,103,5,104));
  push_back(new Adjoint(this,104,6,105));
  push_back(new Adjoint(this,105,7,8));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,112));
  push_back(new Adjoint(this,112,4,113));
  push_back(new Adjoint(this,113,5,114));
  push_back(new Adjoint(this,114,6,115));
  push_back(new Adjoint(this,115,7,8));
  */
  // test: t-channel 2->7 gluons
  // brute force -> 
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,102));
  push_back(new Adjoint(this,102,4,103));
  push_back(new Adjoint(this,103,5,104));
  push_back(new Adjoint(this,104,6,105));
  push_back(new Adjoint(this,105,7,106));
  push_back(new Adjoint(this,106,8,9));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,112));
  push_back(new Adjoint(this,112,4,113));
  push_back(new Adjoint(this,113,5,114));
  push_back(new Adjoint(this,114,6,115));
  push_back(new Adjoint(this,115,7,116));
  push_back(new Adjoint(this,116,8,9));
  */
  m_aindex=10;
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
  size_t terms(size());
  for (size_t i(0);i<(*this)().size();++i) {
    Expression *expression((Expression*)(*this)()[i]);
    terms+=expression->Size();
  }
  return terms;
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
  m_result=1.0;
  double result2(0.0);
  bool treat(false);
  do {
    treat=false;
    for (Term_Vector::iterator tit(begin());tit!=end();++tit) {
      size_t oldsize((*this)().size());
      if ((*tit)->Evaluate(this)) {
	if ((*this)().size()!=oldsize) {
	  Expression *expression((Expression*)(*this)().back());
	  expression->Evaluate();
	  result2+=expression->Result();
	  m_evaluated+=expression->Size()+expression->Evaluated();
	  delete expression;
	  (*this)().pop_back();
	}
	treat=true;
	break;
      }
    }
    if (msg.LevelIsTracking()) {
      Expression *root(this);
      while (--*root) root=(Expression*)--*root;
      msg.Out()<<"Terms evaluated: "<<root->Evaluated()<<"\n"
	       <<"Terms left     : "<<root->Size()
	       <<mm(1,mm::up)<<bm::cr<<std::flush;
    }
  } while (treat==true);
  for (Term_Vector::iterator tit(begin());tit!=end();++tit) {
    if ((*tit)->Type()!=ctt::number) {
      msg.Error()<<"Expression::Evaluate(): Result is nan."<<std::endl;
    }
    else {
      Number *number((Number*)*tit);
      m_result*=(*number)();
    }
  }
  m_result+=result2;
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
