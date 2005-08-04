#include "Color.H"

#include "Exception.H"
#include "Run_Parameter.H"

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
  // test: t-channel gg->qqb -> 16/3
  //  ->    0m0.072s (1.6GHz 32Bit)
  /*
  back() = new Fundamental(this,0,0,1);
  push_back(new Fundamental(this,1,1,2));
  push_back(new Fundamental(this,1,2,3));
  push_back(new Fundamental(this,0,3,0));
  */
  // test: t-u-interference gg->qqb -> -2/3
  //  ->    0m0.074s (1.6GHz 32Bit)
  /*
  back() = new Fundamental(this,0,0,1);
  push_back(new Fundamental(this,1,1,2));
  push_back(new Fundamental(this,0,2,3));
  push_back(new Fundamental(this,1,3,0));
  */
  // test: t-channel 2->2 gluons -> 72 (256 terms)
  //  ->   0m0.081s (1.6GHz 32Bit)
  //  ->   0m0.036s (2GHz 64Bit)
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,4));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,4));
  */
  // test: t-channel 2->3 gluons -> 216 (4096 terms)
  //  ->   0m0.172s (1.6GHz 32Bit)
  //  ->   0m0.059s (2GHz 64Bit)
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,102));
  push_back(new Adjoint(this,102,4,5));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,112));
  push_back(new Adjoint(this,112,4,5));
  */
  // test: t-channel 2->4 gluons -> 648 (65536 terms)
  //  ->   0m1.688s (1.6GHz 32Bit)
  //  ->   0m0.334s (2GHz 64Bit)
  /*
  back() = new Adjoint(this,1,2,101);
  push_back(new Adjoint(this,101,3,102));
  push_back(new Adjoint(this,102,4,103));
  push_back(new Adjoint(this,103,5,6));
  push_back(new Adjoint(this,1,2,111));
  push_back(new Adjoint(this,111,3,112));
  push_back(new Adjoint(this,112,4,113));
  push_back(new Adjoint(this,113,5,6));
  */
  // test: t-channel 2->5 gluons -> 1944 (1048576 terms)
  //  ->   0m31.114s (1.6GHz 32Bit)
  //  ->   0m5.952s  (2GHz 64Bit)
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
  // test: t-channel 2->6 gluons -> 5832 (16777216 terms)
  //  ->   8m51.661s (1.6GHz 32Bit)
  //  ->   1m35.427s (2GHz 64Bit)
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
  // test: t-channel 2->7 gluons -> 17496 (268435456 terms)
  //  ->              (1.6GHz 32Bit)
  //  ->   29m41.262s (2GHz 64Bit)
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
