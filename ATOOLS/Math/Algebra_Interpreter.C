#include "Algebra_Interpreter.H"

#ifdef USING__Calc_only
#include "Tools.H"
#else
#include "MathTools.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "Message.H"
#endif
#ifndef USING__double_only
#include "Vector.H"
#endif
#include <typeinfo>

using namespace ATOOLS;

#define PTS long unsigned int
#define PT(ARG) (PTS)(ARG)

struct TDouble: public Term {
  double m_value;
  TDouble(const double &value): m_value(value) {}
};// end of struct Double

#ifndef USING__double_only
struct TVec4D: public Term {
  Vec4D m_value;
  TVec4D(const Vec4D &value): m_value(value) {}
};// end of struct Vec4D
#endif

Function::Function(const std::string &tag): 
  m_tag(tag), p_interpreter(NULL) {}

Function::~Function() {}

std::string Function::Evaluate(const std::vector<std::string> &args) const
{
  THROW(fatal_error,"No evaluation rule.");
}

Term *Function::Evaluate(const std::vector<Term*> &args) const
{
  THROW(fatal_error,"No evaluation rule.");
}

class Number: public Function {
private:

  Tag_Replacer *p_replacer;

  bool m_replace;

  double m_svalue, m_sign;

  mutable TDouble m_value;

public:

  Number(const std::string &tag,Tag_Replacer *const replacer);

  Term *Evaluate(const std::vector<Term*> &args) const;

};// end of class Number

Number::Number(const std::string &tag,Tag_Replacer *const replacer): 
  Function(tag), p_replacer(replacer), m_replace(false), 
  m_svalue(0.0), m_sign(1.0), m_value(1.0)
{
  std::string value=tag;
  if (tag[0]=='-') {
    m_sign=-1.0;
    value=value.substr(1);
  }
  m_value.m_tag=value;
  p_replacer->ReplaceTags(value);
  if ((tag[0]!='-'?tag:tag.substr(1))!=value) m_replace=true;
  m_value.m_value=m_svalue=ToType<double>(value);
}

Term *Number::Evaluate(const std::vector<Term*> &args) const
{
  if (args.size()!=0) THROW(fatal_error,"Number requires no argument.");
  std::string tag=m_tag;
  if (m_replace) p_replacer->ReplaceTags(&m_value);
  else m_value.m_value=m_svalue;
  m_value.m_value*=m_sign;
  return &m_value;
}

#ifndef USING__double_only
class Vector: public Function {
private:

  Tag_Replacer *p_replacer;

  bool  m_replace;

  Vec4D m_svalue;

  mutable TVec4D m_value;

public:

  Vector(const std::string &tag,Tag_Replacer *const replacer);

  Term *Evaluate(const std::vector<Term*> &args) const;

};// end of class Vector

Vector::Vector(const std::string &tag,Tag_Replacer *const replacer): 
  Function(tag), p_replacer(replacer), m_replace(false), m_value(Vec4D()) 
{
  m_value.m_tag=tag;
  std::string value=tag;
  p_replacer->ReplaceTags(value);
  if (tag!=value) m_replace=true;
  m_value.m_value=m_svalue=ToType<Vec4D>(value);
}

Term *Vector::Evaluate(const std::vector<Term*> &args) const
{
  if (args.size()!=0) THROW(fatal_error,"Vector requires no argument.");
  std::string tag=m_tag;
  if (m_replace) p_replacer->ReplaceTags(&m_value);
  else m_value.m_value=m_svalue;
  return &m_value;
}
#endif

Operator::~Operator() {}

#define DEFINE_BINARY_DOUBLE_OPERATOR(NAME,TAG,OP,PRIORITY,TYPE)    \
  DEFINE_BINARY_OPERATOR(NAME,TAG,PRIORITY)                         \
  {                                                                 \
    TYPE arg0=ToType<TYPE>(args[0]);                                \
    TYPE arg1=ToType<TYPE>(args[1]);                                \
    return ToString(arg0 OP arg1);                                  \
  }                                                                 \
  Term *NAME::Evaluate(const std::vector<Term*> &args) const        \
  {                                                                 \
    ((TDouble*)args[0])->m_value=(TYPE)((TDouble*)args[0])->m_value \
      OP(TYPE)((TDouble*)args[1])->m_value;                         \
    return args[0];                                                 \
  }

DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Plus,"+",+,12,double)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Minus,"-",-,12,double)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Times,"*",*,13,double)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Divide,"/",/,13,double)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Modulus,"%",%,13,long int)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Shift_Left,"<<",<<,11,long int)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Shift_Right,">>",>>,11,long int)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Logical_And,"&&",&&,5,long int)
DEFINE_BINARY_DOUBLE_OPERATOR(Binary_Logical_Or,"||",||,4,long int)

bool IsAlpha(const std::string& expr) 
{
  for (size_t i=0;i<expr.length();++i) 
    if (expr[i]<48 || expr[i]>57) {
      if (expr[i]=='e' || expr[i]=='E')
	if (i>0 && expr[i-1]>=48 && expr[i-1]<=57) continue;
      if (expr[i]=='.' || expr[i]=='-')
	if ((i>0 && expr[i-1]>=48 && expr[i-1]<=57) ||
	    (i<expr.length()-1 && expr[i+1]>=48 && expr[i+1]<=57))
	  continue;
      return true;
    }
  return false;
}

#define DEFINE_BINARY_SORTABLE_OPERATOR(NAME,TAG,OP,PRIORITY,TYPE)  \
  DEFINE_BINARY_OPERATOR(NAME,TAG,PRIORITY)                         \
  {                                                                 \
    if (IsAlpha(args[0]) || IsAlpha(args[1]))                       \
      return ToString(args[0] OP args[1]);                          \
    TYPE arg0=ToType<TYPE>(args[0]);                                \
    TYPE arg1=ToType<TYPE>(args[1]);                                \
    return ToString(arg0 OP arg1);                                  \
  }                                                                 \
  Term *NAME::Evaluate(const std::vector<Term*> &args) const        \
  {                                                                 \
    ((TDouble*)args[0])->m_value=(TYPE)((TDouble*)args[0])->m_value \
      OP(TYPE)((TDouble*)args[1])->m_value;                         \
    return args[0];                                                 \
  }

DEFINE_BINARY_SORTABLE_OPERATOR(Binary_Equal,"==",==,9,double)
DEFINE_BINARY_SORTABLE_OPERATOR(Binary_Not_Equal,"!=",!=,9,double)
DEFINE_BINARY_SORTABLE_OPERATOR(Binary_Less,"<",<,10,double)
DEFINE_BINARY_SORTABLE_OPERATOR(Binary_Greater,">",>,10,double)
DEFINE_BINARY_SORTABLE_OPERATOR(Binary_Less_Equal,"<=",<=,10,double)
DEFINE_BINARY_SORTABLE_OPERATOR(Binary_Greater_Equal,">=",>=,10,double)

DEFINE_UNARY_OPERATOR(Unary_Not,"!",14)
{
  double arg0=ToType<int>(args[0]);
  return ToString(!arg0);
}

Term *Unary_Not::Evaluate(const std::vector<Term*> &args) const
{
  ((TDouble*)args[0])->m_value=!(int)((TDouble*)args[0])->m_value;
  return args[0];
}

#define DEFINE_BINARY_DOUBLE_FUNCTION(NAME,TAG,OP)                  \
  DEFINE_FUNCTION(NAME,TAG)                                         \
  {                                                                 \
    if (args.size()!=2)                                             \
      THROW(fatal_error,std::string(TAG)+" requires 2 arguments."); \
    double arg0=ToType<double>(args[0]);                            \
    double arg1=ToType<double>(args[1]);                            \
    return ToString(OP(arg0,arg1));                                 \
  }                                                                 \
  Term *NAME::Evaluate(const std::vector<Term*> &args) const        \
  {                                                                 \
    ((TDouble*)args[0])->m_value=OP(((TDouble*)args[0])->m_value,   \
				    ((TDouble*)args[1])->m_value);  \
    return args[0];                                                 \
  }

DEFINE_BINARY_DOUBLE_FUNCTION(Power,"pow",pow)

#define DEFINE_ITERATED_DOUBLE_FUNCTION(NAME,TAG,OP)                 \
  DEFINE_FUNCTION(NAME,TAG)                                          \
  {                                                                  \
    if (args.size()<2)                                               \
      THROW(fatal_error,std::string(TAG)+" requires 2 arguments.");  \
    double arg0=ToType<double>(args[0]);                             \
    for (size_t i=1;i<args.size();++i)                               \
      arg0=OP(arg0,ToType<double>(args[i]));                         \
    return ToString(arg0);                                           \
  }                                                                  \
  Term *NAME::Evaluate(const std::vector<Term*> &args) const         \
  {                                                                  \
    for (size_t i=1;i<args.size();++i)                               \
      ((TDouble*)args[0])->m_value=OP(((TDouble*)args[0])->m_value,  \
	  			      ((TDouble*)args[i])->m_value); \
    return args[0];                                                  \
  }

DEFINE_ITERATED_DOUBLE_FUNCTION(Minimum,"min",Min)
DEFINE_ITERATED_DOUBLE_FUNCTION(Maximum,"max",Max)

#define DEFINE_UNARY_DOUBLE_FUNCTION(NAME,TAG,OP)                  \
  DEFINE_FUNCTION(NAME,TAG)                                        \
  {                                                                \
    if (args.size()!=1)                                            \
      THROW(fatal_error,std::string(TAG)+" requires 1 argument."); \
    double arg0=ToType<double>(args[0]);                           \
    return ToString(OP(arg0));                                     \
  }                                                                \
  Term *NAME::Evaluate(const std::vector<Term*> &args) const       \
  {                                                                \
    ((TDouble*)args[0])->m_value=OP(((TDouble*)args[0])->m_value); \
    return args[0];                                                \
  }

DEFINE_UNARY_DOUBLE_FUNCTION(Logarithm,"log",log)
DEFINE_UNARY_DOUBLE_FUNCTION(Logarithm10,"log10",log10)
DEFINE_UNARY_DOUBLE_FUNCTION(Exponential,"exp",exp)
DEFINE_UNARY_DOUBLE_FUNCTION(Absolute_Value,"abs",dabs)
DEFINE_UNARY_DOUBLE_FUNCTION(Prefix,"sgn",Sign)
DEFINE_UNARY_DOUBLE_FUNCTION(Square,"sqr",sqr)
DEFINE_UNARY_DOUBLE_FUNCTION(Square_Root,"sqrt",sqrt)
DEFINE_UNARY_DOUBLE_FUNCTION(Sine,"sin",sin)
DEFINE_UNARY_DOUBLE_FUNCTION(Cosine,"cos",cos)
DEFINE_UNARY_DOUBLE_FUNCTION(Tangent,"tan",tan)
DEFINE_UNARY_DOUBLE_FUNCTION(Arc_Sine,"asin",asin)
DEFINE_UNARY_DOUBLE_FUNCTION(Arc_Cosine,"acos",acos)
DEFINE_UNARY_DOUBLE_FUNCTION(Arc_Tangent,"atan",atan)

#ifndef USING__double_only
DEFINE_FUNCTION(Vec4D_Part,"Part")
{
  if (args.size()!=2) THROW(fatal_error,"Part requires 2 arguments.");
  Vec4D arg0=ToType<Vec4D>(args[0]);
  int arg1=ToType<int>(args[1]);
  return ToString(arg0[arg1]);
}

Term *Vec4D_Part::Evaluate(const std::vector<Term*> &args) const
{
  ((TDouble*)args[1])->m_value=(((TVec4D*)args[0])->m_value)
    [(int)((TDouble*)args[1])->m_value];
  return args[1];
}

#define DEFINE_ONE_VECTOR_FUNCTION(NAME,TAG,OP)				\
  DEFINE_FUNCTION(NAME,TAG)						\
  {									\
    if (args.size()!=1)							\
      THROW(fatal_error,std::string(TAG)+" requires 1 argument.");      \
    Vec4D arg0=ToType<Vec4D>(args[0]);  			        \
    return ToString(arg0.OP());						\
  }									\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const		\
  {									\
    TDouble *res = new TDouble((((TVec4D*)args[0])->m_value).OP());	\
    p_interpreter->AddTerm(res);                                        \
    return res;								\
  }									\

DEFINE_ONE_VECTOR_FUNCTION(Vec4D_Abs2,"Abs2",Abs2)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_Mass,"Mass",Mass)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_PPerp,"PPerp",PPerp)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_PPerp2,"PPerp2",PPerp2)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_MPerp,"MPerp",MPerp)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_MPerp2,"MPerp2",MPerp2)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_Theta,"Theta",Theta)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_Eta,"Eta",Eta)
DEFINE_ONE_VECTOR_FUNCTION(Vec4D_Phi,"Phi",Phi)

#define DEFINE_TWO_VECTOR_FUNCTION(NAME,TAG,OP)				\
  DEFINE_FUNCTION(NAME,TAG)						\
  {									\
    if (args.size()!=2)							\
      THROW(fatal_error,"Operator requires 2 arguments.");		\
    Vec4D arg0=ToType<Vec4D>(args[0]);  			        \
    Vec4D arg1=ToType<Vec4D>(args[1]);	        		        \
    return ToString(arg0.OP(arg1));					\
  }									\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const		\
  {									\
    TDouble *res = new TDouble((((TVec4D*)args[0])->m_value).		\
      OP(((TVec4D*)args[1])->m_value));					\
    p_interpreter->AddTerm(res);                                        \
    return res;								\
  }									\

DEFINE_TWO_VECTOR_FUNCTION(Vec4D_PPerpR,"PPerpR",PPerp)
DEFINE_TWO_VECTOR_FUNCTION(Vec4D_ThetaR,"ThetaR",Theta)
DEFINE_TWO_VECTOR_FUNCTION(Vec4D_DEta,"DEta",DEta)
DEFINE_TWO_VECTOR_FUNCTION(Vec4D_DPhi,"DPhi",DPhi)

#define DEFINE_ITERATED_VECTOR_OPERATOR(NAME,TAG,OP)                 \
  DEFINE_FUNCTION(NAME,TAG)                                          \
  {                                                                  \
    if (args.size()<2)                                               \
      THROW(fatal_error,std::string(TAG)+" requires 2 arguments.");  \
    Vec4D arg0=ToType<Vec4D>(args[0]);				     \
    for (size_t i=1;i<args.size();++i)                               \
      arg0=arg0 OP ToType<Vec4D>(args[i]);			     \
    return ToString(arg0);                                           \
  }                                                                  \
  Term *NAME::Evaluate(const std::vector<Term*> &args) const         \
  {                                                                  \
    for (size_t i=1;i<args.size();++i)                               \
      ((TVec4D*)args[0])->m_value=((TVec4D*)args[0])->m_value OP     \
	((TVec4D*)args[i])->m_value;				     \
    return args[0];                                                  \
  }

DEFINE_ITERATED_VECTOR_OPERATOR(Vec4D_Plus,"Plus",+)
DEFINE_ITERATED_VECTOR_OPERATOR(Vec4D_Minus,"Minus",-)

#define DEFINE_TWO_VECTOR_OPERATOR(NAME,TAG,OP)				\
  DEFINE_FUNCTION(NAME,TAG)						\
  {									\
    if (args.size()!=2)							\
      THROW(fatal_error,"Operator requires 2 arguments.");		\
    Vec4D arg0=ToType<Vec4D>(args[0]);  			        \
    Vec4D arg1=ToType<Vec4D>(args[1]);	        		        \
    return ToString(arg0 OP arg1);					\
  }									\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const		\
  {									\
    TDouble *res = new TDouble((((TVec4D*)args[0])->m_value)		\
			       OP(((TVec4D*)args[1])->m_value));	\
    p_interpreter->AddTerm(res);                                        \
    return res;								\
  }									\

DEFINE_TWO_VECTOR_OPERATOR(Vec4D_Times,"Times",*)
#endif

Interpreter_Function::~Interpreter_Function() 
{
}

DEFINE_INTERPRETER_FUNCTION(Resolve_Bracket)
{
  if (expr.find("(")==std::string::npos ||
      expr.find(")")==std::string::npos) return expr;
  static int cnt=0;
  if (++cnt==100) 
    THROW(critical_error,"Bracket nesting deeper than 100 levels.");
  int open=0, take=1;
  size_t l=std::string::npos, r=std::string::npos;
  for (size_t i=0;i<expr.length();++i) {
    if (expr[i]=='(') {
      if (i>0 && expr[i-1]=='{') {
	--cnt;
	return expr;
      }
      ++open;
      if (l==std::string::npos) {
	Algebra_Interpreter::Function_Map::const_reverse_iterator fit=
	  p_interpreter->Functions().rbegin();
	for (;fit!=p_interpreter->Functions().rend();++fit) {
	  size_t pos=expr.rfind(fit->second->Tag(),i);
	  if (pos!=std::string::npos &&
	      pos+fit->second->Tag().length()==i) break;
	}
	if (fit==p_interpreter->Functions().rend()) {
	  l=i;
	  take=open;
	}
      }    
    }  
    else if (expr[i]==')') {
      if (open==take && l!=std::string::npos && 
	  r==std::string::npos) r=i;
      --open;
    }
  }
  if (open!=0) THROW(fatal_error,"Ambiguous bracket structure .");
  if (l==std::string::npos) {
    --cnt;
    return expr;
  }
  std::string left=expr.substr(0,l);
  std::string right=expr.substr(r+1);
  msg_Debugging()<<"Resolve_Bracket -> '"
		<<left<<"' '"<<expr.substr(l+1,r-l-1)<<"' '"<<right<<"'\n";
  std::string res=p_interpreter->
    Iterate(left+p_interpreter->
	    Iterate(expr.substr(l+1,r-l-1))+right);
  --cnt;
  return res;
}

DEFINE_INTERPRETER_FUNCTION(Extract_Leaf)
{
  if (expr.find("{")!=0 || expr.find("}")!=expr.length()-1) {
    Node<Function*> *leaf=p_interpreter->Leaf();
    Function *func=NULL;
#ifndef USING__double_only
    if (expr.find(",")!=std::string::npos)
      func = new Vector("("+expr+")",p_interpreter->TagReplacer());
    else if (expr.find("[")!=std::string::npos)
      func = new Vector(expr,p_interpreter->TagReplacer());
    else 
#endif
      func = new Number(expr,p_interpreter->TagReplacer());
    p_interpreter->AddLeaf(func);
    (*leaf)[0]=func;
    std::string value=expr;
    value=p_interpreter->TagReplacer()->ReplaceTags(value);
    return value;
  }
  size_t pos=expr.rfind(",");
  Node<Function*> *leaf=p_interpreter->Leaf(), *mother=--*leaf;
  if (mother!=NULL) 
    for (size_t i=0;i<(*mother)->size();++i)
      if ((*mother)()[i]==leaf) {
	delete (*mother)()[i];
	PTS add(ToType<PTS>(expr.substr(pos+1,expr.length()-pos-2)));
	(*mother)()[i]=dynamic_cast<Node<Function*>*>
	  ((Node<Function*>*)add);
	if ((*mother)()[i]==NULL) 
	  THROW(fatal_error,"Cannot recover node pointer.");
	p_interpreter->SetLeaf((*mother)()[i]);
      }
  std::string value=expr.substr(1,pos-1);
  value=p_interpreter->TagReplacer()->ReplaceTags(value);
  return value;
}

DEFINE_INTERPRETER_FUNCTION(Interprete_Function)
{
  if (expr.find("(")==std::string::npos ||
      expr.find(")")==std::string::npos) return expr;
  Function *func=NULL;
  size_t pos=std::string::npos, rem=std::string::npos;
  for (Algebra_Interpreter::Function_Map::const_reverse_iterator 
	 fit=p_interpreter->Functions().rbegin();
       fit!=p_interpreter->Functions().rend();++fit) {
    if ((pos=expr.rfind(fit->second->Tag()))!=std::string::npos &&
	pos<rem) {
      func=fit->second;
      rem=pos;
    }}
  if (func==NULL) return expr;
  pos=rem;
  size_t last=pos+func->Tag().length()+1, open=0;
  std::vector<std::string> args;
  size_t i=last-1;
  for (;i<expr.length();++i) {
    if (expr[i]=='(' || expr[i]=='{') ++open;
    if (expr[i]==')' || expr[i]=='}') --open;
    if (open==0 || (open==1 && expr[i]==',')) {
      args.push_back(expr.substr(last,i-last));
      last=i+1;
      if (open==0) break;
    }
  }
  std::string left=expr.substr(0,pos);
  std::string right=expr.substr(i+1);
#ifndef USING__Calc_only
  if (msg_LevelIsDebugging()) {
    std::string out=args[0];
    for (size_t j=1;j<args.size();++j) out+=","+args[j];
    msg_Debugging()<<"Interprete_Function -> '"<<left
		  <<"' '"<<func->Tag()<<"("<<out
		  <<")' '"<<right<<"'\n";
  }
#endif
  Node<Function*> *leaf = new Node<Function*>(func,true);
  for (size_t j=0;j<args.size();++j) {
    (*leaf)->push_back(new Node<Function*>(NULL,true));
    (*(*leaf)->back())<<leaf;
    args[j]=p_interpreter->Iterate(args[j]);
    p_interpreter->SetLeaf((*leaf)->back());
    args[j]=p_interpreter->Extractor()->Interprete(args[j]);
  }
  p_interpreter->SetLeaf(leaf);
  return p_interpreter->
    Iterate(left+"{"+func->Evaluate(args)+","+
	    ToString(PT(leaf))+"}"+right);
}

size_t FindBinaryPlus(const std::string &expr,const bool fwd,
		      size_t cpos=std::string::npos)
{
  if (cpos==std::string::npos && fwd) cpos=0;
  size_t pos(fwd?expr.find("+",cpos):expr.rfind("+",cpos));
  if (pos==std::string::npos || (pos==0 && !fwd)) return std::string::npos;
  if (pos==0) return FindBinaryPlus(expr,fwd,1);
  if (((expr[pos-1]=='e' || expr[pos-1]=='E') && pos>1 &&
       ((expr[pos-2]>=48 && expr[pos-2]<=57) || expr[pos-2]=='.')) ||
      expr[pos-1]==',' || expr[pos-1]=='(' || expr[pos-1]=='{') 
    return FindBinaryPlus(expr,fwd,fwd?pos+1:pos-1);
  return pos;  
}

size_t FindBinaryMinus(const std::string &expr,const bool fwd,
		       size_t cpos=std::string::npos)
{
  if (cpos==std::string::npos && fwd) cpos=0;
  size_t pos(fwd?expr.find("-",cpos):expr.rfind("-",cpos));
  if (pos==std::string::npos || (pos==0 && !fwd)) return std::string::npos;
  if (pos==0) return FindBinaryMinus(expr,fwd,1);
  if (((expr[pos-1]=='e' || expr[pos-1]=='E') && pos>1 &&
       ((expr[pos-2]>=48 && expr[pos-2]<=57) || expr[pos-2]=='.')) ||
      expr[pos-1]==',' || expr[pos-1]=='(' || expr[pos-1]=='{') 
    return FindBinaryMinus(expr,fwd,fwd?pos+1:pos-1);
  return pos;  
}

size_t FindUnaryNot(const std::string &expr,const bool fwd,
		    size_t cpos=std::string::npos)
{
  if (cpos==std::string::npos && fwd) cpos=0;
  size_t pos(fwd?expr.find("!",cpos):expr.rfind("!",cpos));
  if (pos==std::string::npos || pos+1>=expr.length()) 
    return std::string::npos;
  if (expr[pos+1]=='=') return FindUnaryNot(expr,fwd,fwd?pos+1:pos-1);
  return pos;  
}

DEFINE_INTERPRETER_FUNCTION(Interprete_Binary)
{
  if (expr.find("(")!=std::string::npos ||
      expr.find(")")!=std::string::npos) return expr;
  Operator *op=NULL;
  size_t pos=std::string::npos;
  for (Algebra_Interpreter::Operator_Map::const_reverse_iterator 
	 oit=p_interpreter->Operators().rbegin();
       oit!=p_interpreter->Operators().rend();++oit) {
    if (oit->second->Tag()=="!") pos=FindUnaryNot(expr,true);
    else if (oit->second->Tag()=="-") pos=FindBinaryMinus(expr,true);
    else if (oit->second->Tag()=="+") pos=FindBinaryPlus(expr,true);
    else pos=expr.find(oit->second->Tag());
    if (pos!=std::string::npos) {
      if (oit->second->Binary() && pos!=0) {
	op=oit->second;
	break;
      }
      else {
	return expr;
      }
    }
  }
  if (op==NULL) return expr;
  std::string lrstr, lstr=expr.substr(0,pos);
  size_t lfpos=0;
  for (Algebra_Interpreter::Operator_Map::const_reverse_iterator 
	 oit=p_interpreter->Operators().rbegin();
       oit!=p_interpreter->Operators().rend();++oit) {
    size_t tlfpos=std::string::npos;
    if (oit->second->Tag()=="-") tlfpos=FindBinaryMinus(lstr,false);
    else if (oit->second->Tag()=="+") tlfpos=FindBinaryPlus(lstr,false);
    else tlfpos=lstr.rfind(oit->second->Tag());
    if (tlfpos!=std::string::npos) 
      lfpos=Max(lfpos,tlfpos+oit->second->Tag().length());
  }
  lrstr=lstr.substr(0,lfpos);
  lstr=lstr.substr(lfpos);
  std::string rrstr, rstr=expr.substr(pos+op->Tag().length());
  size_t rfpos=rstr.length();
  for (Algebra_Interpreter::Operator_Map::const_reverse_iterator 
	 oit=p_interpreter->Operators().rbegin();
       oit!=p_interpreter->Operators().rend();++oit) {
    size_t trfpos=std::string::npos;
    if (oit->second->Tag()=="-") trfpos=FindBinaryMinus(rstr,true);
    else if (oit->second->Tag()=="+") trfpos=FindBinaryPlus(rstr,true);
    else trfpos=rstr.find(oit->second->Tag());
    if (trfpos!=std::string::npos) rfpos=Min(rfpos,trfpos);
  }
  rrstr=rstr.substr(rfpos);
  rstr=rstr.substr(0,rfpos);
  Node<Function*> *leaf = new Node<Function*>(op,true);
  std::vector<std::string> args(2);
  (*leaf)->push_back(new Node<Function*>(NULL,true));
  (*(*leaf)->back())<<leaf;
  args[0]=p_interpreter->Iterate(lstr);
  p_interpreter->SetLeaf((*leaf)->back());
  args[0]=p_interpreter->Extractor()->Interprete(args[0]);
  (*leaf)->push_back(new Node<Function*>(NULL,true));
  (*(*leaf)->back())<<leaf;
  args[1]=p_interpreter->Iterate(rstr);
  p_interpreter->SetLeaf((*leaf)->back());
  args[1]=p_interpreter->Extractor()->Interprete(args[1]);
  p_interpreter->SetLeaf(leaf);
  msg_Debugging()<<"Interprete_Binary -> '"
	    <<lrstr<<"' '"<<args[0]<<"' '"<<op->Tag()
	    <<"' '"<<args[1]<<"' '"<<rrstr<<"'\n";
  return p_interpreter->
    Iterate(lrstr+"{"+op->Evaluate(args)+","+
	    ToString(PT(leaf))+"}"+rrstr);
}

DEFINE_INTERPRETER_FUNCTION(Interprete_Unary)
{
  if (expr.find("(")!=std::string::npos ||
      expr.find(")")!=std::string::npos) return expr;
  Operator *op=NULL;
  size_t pos=std::string::npos;
  for (Algebra_Interpreter::Operator_Map::const_reverse_iterator 
	 oit=p_interpreter->Operators().rbegin();
       oit!=p_interpreter->Operators().rend();++oit) {
    if (oit->second->Tag()=="!") pos=FindUnaryNot(expr,false);
    else pos=expr.rfind(oit->second->Tag());
    if (pos!=std::string::npos) {
      if (!oit->second->Binary()) {
	op=oit->second;
	break;
      }
      else {
	return expr;
      }
    }
  }
  if (op==NULL) return expr;
  std::string lrstr=expr.substr(0,pos);
  std::string rrstr, rstr=expr.substr(pos+op->Tag().length()), mrstr=rstr;
  bool negative=rstr.length()>0?rstr[0]=='-':false;
  if (negative) rstr=rstr.substr(1);
  size_t rfpos=rstr.length();
  for (Algebra_Interpreter::Operator_Map::const_reverse_iterator 
	 oit=p_interpreter->Operators().rbegin();
       oit!=p_interpreter->Operators().rend();++oit) {
    size_t trfpos=std::string::npos;
    if (oit->second->Tag()=="-") trfpos=FindBinaryMinus(rstr,true);
    else if (oit->second->Tag()=="+") trfpos=FindBinaryPlus(rstr,true);
    else trfpos=rstr.find(oit->second->Tag());
    if (trfpos!=std::string::npos) rfpos=Min(rfpos,trfpos);
  }
  if (negative) {
    rstr=mrstr;
    ++rfpos;
  }
  rrstr=rstr.substr(rfpos);
  rstr=rstr.substr(0,rfpos);
  Node<Function*> *leaf = new Node<Function*>(op,true);
  std::vector<std::string> args(1);
  (*leaf)->push_back(new Node<Function*>(NULL,true));
  (*(*leaf)->back())<<leaf;
  args[0]=p_interpreter->Iterate(rstr);
  p_interpreter->SetLeaf((*leaf)->back());
  args[0]=p_interpreter->Extractor()->Interprete(args[0]);
  p_interpreter->SetLeaf(leaf);
  msg_Debugging()<<"Interprete_Unary -> '"
		<<lrstr<<"' '"<<op->Tag()<<"' '"<<args[0]<<"' '"<<rrstr<<"'\n";
  return p_interpreter->
    Iterate(lrstr+"{"+op->Evaluate(args)+","+
	    ToString(PT(leaf))+"}"+rrstr);
}

Algebra_Interpreter::Algebra_Interpreter(const bool standard):
  p_replacer(this), p_root(NULL), p_leaf(NULL)
{
  p_extractor = new Extract_Leaf(this);
  m_interpreters[0] = new Interprete_Function(this);
  m_interpreters[1] = new Resolve_Bracket(this);
  m_interpreters[2] = new Interprete_Binary(this);
  m_interpreters[3] = new Interprete_Unary(this);
  if (!standard) return;
  m_tags["M_PI"]=ToString(M_PI);
  m_tags["M_E"]=ToString(exp(1.0));
  AddOperator(new Binary_Plus());
  AddOperator(new Binary_Minus());
  AddOperator(new Binary_Times());
  AddOperator(new Binary_Divide());
  AddOperator(new Binary_Equal());
  AddOperator(new Binary_Not_Equal());
  AddOperator(new Binary_Less());
  AddOperator(new Binary_Greater());
  AddOperator(new Binary_Less_Equal());
  AddOperator(new Binary_Greater_Equal());
  AddOperator(new Binary_Modulus());
  AddOperator(new Binary_Shift_Left());
  AddOperator(new Binary_Shift_Right());
  AddOperator(new Binary_Logical_And());
  AddOperator(new Binary_Logical_Or());
  AddOperator(new Unary_Not());
  AddFunction(new Power());
  AddFunction(new Logarithm());
  AddFunction(new Logarithm10());
  AddFunction(new Exponential());
  AddFunction(new Absolute_Value());
  AddFunction(new Prefix());
  AddFunction(new Square());
  AddFunction(new Square_Root());
  AddFunction(new Sine());
  AddFunction(new Cosine());
  AddFunction(new Tangent());
  AddFunction(new Arc_Sine());
  AddFunction(new Arc_Cosine());
  AddFunction(new Arc_Tangent());
  AddFunction(new Minimum());
  AddFunction(new Maximum());
#ifndef USING__double_only
  AddFunction(new Vec4D_Part());
  AddFunction(new Vec4D_Abs2());
  AddFunction(new Vec4D_Mass());
  AddFunction(new Vec4D_PPerp());
  AddFunction(new Vec4D_PPerp2());
  AddFunction(new Vec4D_MPerp());
  AddFunction(new Vec4D_MPerp2());
  AddFunction(new Vec4D_Theta());
  AddFunction(new Vec4D_Eta());
  AddFunction(new Vec4D_Phi());
  AddFunction(new Vec4D_PPerpR());
  AddFunction(new Vec4D_ThetaR());
  AddFunction(new Vec4D_DEta());
  AddFunction(new Vec4D_DPhi());
  AddFunction(new Vec4D_Plus());
  AddFunction(new Vec4D_Minus());
  AddFunction(new Vec4D_Times());
#endif
}

Algebra_Interpreter::~Algebra_Interpreter()
{
  while (!m_operators.empty()) {
    delete m_operators.begin()->second;
    m_operators.erase(m_operators.begin());
  }
  while (!m_functions.empty()) {
    delete m_functions.begin()->second;
    m_functions.erase(m_functions.begin());
  }
  if (p_root!=NULL) delete p_root;
  while (m_leafs.size()>0) {
    delete m_leafs.begin()->second;
    m_leafs.erase(m_leafs.begin());
  }
  while (m_terms.size()>0) {
    delete *m_terms.begin();
    m_terms.erase(m_terms.begin());
  }
  while (m_interpreters.size()>0) {
    delete m_interpreters.begin()->second;
    m_interpreters.erase(m_interpreters.begin());
  }
  delete p_extractor;
}

std::string Algebra_Interpreter::Interprete(const std::string &expr)
{
#ifdef DEBUG__Interpreter
  msg_Debugging()<<METHOD<<"("<<expr<<") {\n";
#endif
  if (p_root!=NULL) delete p_root;
  p_root=p_leaf=NULL;
  while (m_leafs.size()>0) {
    delete m_leafs.begin()->second;
    m_leafs.erase(m_leafs.begin());
  }
  if (expr.length()==0) return "0";
  std::string res=expr;
  KillBlanks(res);
  std::string result=Iterate(res);
  size_t pos=result.find(",");
  if (pos==std::string::npos) {
    p_root = p_leaf = new Node<Function*>(NULL,false);
#ifndef USING__double_only
    if (result.find(",")!=std::string::npos)
      (*p_leaf)[0] = new Vector("("+result+")",p_replacer);
    else if (expr.find("[")!=std::string::npos)
      (*p_leaf)[0] = new Vector(result,p_replacer);
    else 
#endif
      (*p_leaf)[0] = new Number(result,p_replacer);
    AddLeaf((*p_leaf)[0]);
    result=p_replacer->ReplaceTags(result);
#ifdef DEBUG__Interpreter
    msg_Debugging()<<"} -> "<<result<<std::endl;
#endif
    return result;
  }
#ifdef DEBUG__Interpreter
  msg_Debugging()<<"} -> "<<result.substr(1,pos-1)<<std::endl;
#endif
  p_root = p_leaf;
  return result.substr(1,pos-1);
}

Term *Algebra_Interpreter::Calculate()
{
  if (p_root==NULL) return new TDouble(0.0);
  while (m_terms.size()>0) {
    delete *m_terms.begin();
    m_terms.erase(m_terms.begin());
  }
  return Iterate(p_root);
}

Term *Algebra_Interpreter::Iterate(Node<Function*> *const node)
{
  if (node->operator->()==NULL) {
    std::vector<Term*> args;
    return (*node)[0]->Evaluate(args);
  }
  std::vector<Term*> args((*node)->size());
  for (size_t i=0;i<(*node)->size();++i) {
    args[i]=Iterate((*node)()[i]);
  }
  return (*node)[0]->Evaluate(args);
}

void Algebra_Interpreter::AddFunction(Function *const f)
{
  m_functions.insert(Function_Pair(f->Tag(),f)); 
  f->SetInterpreter(this);
}

void Algebra_Interpreter::AddOperator(Operator *const b)
{ 
  m_operators.insert(Operator_Pair(b->Priority(),b)); 
  b->SetInterpreter(this);
}

void Algebra_Interpreter::AddLeaf(Function *const f)
{
  m_leafs.insert(Function_Pair(ToString(f),f)); 
}

void Algebra_Interpreter::AddTerm(Term *const t)
{
  m_terms.insert(t); 
}

std::string Algebra_Interpreter::Iterate(const std::string &expr)
{
  static size_t depth=0;
  if (++depth>1000) THROW(critical_error,"Max depth reached.");
  msg_Indent();
  std::string res=expr;
  Interpreter_Map::const_iterator iit=m_interpreters.begin();
  for (;iit!=m_interpreters.end();++iit) 
    res=iit->second->Interprete(res);
  --depth;
  return res;
}

void Algebra_Interpreter::AddTag(const std::string &tag,
				   const std::string &value)
{
  m_tags[tag]=value;
}

void Algebra_Interpreter::SetTags(const String_Map &tags)
{
  m_tags=tags;
}

std::string Algebra_Interpreter::ReplaceTags(std::string &expr) const
{
  size_t pos=std::string::npos;
  for (String_Map::const_reverse_iterator sit=m_tags.rbegin();
       sit!=m_tags.rend();++sit) {
    if ((pos=expr.find(sit->first))!=std::string::npos) 
      return ReplaceTags(expr.replace(pos,sit->first.length(),
				      sit->second));}
  return expr;
}

Term *Algebra_Interpreter::ReplaceTags(Term *expr) const
{
  size_t pos=std::string::npos;
  for (String_Map::const_reverse_iterator sit=m_tags.rbegin();
       sit!=m_tags.rend();++sit) {
    if ((pos=expr->m_tag.find(sit->first))!=std::string::npos) {
#ifndef USING__double_only
      if (sit->second[0]=='(' && sit->second[sit->second.length()-1]==')') 
	((TVec4D*)expr)->m_value=ToType<Vec4D>(sit->second);      
      else 
#endif
	((TDouble*)expr)->m_value=ToType<double>(sit->second);      
      return expr;
    }
  }
  return expr;
}

void Algebra_Interpreter::PrintNode(Node<Function*> *const node) const
{
  msg_Info()<<"("<<node<<") ["<<typeid(*(*node)[0]).name()<<"] '"
	    <<((*node)[0]!=NULL?(*node)[0]->Tag():"<NULL>")<<"' {\n";
  {
    msg_Indent();
    if (node->operator->()!=NULL) 
      for (size_t i=0;i<(*node)().size();++i) PrintNode((*node)()[i]);
    else msg_Info()<<"<NULL>\n";
  }
  msg_Info()<<"}\n";
}

Tag_Replacer::~Tag_Replacer() {}

std::string &Tag_Replacer::KillBlanks(std::string& expr) const
{
  for (size_t i=0;i<expr.length();++i) 
    if (expr[i]==32 || expr[i]==9) expr.replace(i--,1,"");
  return expr;
}

std::string Tag_Replacer::ReplaceTags(std::string &expr) const
{
  return expr;
}

Term *Tag_Replacer::ReplaceTags(Term *term) const
{
  return term;
}
