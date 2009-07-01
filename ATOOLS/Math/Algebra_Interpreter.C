#ifdef USING__Calc_only
#include "Algebra_Interpreter.H"
#include "Tools.H"
#else
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
#endif
#include <typeinfo>

using namespace ATOOLS;

#define PTS long unsigned int
#define PT(ARG) (PTS)(ARG)

Function::Function(const std::string &tag): 
  m_tag(tag), p_interpreter(NULL) {}

Function::~Function() {}

Term *Function::Evaluate(const std::vector<Term*> &args) const
{
  THROW(fatal_error,"No evaluation rule.");
}

class Single_Term: public Function {
private:

  Tag_Replacer *p_replacer;

  bool m_replace, m_sign;

  mutable Term *p_value;

public:

  Single_Term(const std::string &tag,Tag_Replacer *const replacer);

  ~Single_Term();

  Term *Evaluate(const std::vector<Term*> &args) const;

};// end of class Single_Term

Single_Term::Single_Term(const std::string &tag,Tag_Replacer *const replacer): 
  Function(tag), p_replacer(replacer), m_replace(false), 
  m_sign(0), p_value(NULL)
{
  std::string value=tag;
  if (tag[0]=='-') {
    m_sign=1;
    value=value.substr(1);
  }
  std::string stag(value);
  p_replacer->ReplaceTags(value);
  if (stag!=value) m_replace=true;
  p_value = Term::New(value);
  p_value->SetTag(stag);
  if (m_replace) p_replacer->AssignId(p_value);
}

Single_Term::~Single_Term()
{
  if (p_value) delete p_value;
}

Term *Single_Term::Evaluate(const std::vector<Term*> &args) const
{
  if (args.size()!=0) THROW(fatal_error,"Single_Term requires no argument.");
  if (m_replace) p_replacer->ReplaceTags(p_value);
  if (m_sign) *p_value=-*p_value;
  return p_value;
}

Operator::~Operator() {}

#define DEFINE_BINARY_TERM_OPERATOR(NAME,TAG,OP,PRIORITY)\
  DEFINE_BINARY_OPERATOR(NAME,TAG,PRIORITY)\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const\
  {\
    Term *res = *args[0] OP *args[1];\
    p_interpreter->AddTerm(res);\
    return res;\
  }

DEFINE_BINARY_TERM_OPERATOR(Binary_Plus,"+",+,12)
DEFINE_BINARY_TERM_OPERATOR(Binary_Minus,"-",-,12)
DEFINE_BINARY_TERM_OPERATOR(Binary_Times,"*",*,13)
DEFINE_BINARY_TERM_OPERATOR(Binary_Divide,"/",/,13)
DEFINE_BINARY_TERM_OPERATOR(Binary_Modulus,"%",%,13)
DEFINE_BINARY_TERM_OPERATOR(Binary_Shift_Left,"<<",<<,11)
DEFINE_BINARY_TERM_OPERATOR(Binary_Shift_Right,">>",>>,11)
DEFINE_BINARY_TERM_OPERATOR(Binary_Logical_And,"&&",&&,5)
DEFINE_BINARY_TERM_OPERATOR(Binary_Logical_Or,"||",||,4)
DEFINE_BINARY_TERM_OPERATOR(Bitwise_And,"&",&,8)
DEFINE_BINARY_TERM_OPERATOR(Bitwise_XOr,"^",^,7)
DEFINE_BINARY_TERM_OPERATOR(Bitwise_Or,"|",|,6)

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

#define DEFINE_BINARY_SORTABLE_TERM_OPERATOR(NAME,TAG,OP,PRIORITY)\
  DEFINE_BINARY_OPERATOR(NAME,TAG,PRIORITY)\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const\
  {\
    Term *res = *args[0] OP *args[1];\
    p_interpreter->AddTerm(res);\
    return res;\
  }

DEFINE_BINARY_SORTABLE_TERM_OPERATOR(Binary_Equal,"==",==,9)
DEFINE_BINARY_SORTABLE_TERM_OPERATOR(Binary_Not_Equal,"!=",!=,9)
DEFINE_BINARY_SORTABLE_TERM_OPERATOR(Binary_Less,"<",<,10)
DEFINE_BINARY_SORTABLE_TERM_OPERATOR(Binary_Greater,">",>,10)
DEFINE_BINARY_SORTABLE_TERM_OPERATOR(Binary_Less_Equal,"<=",<=,10)
DEFINE_BINARY_SORTABLE_TERM_OPERATOR(Binary_Greater_Equal,">=",>=,10)

DEFINE_UNARY_OPERATOR(Unary_Not,"!",14)
Term *Unary_Not::Evaluate(const std::vector<Term*> &args) const
{
  Term *res = !*args[0];
  p_interpreter->AddTerm(res);
  return res;
}

#define DEFINE_BINARY_TERM_FUNCTION(NAME,TAG,OP)\
  DEFINE_FUNCTION(NAME,TAG)\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const\
  {\
    Term *res(OP(*args[0],*args[1]));\
    p_interpreter->AddTerm(res);\
    return res;\
  }

DEFINE_BINARY_TERM_FUNCTION(Power,"pow",TPow)

#define DEFINE_ITERATED_TERM_FUNCTION(NAME,TAG,OP)\
  DEFINE_FUNCTION(NAME,TAG)\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const\
  {\
    Term *res(args[0]);\
    for (size_t i=1;i<args.size();++i) {\
      res = OP(*res,*args[i]);\
      p_interpreter->AddTerm(res);\
    }\
    return res;\
  }

DEFINE_ITERATED_TERM_FUNCTION(Minimum,"min",TMin)
DEFINE_ITERATED_TERM_FUNCTION(Maximum,"max",TMax)

#define DEFINE_UNARY_TERM_FUNCTION(NAME,TAG,OP)\
  DEFINE_FUNCTION(NAME,TAG)\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const\
  {\
    Term *res(OP(*args[0]));\
    p_interpreter->AddTerm(res);\
    return res;\
  }

DEFINE_UNARY_TERM_FUNCTION(Logarithm,"log",TLog)
DEFINE_UNARY_TERM_FUNCTION(Logarithm10,"log10",TLog10)
DEFINE_UNARY_TERM_FUNCTION(Exponential,"exp",TExp)
DEFINE_UNARY_TERM_FUNCTION(Absolute_Value,"abs",TAbs)
DEFINE_UNARY_TERM_FUNCTION(Prefix,"sgn",TSgn)
DEFINE_UNARY_TERM_FUNCTION(Square,"sqr",TSqr)
DEFINE_UNARY_TERM_FUNCTION(Square_Root,"sqrt",TSqrt)
DEFINE_UNARY_TERM_FUNCTION(Sine,"sin",TSin)
DEFINE_UNARY_TERM_FUNCTION(Cosine,"cos",TCos)
DEFINE_UNARY_TERM_FUNCTION(Tangent,"tan",TTan)
DEFINE_UNARY_TERM_FUNCTION(Arc_Sine,"asin",TASin)
DEFINE_UNARY_TERM_FUNCTION(Arc_Cosine,"acos",TACos)
DEFINE_UNARY_TERM_FUNCTION(Arc_Tangent,"atan",TATan)

#define DEFINE_ONE_TERM_FUNCTION(NAME,TAG,OP)\
  DEFINE_FUNCTION(NAME,TAG)\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const\
  {\
    Term *res = args[0]->OP();\
    p_interpreter->AddTerm(res);\
    return res;\
  }

DEFINE_ONE_TERM_FUNCTION(Real,"Real",Real)
DEFINE_ONE_TERM_FUNCTION(Imag,"Imag",Imag)
DEFINE_ONE_TERM_FUNCTION(Conj,"Conj",Conj)

DEFINE_FUNCTION(Vec4D_Vec4D,"Vec4D")
Term *Vec4D_Vec4D::Evaluate(const std::vector<Term*> &args) const
{
  Term *res = TVec4D(*args[0],*args[1],*args[2],*args[3]);
  p_interpreter->AddTerm(res);
  return res;
}

DEFINE_ONE_TERM_FUNCTION(Vec4D_Perp,"Perp",Perp)
DEFINE_ONE_TERM_FUNCTION(Vec4D_Plus,"Plus",Plus)
DEFINE_ONE_TERM_FUNCTION(Vec4D_Minus,"Minus",Minus)
DEFINE_ONE_TERM_FUNCTION(Vec4D_PPlus,"PPlus",PPlus)
DEFINE_ONE_TERM_FUNCTION(Vec4D_PMinus,"PMinus",PMinus)
DEFINE_ONE_TERM_FUNCTION(Vec4D_Abs2,"Abs2",Abs2)
DEFINE_ONE_TERM_FUNCTION(Vec4D_Mass,"Mass",Mass)
DEFINE_ONE_TERM_FUNCTION(Vec4D_PPerp,"PPerp",PPerp)
DEFINE_ONE_TERM_FUNCTION(Vec4D_PPerp2,"PPerp2",PPerp2)
DEFINE_ONE_TERM_FUNCTION(Vec4D_MPerp,"MPerp",MPerp)
DEFINE_ONE_TERM_FUNCTION(Vec4D_MPerp2,"MPerp2",MPerp2)
DEFINE_ONE_TERM_FUNCTION(Vec4D_Theta,"Theta",Theta)
DEFINE_ONE_TERM_FUNCTION(Vec4D_Eta,"Eta",Eta)
DEFINE_ONE_TERM_FUNCTION(Vec4D_Phi,"Phi",Phi)

#define DEFINE_TWO_TERM_FUNCTION(NAME,TAG,OP)\
  DEFINE_FUNCTION(NAME,TAG)\
  Term *NAME::Evaluate(const std::vector<Term*> &args) const\
  {\
    Term *res = args[0]->OP(*args[1]);\
    p_interpreter->AddTerm(res);\
    return res;\
  }

DEFINE_TWO_TERM_FUNCTION(Vec4D_Comp,"Comp",Comp)
DEFINE_TWO_TERM_FUNCTION(Vec4D_PPerpR,"PPerpR",PPerp)
DEFINE_TWO_TERM_FUNCTION(Vec4D_ThetaR,"ThetaR",Theta)
DEFINE_TWO_TERM_FUNCTION(Vec4D_DEta,"DEta",DEta)
DEFINE_TWO_TERM_FUNCTION(Vec4D_DPhi,"DPhi",DPhi)

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
#ifdef DEBUG__Interpreter
  msg_IODebugging()<<"Resolve_Bracket -> '"
		<<left<<"' '"<<expr.substr(l+1,r-l-1)<<"' '"<<right<<"'\n";
#endif
  std::string mid(p_interpreter->Iterate(expr.substr(l+1,r-l-1)));
  std::string res=p_interpreter->Iterate(left+mid+right);
  --cnt;
  return res;
}

DEFINE_INTERPRETER_FUNCTION(Extract_Leaf)
{
  if (expr.find("{")!=0 || expr.find("}")!=expr.length()-1) {
    Node<Function*> *leaf=p_interpreter->Leaf();
    Function *func=NULL;
    std::string value(expr);
    if (expr.find(',')!=std::string::npos) value="("+expr+")";
    func = new Single_Term(value,p_interpreter->TagReplacer());
    p_interpreter->AddLeaf(func);
    (*leaf)[0]=func;
    value=p_interpreter->TagReplacer()->ReplaceTags(value);
    return "{"+ToString(PT(leaf))+"}";
  }
  Node<Function*> *leaf=p_interpreter->Leaf(), *mother=--*leaf;
  size_t pos(expr.rfind('{')); 
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
  return expr.substr(pos+1,expr.length()-pos-2);
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
#ifdef DEBUG__Interpreter
  if (msg_LevelIsIODebugging()) {
    std::string out=args[0];
    for (size_t j=1;j<args.size();++j) out+=","+args[j];
    msg_IODebugging()<<"Interprete_Function -> '"<<left
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
    Iterate(left+"{"+ToString(PT(leaf))+"}"+right);
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
#ifdef DEBUG__Interpreter
  msg_IODebugging()<<"Interprete_Binary -> '"
	    <<lrstr<<"' '"<<args[0]<<"' '"<<op->Tag()
	    <<"' '"<<args[1]<<"' '"<<rrstr<<"'\n";
#endif
  return p_interpreter->
    Iterate(lrstr+"{"+ToString(PT(leaf))+"}"+rrstr);
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
#ifdef DEBUG__Interpreter
  msg_IODebugging()<<"Interprete_Unary -> '"
		<<lrstr<<"' '"<<op->Tag()<<"' '"<<args[0]<<"' '"<<rrstr<<"'\n";
#endif
  return p_interpreter->
    Iterate(lrstr+"{"+ToString(PT(leaf))+"}"+rrstr);
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
  AddFunction(new Real());
  AddFunction(new Imag());
  AddFunction(new Conj());
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
  AddOperator(new Bitwise_And());
  AddOperator(new Bitwise_XOr());
  AddOperator(new Bitwise_Or());
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
  AddFunction(new Vec4D_Vec4D());
  AddFunction(new Vec4D_Comp());
  AddFunction(new Vec4D_Perp());
  AddFunction(new Vec4D_Plus());
  AddFunction(new Vec4D_Minus());
  AddFunction(new Vec4D_PPlus());
  AddFunction(new Vec4D_PMinus());
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
    delete m_terms.back();
    m_terms.pop_back();
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
  msg_IODebugging()<<METHOD<<"("<<expr<<") {\n";
#endif
  m_argvs.clear();
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
  if (result==res) {
    p_root = p_leaf = new Node<Function*>(NULL,false);
    (*p_leaf)[0] = new Single_Term(result,p_replacer);
    AddLeaf((*p_leaf)[0]);
    AddArgs(p_root);
    result=p_replacer->ReplaceTags(result);
#ifdef DEBUG__Interpreter
    msg_IODebugging()<<"} -> "<<result<<std::endl;
#endif
    return result;
  }
#ifdef DEBUG__Interpreter
  msg_IODebugging()<<"} -> "<<result<<std::endl;
#endif
  p_root = p_leaf;
  AddArgs(p_root);
  result=ToString(*Calculate());
  if (msg_LevelIsIODebugging()) {
    msg_IODebugging()<<METHOD<<"("<<expr<<"): {\n";
    {
      msg_Indent();
      PrintEquation();
    }
    msg_IODebugging()<<"} -> "<<result<<"\n";
  }
  return result;
}

Term *Algebra_Interpreter::Calculate()
{
  for (Term_Vector::const_iterator
	 tit(m_terms.begin());tit!=m_terms.end();++tit)
    (*tit)->Delete();
  m_terms.clear();
  if (p_root==NULL) {
    Term *res(Term::New(std::string("0.0")));
    m_terms.push_back(res);
    return res;
  }
  size_t n(0);
  return Iterate(p_root,n);
}

void Algebra_Interpreter::AddArgs(Node<Function*> *const node)
{
  if (node->operator->()==NULL) {
    m_argvs.push_back(Term_Vector(0));
    return;
  }
  m_argvs.push_back(Term_Vector((*node)->size()));
  for (size_t i=0;i<(*node)->size();++i)
    AddArgs((*node)()[i]);
}

Term *Algebra_Interpreter::Iterate
(Node<Function*> *const node,size_t &n)
{
  Term_Vector &args(m_argvs[n++]);
  if (node->operator->()==NULL)
    return (*node)[0]->Evaluate(args);
  for (size_t i=0;i<(*node)->size();++i)
    args[i]=Iterate((*node)()[i],n);
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
  m_terms.push_back(t); 
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
       sit!=m_tags.rend();++sit)
    if ((pos=expr.find(sit->first))!=std::string::npos)
      return ReplaceTags(expr.replace(pos,sit->first.length(),
				      sit->second));
  return expr;
}

Term *Algebra_Interpreter::ReplaceTags(Term *expr) const
{
  std::string cexpr(expr->Tag());
  expr->Set(ReplaceTags(cexpr));
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

void Tag_Replacer::AssignId(Term *term)
{
}
