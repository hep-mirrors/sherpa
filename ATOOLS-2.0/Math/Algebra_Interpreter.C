#include "Algebra_Interpreter.H"

#include "MyStrStream.H"
#include "MathTools.H"
#include "Exception.H"
#include "Message.H"

using namespace ATOOLS;

Function::Function(const std::string &tag): m_tag(tag) {}

Function::~Function() {}

std::string Function::Evaluate(const std::vector<std::string> &args) const
{
  THROW(fatal_error,"No evaluation rule.");
}

Operator::~Operator() {}

DEFINE_BINARY_OPERATOR(Binary_Plus,"+",12)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return ToString(arg0+arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Minus,"-",12)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return ToString(arg0-arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Times,"*",13)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return ToString(arg0*arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Divide,"/",13)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return ToString(arg0/arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Modulus,"%",13)
{
  long int arg0=ToType<long int>(args[0]);
  long int arg1=ToType<long int>(args[1]);
  return ToString(arg0%arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Shift_Left,"<<",11)
{
  long int arg0=ToType<long int>(args[0]);
  long int arg1=ToType<long int>(args[1]);
  return ToString(arg0<<arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Shift_Right,">>",11)
{
  long int arg0=ToType<long int>(args[0]);
  long int arg1=ToType<long int>(args[1]);
  return ToString(arg0>>arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Equal,"==",9)
{
  return args[0]==args[1]?"1":"0";
}

DEFINE_BINARY_OPERATOR(Binary_Not_Equal,"!=",9)
{
  return args[0]!=args[1]?"1":"0";
}

DEFINE_BINARY_OPERATOR(Binary_Less,"<",10)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return arg0<arg1?"1":"0";
}

DEFINE_BINARY_OPERATOR(Binary_Greater,">",10)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return arg0>arg1?"1":"0";
}

DEFINE_BINARY_OPERATOR(Binary_Less_Equal,"<=",10)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return arg0<=arg1?"1":"0";
}

DEFINE_BINARY_OPERATOR(Binary_Greater_Equal,">=",10)
{
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return arg0>=arg1?"1":"0";
}

DEFINE_BINARY_OPERATOR(Binary_Logical_And,"&&",5)
{
  int arg0=ToType<int>(args[0]);
  int arg1=ToType<int>(args[1]);
  return ToString(arg0&&arg1);
}

DEFINE_BINARY_OPERATOR(Binary_Logical_Or,"||",4)
{
  int arg0=ToType<int>(args[0]);
  int arg1=ToType<int>(args[1]);
  return ToString(arg0||arg1);
}

DEFINE_UNARY_OPERATOR(Unary_Not,"!",14)
{
  double arg0=ToType<int>(args[0]);
  return ToString(!arg0);
}

DEFINE_FUNCTION(Power,"pow")
{
  if (args.size()!=2) THROW(fatal_error,"Pow requires 2 arguments.")
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return ToString(pow(arg0,arg1));
}

DEFINE_FUNCTION(Logarithm,"log")
{
  if (args.size()!=1) THROW(fatal_error,"Log requires 1 argument.")
  double arg0=ToType<double>(args[0]);
  return ToString(log(arg0));
}

DEFINE_FUNCTION(Logarithm10,"log10")
{
  if (args.size()!=1) THROW(fatal_error,"Log10 requires 1 argument.")
  double arg0=ToType<double>(args[0]);
  return ToString(log10(arg0));
}

DEFINE_FUNCTION(Exponential,"exp")
{
  if (args.size()!=1) THROW(fatal_error,"Exp requires 1 argument.")
  double arg0=ToType<double>(args[0]);
  return ToString(exp(arg0));
}

DEFINE_FUNCTION(Absolute_Value,"abs")
{
  if (args.size()!=1) THROW(fatal_error,"Abs requires 1 argument.")
  double arg0=ToType<double>(args[0]);
  return ToString(dabs(arg0));
}

DEFINE_FUNCTION(Square,"sqr")
{
  if (args.size()!=1) THROW(fatal_error,"Sqr requires 1 argument.")
  double arg0=ToType<double>(args[0]);
  return ToString(sqr(arg0));
}

DEFINE_FUNCTION(Square_Root,"sqrt")
{
  if (args.size()!=1) THROW(fatal_error,"Sqrt requires 1 argument.")
  double arg0=ToType<double>(args[0]);
  return ToString(sqrt(arg0));
}

DEFINE_FUNCTION(Minimum,"min")
{
  if (args.size()!=2) THROW(fatal_error,"Min requires 2 arguments.")
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return ToString(Min(arg0,arg1));
}

DEFINE_FUNCTION(Maximum,"max")
{
  if (args.size()!=2) THROW(fatal_error,"Max requires 2 arguments.")
  double arg0=ToType<double>(args[0]);
  double arg1=ToType<double>(args[1]);
  return ToString(Max(arg0,arg1));
}

DEFINE_INTERPRETER_FUNCTION(Resolve_Bracket)
{
  if (expr.find("(")==std::string::npos ||
      expr.find(")")==std::string::npos) return expr;
  static int cnt=0;
  if (++cnt==100) exit(0);
  int open=0, take=1;
  size_t l=std::string::npos, r=std::string::npos;
  for (size_t i=0;i<expr.length();++i) {
    if (expr[i]=='(') {
      ++open;
      if (l==std::string::npos) {
	Algebra_Interpreter::Function_Set::const_iterator fit=
	  p_interpreter->Functions().begin();
	for (;fit!=p_interpreter->Functions().end();++fit) {
	  size_t pos=expr.rfind((*fit)->Tag(),i);
	  if (pos!=std::string::npos &&
	      pos+(*fit)->Tag().length()==i) break;
	}
	if (fit==p_interpreter->Functions().end()) {
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
  if (l==std::string::npos) return expr;
  std::string left=expr.substr(0,l);
  std::string right=expr.substr(r+1);
  msg_Tracking()<<"Resolve_Bracket -> '"
		<<left<<"' '"<<expr.substr(l+1,r-l-1)<<"' '"<<right<<"'\n";
  return p_interpreter->
    Iterate(left+p_interpreter->
	    Iterate(expr.substr(l+1,r-l-1))+right);
}

DEFINE_INTERPRETER_FUNCTION(Interprete_Function)
{
  if (expr.find("(")==std::string::npos ||
      expr.find(")")==std::string::npos) return expr;
  Function *func=NULL;
  size_t pos=std::string::npos;
  for (Algebra_Interpreter::Function_Set::const_reverse_iterator 
	 fit=p_interpreter->Functions().rbegin();
       fit!=p_interpreter->Functions().rend();++fit) 
    if ((pos=expr.rfind((*fit)->Tag()))!=std::string::npos) {
      func=*fit;
      break;
    }
  if (func==NULL) return expr;
  size_t last=pos+func->Tag().length()+1, open=0;
  std::vector<std::string> args;
  size_t i=last-1;
  for (;i<expr.length();++i) {
    if (expr[i]=='(') ++open;
    if (expr[i]==')') --open;
    if (open==0 || expr[i]==',') {
      args.push_back(expr.substr(last,i-last));
      last=i+1;
      if (open==0) break;
    }
  }
  std::string left=expr.substr(0,pos);
  std::string right=expr.substr(i+1);
  if (msg.LevelIsTracking()) {
    std::string out=args[0];
    for (size_t j=1;j<args.size();++j) out+=","+args[j];
    msg_Tracking()<<"Interprete_Function -> '"<<left
		  <<"' '"<<func->Tag()<<"("<<out
		  <<")' '"<<right<<"'\n";
  }
  for (size_t j=0;j<args.size();++j) 
    args[j]=p_interpreter->Iterate(args[j]);
  return p_interpreter->Iterate(left+func->Evaluate(args)+right);
}

size_t FindBinaryPlus(const std::string &expr,const bool fwd,
		      size_t cpos=std::string::npos)
{
  if (cpos==std::string::npos && fwd) cpos=0;
  size_t pos(fwd?expr.find("+",cpos):expr.rfind("+",cpos));
  if (pos==std::string::npos || (pos==0 && !fwd)) return std::string::npos;
  if (pos==0) return FindBinaryPlus(expr,fwd,1);
  if (expr[pos-1]=='e' || expr[pos-1]=='E' ||
      expr[pos-1]<48 || expr[pos-1]>57) 
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
  if (expr[pos-1]=='e' || expr[pos-1]=='E' ||
      expr[pos-1]<48 || expr[pos-1]>57) 
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
  size_t pos=std::string::npos, tpos=std::string::npos;
  for (Algebra_Interpreter::Operator_Map::const_reverse_iterator 
	 oit=p_interpreter->Operators().rbegin();
       oit!=p_interpreter->Operators().rend();++oit) {
    if (oit->second->Tag()=="!") pos=FindUnaryNot(expr,false);
    else if (oit->second->Tag()=="-") pos=FindBinaryMinus(expr,false);
    else if (oit->second->Tag()=="+") pos=FindBinaryPlus(expr,false);
    else pos=expr.find(oit->second->Tag());
    if (pos!=std::string::npos) {
      if (oit->second->Binary() && pos!=0) {
	size_t mpos=pos+oit->second->Tag().length();
	if (mpos<expr.length() && expr[mpos]=='-') tpos=pos;
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
  std::vector<std::string> args(2);
  args[0]=p_interpreter->Iterate(lstr);
  args[1]=p_interpreter->Iterate(rstr);
  msg_Tracking()<<"Interprete_Binary -> '"
	    <<lrstr<<"' '"<<lstr<<"' '"<<op->Tag()
	    <<"' '"<<rstr<<"' '"<<rrstr<<"'\n";
  return p_interpreter->
    Iterate(lrstr+op->Evaluate(args)+rrstr);
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
    if (pos!=std::string::npos)
      if (!oit->second->Binary()) {
	op=oit->second;
	break;
      }
      else {
	return expr;
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
  std::vector<std::string> args(1);
  args[0]=p_interpreter->Iterate(rstr);
  msg_Tracking()<<"Interprete_Unary -> '"
		<<lrstr<<"' '"<<op->Tag()<<"' '"<<rstr<<"' '"<<rrstr<<"'\n";
  return p_interpreter->
    Iterate(lrstr+op->Evaluate(args)+rrstr);
}

Algebra_Interpreter::Algebra_Interpreter(const bool standard):
  p_replacer(this)
{
  m_interpreters.insert(new Resolve_Bracket(this));
  m_interpreters.insert(new Interprete_Binary(this));
  m_interpreters.insert(new Interprete_Unary(this));
  m_interpreters.insert(new Interprete_Function(this));
  if (!standard) return;
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
  AddOperator(new Unary_Not());
  AddFunction(new Power());
  AddFunction(new Logarithm());
  AddFunction(new Logarithm10());
  AddFunction(new Exponential());
  AddFunction(new Absolute_Value());
  AddFunction(new Square());
  AddFunction(new Square_Root());
  AddFunction(new Minimum());
  AddFunction(new Maximum());
}

Algebra_Interpreter::~Algebra_Interpreter()
{
  while (!m_operators.empty()) {
    delete m_operators.begin()->second;
    m_operators.erase(m_operators.begin());
  }
  while (!m_functions.empty()) {
    delete *m_functions.begin();
    m_functions.erase(m_functions.begin());
  }
}

std::string Algebra_Interpreter::Interprete(const std::string &expr)
{
  msg_Tracking()<<"Algebra_Interpreter::Interprete("<<expr<<")\n";
  std::string res=expr;
  KillBlanks(res);
  res=p_replacer->ReplaceTags(res);
  return Iterate(res);
}

void Algebra_Interpreter::AddFunction(Function *const f)
{
  m_functions.insert(f); 
}

void Algebra_Interpreter::AddOperator(Operator *const b)
{ 
  m_operators.insert(Operator_Pair(b->Priority(),b)); 
}

std::string Algebra_Interpreter::Iterate(const std::string &expr)
{
  static size_t depth=0;
  if (++depth>1000) THROW(critical_error,"Max depth reached.");
  msg_Indent();
  std::string res=expr;
  Interpreter_Set::const_iterator iit=m_interpreters.begin();
  for (;iit!=m_interpreters.end();++iit) 
    res=(*iit)->Interprete(res);
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
  msg_Tracking()<<"Algebra_Interpreter::ReplaceTags("<<expr<<")\n";
  size_t pos=std::string::npos;
  for (String_Map::const_iterator sit=m_tags.begin();
       sit!=m_tags.end();++sit) {
    if ((pos=expr.find(sit->first))!=std::string::npos) 
      return ReplaceTags(expr.replace(pos,sit->first.length(),sit->second));}
  return expr;
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
