#ifdef COMPILE__Getter_Function
#ifndef OBJECT_TYPE
#error object type undefined 
#error specify an object type using #define OBJECT_TYPE Object_Type 
#endif
#ifndef PARAMETER_TYPE
#error parameter type undefined
#error specify a parameter type using #define PARAMETER_TYPE Parameter_Type 
#endif
#ifndef EXACTMATCH
#define EXACTMATCH true
#endif
#ifndef SORT_CRITERION
#define SORT_CRITERION std::less<std::string>
#endif

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <typeinfo>
#include <cstdlib>

using namespace ATOOLS;

template<class ObjectType,class ParameterType,class SortCriterion>
typename Getter_Function<ObjectType,ParameterType,SortCriterion>::
String_Getter_Map *Getter_Function<ObjectType,ParameterType,SortCriterion>::
s_getters=NULL;

template<class ObjectType,class ParameterType,class SortCriterion>
bool Getter_Function<ObjectType,ParameterType,SortCriterion>::
s_exactmatch=EXACTMATCH;

template<class ObjectType,class ParameterType,class SortCriterion>
Getter_Function<ObjectType,ParameterType,SortCriterion>::
Getter_Function(const std::string &name):
  m_display(true)
{
  static bool initialized=false;
  if (!initialized || s_getters==NULL) {
    s_getters = new String_Getter_Map();
    initialized=true;
  }
#ifdef DEBUG__Getter_Function
  std::cout<<"Getter_Function::Getter_Function(..): "
	   <<"Added getter '"<<this<<"'("<<Demangle(typeid(*this).name())
	   <<") -> \""<<name<<"\"."<<std::endl;
#endif
  typename String_Getter_Map::iterator git(s_getters->find(name));
  if (git!=s_getters->end()) {
    std::cout<<std::string(80,'#')<<std::endl;
    std::cout<<"Getter_Function<"<<Demangle(typeid(ObjectType*).name())
	     <<","<<Demangle(typeid(ParameterType*).name())<<"> {\n"
	     <<"  Doubled identifier \""<<name<<"\"!\n  Now replacing '"
	     <<Demangle(typeid(*git->second).name())<<"'.\n  "
	     <<"This operation may lead to wrong results "
	     <<"or a program crash.\n}"<<std::endl;
    std::cout<<std::string(80,'#')<<std::endl;
    s_getters->erase(git);
  }
  s_getters->insert(std::pair<const std::string,
		    Getter_Function *const>(name,this));
}

template<class ObjectType,class ParameterType,class SortCriterion>
Getter_Function<ObjectType,ParameterType,SortCriterion>::~Getter_Function()
{
  if (s_getters==NULL) return;
  for (typename String_Getter_Map::iterator git=s_getters->begin();
       git!=s_getters->end();++git) {
    if (git->second==this) {
      s_getters->erase(git);
      break;
    }
  }
  if (s_getters->empty()) {
    delete s_getters;
    s_getters=NULL;
  }
}

template<class ObjectType,class ParameterType,class SortCriterion>
void Getter_Function<ObjectType,ParameterType,SortCriterion>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<Demangle(typeid(*this).name());
}

template<class ObjectType,class ParameterType,class SortCriterion>
ObjectType * Getter_Function<ObjectType,ParameterType,SortCriterion>::
operator()(const Parameter_Type &parameters) const
{
  std::cout<<"Getter_Function::operator(): "
	   <<"Virtual function called."<<std::endl;
  return NULL;
}

template<class ObjectType,class ParameterType,class SortCriterion>
void Getter_Function<ObjectType,ParameterType,SortCriterion>::
PrintGetterInfo(std::ostream &str,const size_t width, const std::string &indent,
                const std::string &separator, const std::string &lineend,
                const std::string &replacefrom, const std::string &replaceto)
{
  if (s_getters==NULL) return;
  const IOS_BASE::fmtflags def=str.flags();
  str.setf(IOS_BASE::left,IOS_BASE::adjustfield);
  for (typename String_Getter_Map::const_iterator git=s_getters->begin();
       git!=s_getters->end();++git)
  {
    if (!git->second->m_display) continue;
    std::string escapedname=StringReplace(git->first,replacefrom,replaceto);
    str<<indent<<std::setw(width)<<escapedname<<separator;
    git->second->PrintInfo(str,width);
    str<<lineend;
  }
  str.setf(def);
}

template<class ObjectType,class ParameterType,class SortCriterion>
ObjectType *Getter_Function<ObjectType,ParameterType,SortCriterion>::
GetObject(const Parameter_Type &parameters) const
{
  return (*this)(parameters);
}

template<class ObjectType,class ParameterType,class SortCriterion>
ObjectType *Getter_Function<ObjectType,ParameterType,SortCriterion>::
GetObject(const std::string &name,const Parameter_Type &parameters)
{
  if (s_getters==NULL) return NULL;
  if (!s_exactmatch) {
    for (typename String_Getter_Map::reverse_iterator git=s_getters->rbegin();
	 git!=s_getters->rend();++git) {
      if ((name.length()==0 && git->first.length()==0) ||
	  (git->first.length()>0 && name.find(git->first)==0))
	return (*git->second)(parameters);
    }
    return NULL;
  }
  typename String_Getter_Map::iterator git=s_getters->find(name);
  if (git!=s_getters->end()) return (*git->second)(parameters);
  return NULL;
}

template<class ObjectType,class ParameterType,class SortCriterion>
std::vector<const Getter_Function<ObjectType,ParameterType,SortCriterion> *>
Getter_Function<ObjectType,ParameterType,SortCriterion>::
GetGetters(const std::string &name)
{
  Getter_List list;
  if (s_getters==NULL) return list;
  for (typename String_Getter_Map::reverse_iterator git=s_getters->rbegin();
       git!=s_getters->rend();++git) {
    if (name.length()==0 ||
	(git->first.length()>0 && git->first.find(name)!=std::string::npos))
    list.push_back(git->second);
  }
  return list;
}

template class Getter_Function<OBJECT_TYPE,PARAMETER_TYPE,SORT_CRITERION>;

#endif
