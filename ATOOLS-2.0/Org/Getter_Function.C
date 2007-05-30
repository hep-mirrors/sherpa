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

#include "Getter_Function.H"
#include "STL_Tools.H"
#include "Exception.H"
#include "Message.H"
#include "MyStrStream.H"
#include <typeinfo>

using namespace ATOOLS;

template<class ObjectType,class ParameterType>
typename Getter_Function<ObjectType,ParameterType>::String_Getter_Map *
Getter_Function<ObjectType,ParameterType>::s_getters=NULL;

template<class ObjectType,class ParameterType>
bool Getter_Function<ObjectType,ParameterType>::s_exactmatch=EXACTMATCH;

template<class ObjectType,class ParameterType>
Getter_Function<ObjectType,ParameterType>::
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
	   <<"Added getter '"<<this<<"' -> \""<<name<<"\"."<<std::endl;
#endif
  if (s_getters->find(name)==s_getters->end()) {
    s_getters->insert(std::pair<const std::string,
		      Getter_Function *const>(name,this));
  }
  else {
    std::cout<<"Getter_Function<"<<typeid(ObjectType*).name()
	     <<","<<typeid(ParameterType*).name()<<">::"
	     <<"Getter_Function(\""<<name<<"\"): "
	     <<"Doubled identifier. Abort."<<std::endl;
    abort();
  }
}

template<class ObjectType,class ParameterType>
Getter_Function<ObjectType,ParameterType>::~Getter_Function()
{
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

template<class ObjectType,class ParameterType>
void Getter_Function<ObjectType,ParameterType>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"No Information";
}

template<class ObjectType,class ParameterType>
ObjectType * Getter_Function<ObjectType,ParameterType>::
operator()(const Parameter_Type &parameters) const
{
  std::cout<<"Getter_Function::operator(): "
	   <<"Virtual function called."<<std::endl;
  return NULL;
}

template<class ObjectType,class ParameterType>
void Getter_Function<ObjectType,ParameterType>::
PrintGetterInfo(std::ostream &str,const size_t width)
{
  const IOS_BASE::fmtflags def=str.flags();
  str.setf(IOS_BASE::left,IOS_BASE::adjustfield);
  for (typename String_Getter_Map::const_iterator git=s_getters->begin();
       git!=s_getters->end();++git) {
    if (!git->second->m_display) continue;
    str<<"   "<<std::setw(width)<<git->first<<" ";
    git->second->PrintInfo(str,width);
    str<<"\n";
  }
  str.setf(def);
}

template<class ObjectType,class ParameterType>
ObjectType *Getter_Function<ObjectType,ParameterType>::
GetObject(const std::string &name,const Parameter_Type &parameters)
{
  if (!s_exactmatch) {
    for (typename String_Getter_Map::reverse_iterator git=s_getters->rbegin();
	 git!=s_getters->rend();++git) {
      if ((name.length()==0 && git->first.length()==0) ||
	  (git->first.length()>0 && name.find(git->first)!=std::string::npos))
	return (*git->second)(parameters);
    }
    return NULL;
  }
  typename String_Getter_Map::iterator git=s_getters->find(name);
  if (git!=s_getters->end()) return (*git->second)(parameters);
  return NULL;
}

template class Getter_Function<OBJECT_TYPE,PARAMETER_TYPE>;

#endif
