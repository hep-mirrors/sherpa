#ifdef COMPILE__Getter_Function
#ifndef OBJECT_TYPE
#error object type undefined 
#error specify an object type using #define OBJECT_TYPE Object_Type 
#endif
#ifndef PARAMETER_TYPE
#error parameter type undefined
#error specify a parameter type using #define PARAMETER_TYPE Parameter_Type 
#endif

#include "Getter_Function.H"
#include "Exception.H"
#include "Message.H"
#include <iomanip>

using namespace ATOOLS;

template<class ObjectType,class ParameterType>
typename Getter_Function<ObjectType,ParameterType>::String_Getter_Map *
Getter_Function<ObjectType,ParameterType>::s_getters=NULL;

template<class ObjectType,class ParameterType>
Getter_Function<ObjectType,ParameterType>::
Getter_Function(const std::string &name)
{
  static bool initialized=false;
  if (!initialized) {
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
    std::cout<<"Getter_Function::Getter_Function(\""<<name<<"\"): "
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
}

template<class ObjectType,class ParameterType>
void Getter_Function<ObjectType,ParameterType>::
PrintInfo(std::ostream &str) const
{
  str<<"No Information";
}

template<class ObjectType,class ParameterType>
ObjectType *const Getter_Function<ObjectType,ParameterType>::
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
  const std::ios_base::fmtflags def=str.flags();
  str.setf(std::ios_base::left,std::ios_base::adjustfield);
  for (typename String_Getter_Map::const_iterator git=s_getters->begin();
       git!=s_getters->end();++git) {
    str<<"   "<<std::setw(width)<<git->first<<" ";
    git->second->PrintInfo(str);
    str<<"\n";
  }
  str.setf(def);
}

template<class ObjectType,class ParameterType>
ObjectType *const Getter_Function<ObjectType,ParameterType>::
GetObject(const std::string &name,const Parameter_Type &parameters)
{
  if (s_getters->find(name)!=s_getters->end()) return (*(*s_getters)[name])(parameters);
  return NULL;
}

template Getter_Function<OBJECT_TYPE,PARAMETER_TYPE>;

#endif
