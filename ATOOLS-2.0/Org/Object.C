#include "Object.H"

#include "Scaling.H"

using namespace ATOOLS;

Object::String_Object_Map Object::s_objects;

Object::Object()
{
  s_objects.insert(std::pair<const std::string,
		   Object *const>(ToString(this),this));
}

Object::Object(const Object &reference)
{
  String_Object_Map::iterator oit=s_objects.begin();
  for (;oit!=s_objects.end();++oit) {
    if (oit->second==&reference) break;
  }
  size_t i=0;
  for (;i<std::string::npos;++i) {
    if (s_objects.find(oit->first+" "+ToString(i))!=s_objects.end()) break;
  }
  s_objects.insert(std::pair<const std::string,
		   Object *const>(oit->first+" "+ToString(i),this));
}

Object::Object(const std::string name)
{
  if (s_objects.find(name)==s_objects.end()) {
    s_objects.insert(std::pair<const std::string,
		     Object *const>(name,this));
    return;
  }
  size_t i=0;
  for (;i<std::string::npos;++i) {
    if (s_objects.find(name+" "+ToString(i))!=s_objects.end()) break;
  }
  s_objects.insert(std::pair<const std::string,
		   Object *const>(name+" "+ToString(i),this));
}

Object::~Object()
{
  for (String_Object_Map::iterator oit=s_objects.begin();
       oit!=s_objects.end();++oit) {
    if (oit->second==this) {
      s_objects.erase(oit);
      break;
    }
  }  
}

void Object::PrintObjects(std::ofstream &str)
{
  str<<"Object::PrintObjects(..): {\n";
  for (String_Object_Map::iterator oit=s_objects.begin();
       oit!=s_objects.end();++oit) {
    str<<"   '"<<oit->first<<"' -> ("<<oit->second<<")\n";
  }  
  str<<"}"<<std::endl;
}
