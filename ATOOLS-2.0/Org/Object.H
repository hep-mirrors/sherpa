#ifndef Object_H
#define Object_H

#include <map>

namespace ATOOLS {

  class Object {
  public:
    
    typedef std::map<const std::string,Object *const> String_Object_Map;

  private:
    
    static String_Object_Map s_objects;

  public:

    // constructors
    Object(const std::string name);
    Object(const Object &reference);
    Object();
    
    // destructor
    virtual ~Object();

    // member functions
    template <class Type> static Type *const GetObject(const std::string name)
    { 
      String_Object_Map::const_iterator pos=s_objects.find(name);
      if (pos!=s_objects.end())	return static_cast<Type *const>(pos->second); 
      return NULL;
    }

    template <class Type> static std::string &Name(const Type *object)
    { 
      for (String_Object_Map::const_iterator oit=s_objects.begin();
	   oit!=s_objects.end();++oit) if (oit->second==object) return oit->first;
      return "";
    }

  };//end of class Object

#define GET_OBJECT(TYPE,NAME) ATOOLS::Object::GetObject<TYPE>(NAME)

  /*!
    \class Object
    \brief Provides easy access to all kinds of objects
    
    The aim of this class is, to provide access to objects, 
    whose definitions do not reside in the module, where the access is needed.
    For example if some Remnant_Base, which belongs to the SHERPA
    module is needed by Simple_Chain in AMISIC, a pointer to it can easily 
    be obtained through the following procedure

    \code
    namespace SHERPA {

       class Remnant_Base: public ATOOLS::Object {
       private:

          // constructor
          Remnant_Base(const int beam);

          ...

       };// end of class Remnant_Base

    }// end of namespace SHERPA
    \endcode

    \code
    using namespace ATOOLS;

    Remnant_Base::Remnant_Base(const int beam):
       // assign each remnant a unique name
       Object(std::string("Remnant_Base_")+ToString(beam)))
    {
    
       ...
    
    }
    \endcode

    \code
    using namespace SHERPA;

    Simple_Chain::Initialize()
    {
       Remnant_Base remnants[2];
       // get pointer to first remnant using its name
       remnants[0]=GET_OBJECT(Remnant_Base,"Remnant_Base_1")
       // get pointer to second remnant using its name
       remnants[1]=GET_OBJECT(Remnant_Base,"Remnant_Base_2")

       ...

    }
    \endcode
  */

}// end of namespace ATOOLS

#endif
