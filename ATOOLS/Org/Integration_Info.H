#ifndef ATOOLS_Org_Integration_Info_H
#define ATOOLS_Org_Integration_Info_H

#ifndef ATOOLS_Org_Info_Key_H
#error The header 'Integration_Info.H' must not be included directly. \
Please include the header 'Info_Key.H' instead.
#else

#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <map>
#include <vector>

namespace ATOOLS {

  const double UNDEFINED_LOWER=-std::numeric_limits<double>::max();
  const double UNDEFINED_UPPER=std::numeric_limits<double>::max();
  const double UNDEFINED_DOUBLE=std::numeric_limits<double>::max();
  const double UNDEFINED_WEIGHT=0.;

  const ATOOLS::Vec4D UNDEFINED_VECTOR=ATOOLS::Vec4D();

  class Integration_Info {
  public:

    typedef std::vector<Info_Key*>                     Key_Vector;
    typedef std::pair<size_t,Key_Vector>               SizeT_KeyVector_Pair;
    typedef std::map<std::string,SizeT_KeyVector_Pair> String_KeyPair_Map;
    typedef std::pair<size_t,String_KeyPair_Map>       SizeT_KeyMap_Pair;
    typedef std::map<std::string,SizeT_KeyMap_Pair>    String_MapPair_Map;

  private:

    String_MapPair_Map m_keymap;

    std::vector<Double_Container> m_doubles;
    std::vector<Vector_Container> m_vectors;
    std::vector<Double_Container> m_weights;
    std::vector<si::code>         m_status;

    void AssignKey(Info_Key &key,const size_t doubles,const size_t vectors);
    void ReleaseKey(Info_Key &key);

    Double_Container &Doubles(const size_t valuekey);
    Vector_Container &Vectors(const size_t valuekey);

    double        &Double(const size_t valuekey,const size_t i);
    ATOOLS::Vec4D &Vector(const size_t valuekey,const size_t i);

    double        Double(const size_t valuekey,const size_t i) const;
    ATOOLS::Vec4D Vector(const size_t valuekey,const size_t i) const;

    void   SetWeight(const size_t valuekey,const size_t weightkey,
		     const double weight);
    double Weight(const size_t valuekey,const size_t weightkey) const;
    
    si::code SetStatus(const size_t statuskey,const si::code status);
    si::code Status(const size_t statuskey);

    friend class Info_Key;

    friend std::ostream &operator<<(std::ostream &str,
				    const Integration_Info &info);
    friend std::ostream &operator<<(std::ostream &str,const Info_Key &key);

  public:

    // constructor
    Integration_Info();

    // destructor
    ~Integration_Info();

    // member functions
    void ResetAll();

  };// end of class Integration_Info

  std::ostream &operator<<(std::ostream &str,const Double_Container &doubles);
  std::ostream &operator<<(std::ostream &str,const Vector_Container &vectors);

  /*!
    \class Integration_Info
    \brief Provides fast exchange of information between channels

    This class provides fast and easy information exchange between the channels 
    of the phase space integrator.
    Information is stored in a std::vector<Double_Container> and a 
    std::vector<Vector_Container> respectively to provide constant time access.
    Weight storage and status information is provided too, 
    such that multiple calculation of weights becomes obsolete.

    Each Channel evaluating a certain variable may access Integration_Info 
    through a corresponding Info_Key which communicates with Integration_Info 
    automatically.
    The appropriate way is to create a key like in the example below
    \code
       class Example_Channel: public Single_Channel {
       private:
         Info_Key m_key;
       public:
         ...
       };

       Example_Channel::Example_Channel(Integration_Info *info)
       {
         // <variable_info> is the exponent for s', for example
         m_key.SetInfo("<variable_info>")
	 // <variable_name> is "s'", for example
         m_key.Assign("<variable_name>",<value_size>,<vector_size>,info);
       }
    \endcode
    The usage of keys is described in the documentation of Info_Key.

    Direct access to Integration_Info is prohibited to preserve 
    the corruption of data.
  */

}// end of namespace ATOOLS

#include "ATOOLS/Org/Integration_Info.inl.H"

#endif

#endif
