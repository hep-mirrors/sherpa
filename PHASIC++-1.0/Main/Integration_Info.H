#ifndef Integration_Info_H
#define Integration_Info_H

#include "Vector.H"
#include <map>
#include <vector>
#include <limits>
#include <ostream>

namespace PHASIC {

  const double UNDEFINED_LOWER=-std::numeric_limits<double>::max();
  const double UNDEFINED_UPPER=std::numeric_limits<double>::max();
  const double UNDEFINED_DOUBLE=std::numeric_limits<double>::max();
  const double UNDEFINED_WEIGHT=0.;

  const ATOOLS::Vec4D UNDEFINED_VECTOR=ATOOLS::Vec4D();

  typedef std::vector<double>        Double_Container;
  typedef std::vector<ATOOLS::Vec4D> Vector_Container;

  namespace si {

    enum code { idle      =  0,
		reset     =  1,
		diced     =  2,
		weighted  =  3,
		error     =  4  };

  }// end of namespace si

  class Integration_Info;

  class Info_Key {
  private:

    Integration_Info *p_info;
    std::string       m_name, m_info;
    
    size_t m_valuekey, m_weightkey;

    friend class Integration_Info;

  public:

    // constructors
    Info_Key();

    // destructor
    ~Info_Key();

    // setters
    void     SetInfo(const std::string info);
    si::code SetStatus(const si::code sicode);

    // getters
    const std::string &Name() const;
    const std::string &Info() const;

    si::code Status() const;

    // member functions
    void Assign(const std::string name,const size_t doubles,
		const size_t vectors,Integration_Info *info);

    Double_Container &Doubles();
    Vector_Container &Vectors();

    double        Double(const size_t i);
    ATOOLS::Vec4D Vector(const size_t i);
    double        Weight() const;

    double        &operator[](const size_t i);
    ATOOLS::Vec4D &operator()(const size_t i);

    void operator<<(const double weight);

    friend std::ostream &operator<<(std::ostream &str,const Info_Key &key);

  };// end of class Info_Key

  class Integration_Info {
  private:

    std::map<std::string,std::pair<size_t,std::vector<Info_Key*> > > m_keymap;

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

    void   SetWeight(const size_t valuekey,const size_t weightkey,const double weight);
    double Weight(const size_t valuekey,const size_t weightkey) const;
    
    si::code SetStatus(const size_t statuskey,const si::code status);
    si::code Status(const size_t statuskey);

    friend class Info_Key;

    friend std::ostream &operator<<(std::ostream &str,const Integration_Info &info);
    friend std::ostream &operator<<(std::ostream &str,const Info_Key &key);

  public:

    // constructor
    Integration_Info();

    // member functions
    void ResetAll();

  };// end of class Integration_Info

  std::ostream &operator<<(std::ostream &str,const Double_Container &doubles);
  std::ostream &operator<<(std::ostream &str,const Vector_Container &vectors);

  /*!
    \class Info_Key
    \brief Handles access to Integration_Info

    This class provides access to Integration_Info which handles value and 
    weight storage. Each Info_Key is assigned a reference point in Integration_Info 
    during Assign(<name>,<double_size>,<vector_size>,<info>).
    In the following the key calls Integration_Info to obtain and to set values
    and weights and to test for the status of the corresponding variable (see 
    the example below).
    \code
       class Example_ISR_Channel: public Single_Channel {
       private:
         double   m_exponent
         Info_Key m_spkey, ...;
       public:
         ...
       };

       void Example_ISR_Channel::GeneratePoint(const double *rans)
       {
         // m_spkey owns three values where m_spkey[0] -> lower limit
	 //                                 m_spkey[1] -> upper limit
	 //                                 m_spkey[2] -> value
         m_spkey[2]=CE::MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],ran[0]);
	 ...
       }

       void Example_ISR_Channel::GenerateWeight()
       {
         // operator<<(..) sets the unique weight
         m_spkey<<1./CE::MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[2]);
	 ...
       }
    \endcode
  */

  /*!
    \class Integration_Info
    \brief Provides fast exchange of information between channels

    This class provides fast and easy information exchange between the channels 
    of the phase space integrator.
    Information is stored in a std::vector<Double_Container> and a 
    std::vector<Vector_Container> respectively to provide constant time access.
    Weight storage and status information is provided too, such that multiple calculation 
    of weights becomes obsolete.

    Each Channel evaluating a certain variable may access Integration_Info through
    a corresponding Info_Key which communicates with Integration_Info automatically.
    The appropriate way is to create a key like in the example below
    \code
       class Example_Channel: public Single_Channel {
       private:
         Info_Key m_key;
       public:
         ...
       };

       Example_Channel::Example_Cahnnel(Integration_Info *info)
       {
         // <variable_info> is the exponent for s', for example
         m_key.SetInfo("<variable_info>")
	 // <variable_name> is "s'", for example
         m_key.Assign("<variable_name>",<value_size>,<vector_size>,info);
       }
    \endcode
    The usage of keys is described in the documentation of Info_Key.

    Direct access to Integration_Info is prohibited to preserve the corruption of 
    data.
  */

}// end of namespace PHASIC

#include "Integration_Info.IC"

#endif
