#ifndef Data_Read_H
#define Data_Read_H

#include <string>
#include <map>
#include <fstream>

#include "Flavour.H"

namespace AORGTOOLS {  
  //! The On/Off switch.
  class Switch {
  public:
    enum code { Off     = 0, 
		On      = 1,
		Unknown = 99}; 
  };

  
  //! The different String types
  class String_Type {
  public:
    enum code {NoString = 0, 
	       String   = 1, 
	       Library  = 2, 
	       Unknown  = 99
    };
  };

  //! The different Model types
  class Model_Type {
  public:
    enum code {pure_QCD = 0, 
	       QCD      = 1, 
	       EW       = 2, 
	       SM       = 3, 
	       MSSM     = 4, 
	       THDM     = 5,
               NCQED    = 6,
	       Unknown  = 99
    };
  };

  // different (QED) - ISR realization  
  class ISR_Type {
  public:
    enum code {No             = 0, 
	       Simple_Struc   = 1,
	       Extended_Struc = 2,
	       PA             = 3,
	       EA             = 4,
	       MPA            = 5,
	       KoralZ         = 6,
	       Pythia         = 7,
	       KKMC           = 8,
	       Unknown        = 99
    };
  };


  // different Beam realization  
  class Beam_Type {
  public:
    enum code {No             = 0, 
	       Laser_Back     = 1,
	       Unknown        = 99
    };
  };


  //! The Data file reading class; it can be used for all different kind of parameters.
  /*! 
    It is capable of reading standard types (e.g. int, double, string,...) as well
    as a few user defined types (e.g. Switch::code, Flavour,...) from a data file like:
      ! comment line
      BEAM1      =  P+
      BEAM2      =  P-
      CMSENERGY  = 2000.0  ! note: in GeV
    cf. Run_Parameter.C for example usage                                  
  */
  typedef std::map<std::string,std::string> Parameter_Map;
  typedef Parameter_Map::iterator           Parameter_Iterator;

  class Data_Read {
    //associative map (parameter name and coresponding value)
    Parameter_Map parameters;
    
    void FillIn(char * dummy) ;
    void Shorten(std::string& str);
  public:
    // constructor ; reads a datafile via ReadIn()
    Data_Read(std::string filename);
    
    // reads a datafile in buffer
    /* can be used to read more than one data file, e.g. "default.dat" and "user.dat"
       the last apearance of a key word actually defines the value of it */
    void ReadIn(std::string filename);

    // write out or append to file (ios_base::out,ios_base::app)
    void WriteOut(std::string filename,int flag=std::ios::out);

    // returns value for given name or zero in case it did not appear in the read data file
    /* (template function needs to be specialised for user defined types/classes) 
       or the in and output operators of the iostream might be overloaded         */   
    template <class Type> 
      Type GetValue(std::string name);

    template <class Type> 
      Type GetValue(std::string name,Type default_value);

    template <class Type> 
      void SetValue(std::string name, Type value);

    // generates a (hopefully) unique identifier
    std::string GenerateKey();
    // returns sum of digits
    int Crossfoot(std::string name);

  };  

  // specialisation

  template <>  Switch::code       Data_Read::GetValue<Switch::code>(std::string);
  template <>  Model_Type::code   Data_Read::GetValue<Model_Type::code>(std::string);
  template <>  APHYTOOLS::Flavour Data_Read::GetValue<APHYTOOLS::Flavour>(std::string);
  template <>  Beam_Type::code    Data_Read::GetValue<Beam_Type::code>(std::string);
  template <>  ISR_Type::code     Data_Read::GetValue<ISR_Type::code>(std::string);
  template <>  String_Type::code  Data_Read::GetValue<String_Type::code>(std::string);

  /*
  template <>  void  Data_Read::SetValue<Switch::code>(std::string,Switch::code);
  template <>  void  Data_Read::SetValue<Model_Type::code>(std::string,Model_Type::code);
  template <>  void  Data_Read::SetValue<APHYTOOLS::Flavour>(std::string,APHYTOOLS::Flavour);
  template <>  void  Data_Read::SetValue<ISR_Type::code>(std::string,ISR_Type::code);
  template <>  void  Data_Read::SetValue<String_Type::code>(std::string,String_Type::code);
  */


} // namespace AORGTOOLS

#endif  // Data_Read_H
