#ifndef Data_Read_H
#define Data_Read_H

#include <string>
#include <map>
#include <fstream>

#include "Flavour.H"
#include "Algebra_Interpreter.H"

namespace ATOOLS {  
  //! The On/Off switch.
  class Switch {
  public:
    enum code { Off     = 0, 
		On      = 1,
		Unknown = 99}; 
  };

  //! different Beam realization  
  class Beam_Type {
  public:
    enum code { Monochromatic  = 0,
		Gaussian       = 1,
	        Laser_Back     = 2,
	        Simple_Compton = 3,
		Spec_Read      = 4,
	        Unknown        = 99
    };
  };

  class Beam_Shape {
  public:
    enum code { Cylinder          = 0,
                Gaussian_Cylinder = 1,
		Unknown           = 99
    };
  };

  class Beam_Generator {
  public:
    enum code { Internal          = 0,
		Unknown           = 99
    };
  };

  //! The different String types
  class String_Type {
  public:
    enum code { NoString = 0, 
	        String   = 1, 
	        Library  = 2, 
	        Unknown  = 99
    };
  };

  //! The different Model types
  class Model_Type {
  public:
    enum code { pure_QCD      = 0, 
	        QCD           = 1, 
	        pure_EW       = 2, 
	        SM            = 3, 
	        MSSM          = 4, 
	        THDM          = 5,
	        MSSM_EHC      = 4, 
	        ADD           = 7,
	        SMEHC         = 8,
		SM_ZPrime     = 9,
		SM_AGC        = 10,
	        SM_Phantom_U1 = 11, 
	        Unknown       = 99
    };
  };

  //! different (QED) - ISR realization  
  class ISR_Type {
  public:
    enum code { No             = 0, 
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

  typedef std::map<std::string,std::string> Parameter_Map;
  typedef Parameter_Map::iterator           Parameter_Iterator;

  class Data_Read: public Tag_Replacer {
    //associative map (parameter name and coresponding value)
    Parameter_Map m_parameters;
    static Parameter_Map s_commandlineparameters, s_tags;

    bool m_fileexists;
    
    std::string m_filename;

    void FillIn(std::string buffer) ;
    static void Shorten(std::string& str);
    // returns sum of digits
    int Crossfoot(std::string name);
    void AddCommandLine();
    template <class Type>
    const Type ReturnData(const std::string &name,const Type type);
    std::string ReplaceTags(std::string &expr) const;    
  public:
    Data_Read(std::string filename, bool ignoremissingfile=false);
    void ReadIn(std::string filename, bool ignoremissingfile=false);
    void WriteOut(std::string filename,int flag=std::ios::out);

    template <class Type> 
      Type GetValue(std::string name);

    template <class Type> 
      Type GetValue(std::string name,Type default_value);

    //    template <class Type> 
    void SetValue(std::string name, std::string value);

    static void SetCommandLine(std::string name, std::string value);

    std::string GenerateKey();
    
    bool FileExists() { return m_fileexists; }
    inline std::string FileName() const { return m_filename; }
    inline static void SetTags(const Parameter_Map &tags) { s_tags=tags; }
    inline static void ResetTags() { s_tags.clear(); }
  };  

  // specialisation

  template <>  std::string        Data_Read::GetValue<std::string>(std::string);
  template <>  Switch::code       Data_Read::GetValue<Switch::code>(std::string);
  template <>  Model_Type::code   Data_Read::GetValue<Model_Type::code>(std::string);
  template <>  Flavour Data_Read::GetValue<Flavour>(std::string);
  template <>  Beam_Type::code    Data_Read::GetValue<Beam_Type::code>(std::string);
  template <>  Beam_Generator::code    Data_Read::GetValue<Beam_Generator::code>(std::string);
  template <>  Beam_Shape::code   Data_Read::GetValue<Beam_Shape::code>(std::string);
  template <>  ISR_Type::code     Data_Read::GetValue<ISR_Type::code>(std::string);
  template <>  String_Type::code  Data_Read::GetValue<String_Type::code>(std::string);


  // do we need those:  
  /*
  template <>  Switch::code       Data_Read::GetValue<Switch::code>(std::string,Switch::code);
  template <>  Model_Type::code   Data_Read::GetValue<Model_Type::code>(std::string,Model_Type::code);
  template <>  Flavour Data_Read::GetValue<Flavour>(std::string,Flavour);
  template <>  Beam_Type::code    Data_Read::GetValue<Beam_Type::code>(std::string,Beam_Type::code);
  template <>  ISR_Type::code     Data_Read::GetValue<ISR_Type::code>(std::string,ISR_Type::code);
  template <>  String_Type::code  Data_Read::GetValue<String_Type::code>(std::string,String_Type::code);
  */

  /*!
    \file contains the class Data_Read
  */

  /*!
    \class Data_Read
    \brief  is an all purpose class for reading in parameters.

    This is the data file reading class. It can be used for all different kind of parameters.
    It is capable of reading standard types (e.g. int, double, string,...) as well
    as a few user defined types (e.g. Switch::code, Flavour,...) from a data file like:
    \verbatim
      ! comment line
      BEAM1      =  P+
      BEAM2      =  P-
      CMSENERGY  = 2000.0  ! note: in GeV
    \endverbatim  
    cf. Run_Parameter.C for example usage                                  
  */

  /*!
    \fn    Data_Read::Data_Read(std::string filename)
    \brief Constructor, reads a datafile via ReadIn()
  */

  /*!
    \fn    void Data_Read::ReadIn(std::string filename)    
    \brief reads a datafile in buffer
    
    This routine reads a parameter file in buffer. It can be called more than once to read 
    a number of data files, e.g. "default.dat" and "user.dat".
    The last appearance of a key word actually defines the value of it
  */

  /*!
    \fn    void Data_Read::WriteOut(std::string filename,int flag=std::ios::out)
    \brief write out or append to file (ios_base::out,ios_base::app)
  */

  /*!
    \fn    Type Data_Read::GetValue(std::string name);
    \brief returns value for given name or zero in case it did not appear in the read data file
    
    \internal
    (template function needs to be specialised for user defined types/classes) 
       or the in and output operators of the iostream might be overloaded      
  */

  /*!
    \fn    Type Data_Read::GetValue(std::string name,Type default_value)
    \brief returns value for given parameter name or a given default value in case
        it did not appear in the data file
  */

  /*!
    \fn    void Data_Read::SetValue(std::string name, Type value)
    \brief set a value for a parameter
  */

  /*!
    \fn    std::string Data_Read::GenerateKey()
    \brief generates a (hopefully) unique identifier
  */

  /*!   
    \fn    int Crossfoot(std::string name)
    \brief returns sum of digits

  */
} // namespace ATOOLS

#endif  // Data_Read_H
