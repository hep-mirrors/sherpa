#ifndef Color_Function_h
#define Color_Function_h

namespace AMEGIC {

  class cf {
  public:
    enum code {
      T=0, 
      F=1, 
      D=2, 
      None=3, 
      G=4, 
      Unknown=99}; 
  };
  

  class Color_Function {
  public:
    Color_Function()  {
      type = cf::Unknown;
      Next = 0;
      for (short int i=0;i<3;++i) partarg[i] = -1;
      for (short int i=0;i<3;++i) strarg[i] = '?';
    }
    
    Color_Function(cf::code _type, 
		   int _partarg0 = -1, int _partarg1 = -1, int _partarg2 = -1, 
		   char _strarg0 = '?', char _strarg1 = '?', char _strarg2 = '?') 
      :type(_type) { 
      
      Next = 0;
      partarg[0] = _partarg0;
      partarg[1] = _partarg1;
      partarg[2] = _partarg2;
      strarg[0]  = _strarg0;
      strarg[1]  = _strarg1;
      strarg[2]  = _strarg2;
      }
    
    cf::code  type; // Number of args, multi gluons etc.
    int  partarg[3];
    char strarg[3];
    Color_Function* Next;

    void SetParticleArg(int a, int b,int c=-1) {
      partarg[0] = a;partarg[1] = b;
      if (type!=cf::D && type!=cf::G) partarg[2] = c;
    }

    void SetStringArg(char a, char b, char c='?') {
      strarg[0] = a;strarg[1] = b;
      if (type!=cf::D && type!=cf::G) strarg[2] = c;
    }

    std::string String() {
      
      if (type==cf::None) return std::string("1");
      
      std::string help;
      switch (type) {
      case 0 :  help = std::string("T[ , , ]");break;
      case 1 :  help = std::string("F[ , , ]");break;
      case 2 :  help = std::string("D[ , ]");break;
      case 4 :  help = std::string("G[ , ]");break;
      default : return std::string("1");
      }
      for (short int i=0;i<3;i++) {
	if ((type==cf::D || type==cf::G) && i==2) break;
	help[2+i*2] = strarg[i];
      }
      return help;
    }
   
  };  

  // Format of a Color_Function to transport via MPI:
  struct MPI_Color_Function {
    
    int    m_type;
    int    m_partarg[3];         
    char   m_strarg[4];         
  };
  
  //for MPI usage
  
  void Color_Function2MPI(const Color_Function * , MPI_Color_Function &);
  
  Color_Function * MPI2Color_Function(const MPI_Color_Function &);

  std::ostream & operator<<(std::ostream &, const MPI_Color_Function &);

/*! 
    \file
    \brief this file contains the classes AMEGIC::Color_Function
           and AMEGIC::cf   
  */

  /*! 
    \class cf 
    \brief Helper class that collects the existing types of possible Color_Functions

  */

  /*! 
    \class Color_Function 
    \brief In this class the basic ingredients of \f$ SU(3)_C \f$ algebra are defined 

    In order to attach the underlying color information to a Vertex object 
    Color_Function objects are defined. A Color_Function possesses a type, 
    defined in the class AMEGIC::cf, and a number of arguments which come 
    in two types, either as an int or a char variable. Per default the integer 
    arguments of a Color_Function, called partarg, are set to '-1' and the char 
    arguments, called strarg, are set to '?'. In addition attached to a 
    Color_Function is a string expression that contains the information of the 
    Color_Function type and charakter arguments. Tab. ?? list the Color_Function 
    types defined in AMEGIC++ and there associated number of arguments the 
    \f$ SU(3)_C\f$ analoga and the corresponding string expression.
        
    <table width=300 align=center>
    <tr>
    <td>Color_Funtion</td>
    <td>args</td>
    <td>\f$ SU(3)_C \f$</td>
    <td>String()</td>
    </tr>
    <tr>
    <td>None</td>
    <td>0</td>
    <td>\f$ 1 \f$</td>
    <td>1</td>
    </tr>
    <tr>
    <td>T</td>
    <td>3</td>
    <td>\f$ T^A_{ij} \f$</td>
    <td>T[A,i,j]</td>
    </tr>
    <tr>
    <td>F</td>
    <td>3</td>
    <td>\f$ f^{ABC} \f$</td>
    <td>F[A,B,C]</td>
    </tr>
    <tr>
    <td>D</td>
    <td>2</td>
    <td>\f$ \delta_{i,j} \f$</td>
    <td>D[i,j]</td>
    </tr>
    <tr>
    <td>G</td>
    <td>2</td>
    <td>\f$ \delta_{A,B} \f$</td>
    <td>G[A,B]</td>
    </tr>
    </table>
        
    The indizes \f$ i,j \f$ denote quark color degrees of freedom and  
    \f$ A,B,C \f$ label the \f$ SU(3)_C \f$ gluon octet. Apart from its type, 
    arguments and string each Color_Function owns a pointer to a further
    Color_Function called Next, a priori this pointer is initialized as nil.

  */

   /*!
     \fn Color_Function::Color_Function() 
    \brief Constructs a Color_Function object of type cf::Unknown and all partarg 
    initialized as '-1' and all strarg as '?'. The pointer Next (also a Color_Function) 
    is initialized as nil.

  */

  /*!
     \fn Color_Function::Color_Function(cf::code _type, 
		   int _partarg0 = -1, int _partarg1 = -1, int _partarg2 = -1, 
		   char _strarg0 = '?', char _strarg1 = '?', char _strarg2 = '?') 
    \brief Can be used to contruct a Color_Function with a definite type and definite arguments.
    Again the pointer Next is initialized as nil.

  */

  /*!
    \fn void Color_Function::SetParticleArg(int a, int b,int c=-1)   
    \brief The setting of the integer arguments of a Color_Function object is 
    performed. Since the third argument is per default '-1' it is sufficient to 
    set only two arguments for a Color_Function possessing only two physical
    meaningful arguments.
  */
  
/*!
    \fn void Color_Function::SetStringArg(char a, char b, char c='?')   
    \brief The setting of the charakter arguments of a Color_Function object is 
    performed. Since the third argument is per default set to '?' it is sufficient 
    to initialize only two arguments for a Color_Function possessing only two physical
    meaningful arguments.
  */
    
  /*!
    \fn std::string Color_Function::String()
    \brief Returns the string expression of the associated Color_Function in terms 
    of the Color_Function type and its charakter arguments.
  */
    }
#endif
