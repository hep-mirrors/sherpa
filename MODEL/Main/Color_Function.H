#ifndef MODEL__Main__Color_Function_H
#define MODEL__Main__Color_Function_H

#include <string>

namespace MODEL {

  class cf {
  public:
    enum code {
      T=0, 
      F=1, 
      D=2, 
      None=3, 
      G=4,  
      UFO=5,
      Unknown=99}; 
  };

  std::ostream &operator<<(std::ostream &str,const cf::code &c);

  class Color_Function {
  public:
    cf::code  m_type; // Number of args, multi gluons etc.
    int  m_partarg[3];
    char m_strarg[3];
    std::string m_string;
    Color_Function* p_next;
  public:
    Color_Function()  {
      m_type = cf::Unknown;
      p_next = 0;
      for (short int i=0;i<3;++i) m_partarg[i] = -1;
      for (short int i=0;i<3;++i) m_strarg[i] = '?';
    }
    
    Color_Function(cf::code _type, 
		   int _partarg0 = -1, int _partarg1 = -1, int _partarg2 = -1, 
		   Color_Function* n=NULL) 
      : m_type(_type), p_next(n) 
    {       
      m_partarg[0] = _partarg0;
      m_partarg[1] = _partarg1;
      m_partarg[2] = _partarg2;
      m_strarg[0]  = _partarg0>0?_partarg0+47:52;
      m_strarg[1]  = _partarg1>0?_partarg1+47:52;
      m_strarg[2]  = _partarg2>0?_partarg2+47:52;
    }

    Color_Function(const Color_Function & c) : p_next(NULL)
    {
      *this=c;
    }

    Color_Function & operator=(const Color_Function & c);
    ~Color_Function();

    bool operator==(const Color_Function &c) const;
    
    void SetParticleArg(int a, int b,int c=-1) {
      m_partarg[0] = a;m_partarg[1] = b;
      if (m_type!=cf::D && m_type!=cf::G) m_partarg[2] = c;
    }

    void SetStringArg(char a, char b, char c='?') {
      m_strarg[0] = a;m_strarg[1] = b;
      if (m_type!=cf::D && m_type!=cf::G) m_strarg[2] = c;
    }

    int ParticleArg(int i) const
    { 
      return m_partarg[i];
    }
    char StringArg(int i) const
    { 
      return m_strarg[i];
    }
    cf::code Type() const
    {
      return m_type;
    }
    Color_Function * Next() { 
      return p_next; 
    }
    void Conjugate() 
    {
      if (Type()!=cf::T) return;

      int help   = m_partarg[1];
      m_partarg[1] = m_partarg[2];
      m_partarg[2] = help;

      char help2 = m_strarg[1];
      m_strarg[1] = m_strarg[2];
      m_strarg[2] = help2;  
    }
    
    void Replace(int hit,int prop) {
      for (short int i=0;i<3;i++) {
	if ((Type()==cf::D || Type()==cf::G) && i==2) break;    
	if (ParticleArg(i)==hit) m_partarg[i] = prop;
      }
    }
    void Replace(int hit,char repl) {
      for (short int j=0;j<3;j++) {
	if ((Type()==cf::D || Type()==cf::G) && j==2) break;
	if (ParticleArg(j)==hit) m_strarg[j] = repl;
      }
    }
    void Append(Color_Function * last) 
    {
      Color_Function* help = this;
      while (help->p_next) help = help->p_next;
      help->p_next = last;  
    }
    Color_Function * Erase(Color_Function * granny) 
    {
      Color_Function * mo = p_next;
      if (granny)  granny->p_next=mo;
      p_next=0;
      delete this;
      return mo;
    }

    std::string String() const;
    std::string PID() const;   
    std::string FullString() const;
  };  

/*! 
    \file
    \brief this file contains the classes MODEL::Color_Function
           and MODEL::cf   
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
    defined in the class MODEL::cf, and a number of arguments which come 
    in two types, either as an int or a char variable. Per default the integer 
    arguments of a Color_Function, called partarg, are set to '-1' and the char 
    arguments, called strarg, are set to '?'. In addition attached to a 
    Color_Function is a string expression that contains the information of the 
    Color_Function type and charakter arguments. Tab. ?? list the Color_Function 
    types defined in MODEL++ and there associated number of arguments the 
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
