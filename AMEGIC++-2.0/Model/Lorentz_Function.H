#ifndef Lorentz_Function_h
#define Lorentz_Function_h

#include "MyStrStream.H"
//#include <strstream>
//#include <sstream>
#include <vector>
#include "MyString.H"

namespace AMEGIC {

  class lf {
  public:
    enum code {None   = 0, 
	       Gamma  = 1, 
	       Gab    = 2, 
	       Gauge3 = 3, 
	       Gauge4 = 4, 
	       Gluon4 = 5, 
	       SSV    = 6, 
	       SSS    = 7, 
	       FFS    = 8, 
	       Pol    =10,
	       VVSS   =11,//the generic interaction of 2 vectors and 2 scalars (Gab)
	       SSSS   =12,//the generic interaction of 4 scalars
	       FFT    =60,
	       VVT    =61,
	       SST    =62,
	       FFVT   =63,
	       VVVT   =64,
	       SSST   =65,
               FFGS   =70,
	       VVGS   =71,
	       SSGS   =72,
	       FFVGS  =73,
	       Unknown=99}; 
  };
  
  class Lorentz_Function {
    int m_permcount;
    std::vector<int*> m_permlist;
    std::vector<int>  m_signlist;

    std::string Str(int a) const
    {
      MyStrStream sstr;
      sstr<<partarg[a];
      std::string help;
      sstr>>help;
      return help;
    } 
  public:
    Lorentz_Function() {
      type = lf::Unknown;
      Next = 0;
      for (short int i=0;i<4;i++) partarg[i] = -1;
    }

    Lorentz_Function(lf::code _type) :type(_type) {
      Next = 0;
      for (short int i=0;i<4;i++) partarg[i] = -1;
    }
    
    lf::code  type; 
    int  partarg[4];
    Lorentz_Function* Next;
    int NofIndex() {
      switch (type) {
      case lf::Gamma : return 1;
      case lf::Gab   : return 2;
      case lf::Gauge3: return 3;
      case lf::Gauge4: return 4;
      case lf::Gluon4: return 4;
      case lf::SSV   : return 3;	
      case lf::Pol   : return 1;	
      case lf::VVSS  : return 2;	
      case lf::FFT   : return 1;
      case lf::VVT   : return 3;
      case lf::SST   : return 3;
      case lf::FFVT  : return 2;
      case lf::VVVT  : return 4;
      case lf::SSST  : return 1;
      case lf::VVGS  : return 3;
      case lf::FFVGS : return 1;
      case lf::SSGS  : return 2;
      default: return 0;
      }
      //default (includes FFS, SSS and SSSS )
      return 0;
    }

    void SetParticleArg(int a=-1, int b=-1,int c=-1,int d=-1) {
      partarg[0] = a;partarg[1] = b;
      partarg[2] = c;partarg[3] = d;
    }

    int  GetSign();
    void AddPermutation(int,int,int,int,int);
    void InitPermutation();
    int  ResetPermutation();
    int  NextPermutation();

    std::string String(int shortversion=0) const {
      if (type==lf::SSS)  return std::string("1");
      if (type==lf::FFS)  return std::string("1");
      if (type==lf::SSSS) return std::string("1");
      std::string help;
      switch (type) {
      case lf::Gamma:  
	// Gam[0]
	help = std::string("Gam[") + Str(0) + std::string("]");break;
      case lf::Pol:  
	// Eps[0]
	help = std::string("Eps[") + Str(0) + std::string("]");break;
      case lf::Gab:   
	// G[0,1]
	help = std::string("G[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::VVSS:   
	// G(2V2S)[0,1]
	help = std::string("G(2V2S)[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::Gauge3: 
	// (P[0,2]-P[1,2])*G(0,1)+(P[1,0]-P[2,0])*G(1,2)+(P[2,1]-P[0,1])*G(2,0)
	if (shortversion) {
	  help += std::string("V3[") + Str(0) + std::string(",") + 
	                          Str(1) + std::string(",") + 
                                  Str(2) + std::string("]");
	}
	else {
	  help  = std::string("(P[") + Str(0) + std::string(",") + Str(2) + std::string("]-");
	  help += std::string("P[")  + Str(1) + std::string(",") + Str(2) + std::string("])*");
	  help += std::string("G[")  + Str(0) + std::string(",") + Str(1) + std::string("]");
	  
	  help += std::string("+");
	  
	  help += std::string("(P[") + Str(1) + std::string(",") + Str(0) + std::string("]-");
	  help += std::string("P[")  + Str(2) + std::string(",") + Str(0) + std::string("])*");
	  help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("]");
	  
	  help += std::string("+");
	  
	  help += std::string("(P[") + Str(2) + std::string(",") + Str(1) + std::string("]-");
	  help += std::string("P[")  + Str(0) + std::string(",") + Str(1) + std::string("])*");
	  help += std::string("G[")  + Str(2) + std::string(",") + Str(0) + std::string("]");
	}
	break;
      case lf::SSV:
	//P[0,2]-P[1,2]
	help = std::string("P[") + Str(0) + std::string(",") + Str(2) +std::string("]-"); 
	help += std::string("P[") + Str(1) + std::string(",") + Str(2) +std::string("]");
	break;
      case lf::Gauge4: 
	//(2G(0,1)*G(2,3)-G(0,2)*G(1,3)-G(0,3)*G(1,2))
	help  = std::string("(2*G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
	help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
	help += std::string("G[")  + Str(0) + std::string(",") + Str(2) + std::string("]*");
	help += std::string("G[")  + Str(1) + std::string(",") + Str(3) + std::string("]-");
	help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
	help += std::string("G[")  + Str(1) + std::string(",") + Str(2) + std::string("])");
	break;
      case lf::Gluon4:
	//G(0,1)*G(2,3)-G(0,3)*G(2,1)
	if (shortversion) {
	  help += std::string("G4[") + Str(0) + std::string(",") + 
	                          Str(1) + std::string(",") + 
	                          Str(2) + std::string(",") + 
                                  Str(3) + std::string("]");
	}
	else {
	  help  = std::string("(G[")  + Str(0) + std::string(",") + Str(1) + std::string("]*");
	  help += std::string("G[")  + Str(2) + std::string(",") + Str(3) + std::string("]-");
	  help += std::string("G[")  + Str(0) + std::string(",") + Str(3) + std::string("]*");
	  help += std::string("G[")  + Str(2) + std::string(",") + Str(1) + std::string("])");
	}
	break;
       case lf::FFT:
	help = std::string("FFT[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::FFVT:
	help = std::string("FFVT[") + Str(0) + std::string(",") + Str(1) + std::string("]");break;
      case lf::FFVGS:
	help = std::string("FFVGS[") + Str(0) + std::string("]");break;
      case lf::VVT:
	help = std::string("VVT[") + Str(0) + std::string(",") + Str(1) 
	     + std::string(",") + Str(2) + std::string("]");break;    
      case lf::VVGS:
	help = std::string("VVGS[") + Str(0) + std::string(",") + Str(1) 
	     + std::string(",") + Str(2) + std::string("]");break;
      }     
      return help;
    }
  };

  // Format of a Lorentz_Function to transport via MPI:
  struct MPI_Lorentz_Function {
    
    int m_type;
    int    m_partarg[4];         

   };

  //for MPI usage
  
  void Lorentz_Function2MPI(const Lorentz_Function * , MPI_Lorentz_Function &);
  
  Lorentz_Function * MPI2Lorentz_Function(const MPI_Lorentz_Function &);
  
  std::ostream & operator<<(std::ostream &, const MPI_Lorentz_Function &);

 /*! 
    \file
    \brief this file contains the classes AMEGIC::Lorentz_Function
           and AMEGIC::lf   
  */

  /*! 
    \class lf 
    \brief Helper class to collect the existing types of possible Lorentz_Functions

    
  */
 
  /*! 
    \class Lorentz_Function 
    \brief In this class basic Lorentz structures as present in the SM, MSSM and LED
           are defined. According to these Lorentz structures generated amplitudes can
	   uniquely be translated into helicity amplitudes.

    In order to attach the Lorentz structure to a Single_Vertex object 
    Lorentz_Function objects are defined. A Lorentz_Function possesses a type, 
    defined in the class AMEGIC::lf, and a number of arguments. These arguments called 
    partarg are of the integer type. Per default they are set to '-1'. Attached to a 
    Lorentz_Function is a string expression that contains the information of the 
    Lorentz_Function type and its arguments, see Tab. ?? for a complete list of 
    Lorentz_Function types defined in AMEGIC++, there associated number of arguments 
    and the corresponding string expression. 
    
    <table width=450 align=center>
    <tr>
    <td>Lorentz_Function</td>
    <td>args</td>
    <td>Lorentz structure </td>
    <td>string</td>
    </tr>
    <tr>
    <td>Pol</td>
    <td>1</td>
    <td>\f$ \epsilon^{\mu} \f$</td>
    <td>EPS[0]</td>
    </tr>
    <tr>
    <td>SSS</td>
    <td>0</td>
    <td>\f$ 1 \f$</td>
    <td>1</td>
    </tr>
    <tr>
    <td>FFS</td>
    <td>0</td>
    <td>\f$ 1 \f$</td>
    <td>1</td>
    </tr>
    <tr>
    <td>SSV</td>
    <td>3</td>
    <td>\f$ (p_0-p_1)^{\mu} \f$</td>
    <td>P[0,2]-P[1,2]</td>
    </tr>
    <tr>
    <td>Gamma</td>
    <td>1</td>
    <td>\f$ \gamma^{\mu} \f$</td>
    <td>Gam[0]</td>
    </tr>
    <tr>
    <td>Gab</td>
    <td>2</td>
    <td>\f$ g^{\mu\nu} \f$</td>
    <td>G[0,1]</td>
    </tr>
    <tr>
    <td>Gauge3</td>
    <td>3</td>
    <td>\f$ (p_0-p_1)^{\mu}g^{\nu\lambda} + (p_1-p_2)^{\nu} g^{\lambda\mu}+ (p_2-p_0)^{\lambda}g^{\mu\nu} \f$</td>
    <td>(P[0,2]-P[1,2])*G[0,1]+(P[1,0]-P[2,0])*G[1,2]+(P[2,1]-P[0,1])*G[2,0]</td>
    </tr>
    <tr>
    <td>Gauge4</td>
    <td>4</td>
    <td>\f$ 2g^{\alpha\beta}g^{\mu\nu} - g^{\alpha\mu}g^{\beta\nu} - g^{\alpha\nu}g^{\beta\mu}\f$</td>
    <td>2*G(0,1)*G(2,3) - G(0,2)*G(1,3) -  G(0,3)*G(1,2)</td>
    </tr>
    <tr>
    <td>Gluon4</td>
    <td>4</td>
    <td>\f$ g^{\alpha\beta}g^{\mu\nu} - g^{\alpha\nu}g^{\beta\mu}\f$</td>
    <td>G(0,1)*G(2,3) - G(0,3)*G(1,2)</td>
    </tr>
    <tr>
    <td>VVSS</td>
    <td>2</td>
    <td>\f$ g^{\mu\nu}\f$</td>
    <td>G(2V2S)[0,1]</td>
    </tr>
    <tr>
    <td>SSSS</td>
    <td>0</td>
    <td>\f$ 1\f$</td>
    <td>1</td>
    </tr>
    <tr>
    <td>FFT</td>
    <td>1</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>VVT</td>
    <td>3</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>SST</td>
    <td>3</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>FFVT</td>
    <td>2</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>VVVT</td>
    <td>4</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>SSST</td>
    <td>1</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>FFGS</td>
    <td>0</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>VVGS</td>
    <td>3</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>SSGS</td>
    <td>2</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    <tr>
    <td>FFVGS</td>
    <td>1</td>
    <td>\f$ ??\f$</td>
    <td>??</td>
    </tr>
    </table>

    In addition to its type and arguments a Lorentz_Function owns a pointer to a further
    Lorentz_Function called Next, a priori this pointer is initialized as nil.
  */
 
   /*!
     \fn Lorentz_Function::Lorentz_Function() 
    \brief Constructs a Lorentz_Function object of type lf::Unknown and all partarg 
    initialized as '-1'. The Next pointer (which represents a Lorentz_Functionas well) 
    is initialized as nil.

  */

  /*!
     \fn Lorentz_Function::Lorentz_Function(lf::code _type) 
     \brief Contructs a Lorentz_Function object of type _type. The arguments are initialized 
     as '-1' and the pointer Next is set to nil.
  */


  /*!
    \fn std::string Lorentz_Function::String(int shortversion=0)
    \brief Returns the string expression of the associated Lorentz_Function in terms 
    of the underlying Lorentz structure or a shorthand notation. 
  */

  /*!
    \fn int Lorentz_Function::NofIndex()
    \brief Returns the number of arguments the Lorentz_Function owns.
  */

  /*!
    \fn void Lorentz_Function::SetParticleArg(int a=-1, int b=-1,int c=-1,int d=-1)
    \brief The setting of the integer arguments of a Lorentz_Function object is 
    performed. Since all arguments are per default set to '-1' it is sufficient to 
    specify only the number of physical meaningful arguments even this number is smaller
    than four.
  */

  /*!
    \fn int Lorentz_Function::GetSign()
    \brief Returns the sign according to a specific permutation of the arguments of 
    a Lorentz_Function.
  */
  
  /*!
    \fn void Lorentz_Function::AddPermutation(int,int,int,int,int);
    \brief A new permution can be added to the list of permutations associated
    to a Lorentz_Function.

    The first argument determines the sign of the permutation and therefore
    can have the two values '\f$\pm 1\f$'. The remaining arguments determine the 
    exchange of the original order of arguments of the Lorentz_Function under consideration.
    Taking a Lorentz_Function with four arguments (0,1,2,3) the call AddPermutation(-1,2,1,0,3) 
    would result in the argument tuple (2,1,0,3) and the sign '-1' is added to the m_signlist. 

  */
    
  /*!
    \fn void Lorentz_Function::InitPermutation();
    \brief The predefined permutations of all Lorentz_Function objets are initialized.

    After clearing the lists m_permlist and m_signlist of the Lorentz_Function the existing argument 
    permutations are added to m_permlist and the corresponding signs to m_signlist. The counter 
    m_percount is initialized as zero.
  */
  /*!
    \fn int Lorentz_Function::ResetPermutation();
    \brief Resets the arguments of a Lorentz_Function to the original order, resp. (0,1,2,3). 
    The counter m_percount is set to zero.
  */
  /*!
    \fn int Lorentz_Function::NextPermutation();
    \brief This method trys to set up a new permutation of the arguments be going through 
    the m_permlist of the object where the current position is determined by the value of
    m_permcount.

    For Lorentz_Functions possessing less than two arguments clearly there exists no nontrivial
    permutation and therefore the method returns a zero. For all other cases the method rearanges 
    the arguments of the Lorentz_Function according to a permutation out list of m_permlist and 
    increases the value of m_permcount by one. For each call of the method one new permutation is 
    set up and in case the end of list is reached a zero is returned.
  */
  
  /*!
    \var int Lorentz_Function::m_permcount;
    \brief Variable that stores the current position in the list of permutations m_permlist.
  */
  
   /*!
     \var std::vector<int*> Lorentz_Function::m_permlist;
     \brief List of allowed argument permutations belonging to a Lorentz_Function.
   */
   /*!
     \var std::vector<int> Lorentz_Function::m_signlist;
     \brief Associated list of signs corresponding to the permutations in m_permlist.
   */
}
#endif

