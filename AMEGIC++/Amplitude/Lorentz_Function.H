#ifndef MODEL_Interaction_Models_Lorentz_Function_H
#define MODEL_Interaction_Models_Lorentz_Function_H

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/STL_Tools.H"
#include <vector>

namespace MODEL {

  struct LF_Key {};

  class Lorentz_Function {
  private:

    Lorentz_Function(const Lorentz_Function &l): p_next(NULL) { *this=l; }

  protected:

    std::string m_type; 

    int m_partarg[4], m_permcount;

    std::vector<int*> m_permlist;
    std::vector<int>  m_signlist;
    Lorentz_Function* p_next;

    std::string Str(int a) const;

    Lorentz_Function(const std::string &type);

    virtual ~Lorentz_Function();

  public:

    virtual void Delete() = 0;

    virtual int NofIndex() const = 0;

    virtual std::string String(int shortversion=0) const = 0;

    virtual void InitPermutation();

    virtual bool CutVectors();

    virtual Lorentz_Function *GetCopy() const = 0;

    void SetParticleArg(int a=-1,int b=-1,int c=-1,int d=-1);
    void AddPermutation(int sign,int a,int b=-1,int c=-1,int d=-1);

    int  GetSign();
    int  ResetPermutation();
    int  NextPermutation();

    Lorentz_Function &operator=(const Lorentz_Function &l);

    bool operator==(const Lorentz_Function &l) const;

    inline int ParticleArg(int i) const { return m_partarg[i]; }
    inline const int* ParticleArg() const { return m_partarg; }

    inline std::string Type() const { return m_type; }

  };// end of class Lorentz_Function

  typedef ATOOLS::Getter_Function<Lorentz_Function,LF_Key> LF_Getter;

  class LF_None: public Lorentz_Function {
  public:
    LF_None();
    int NofIndex() const;
    std::string String(int shortversion) const;
    Lorentz_Function *GetCopy() const;
    static ATOOLS::AutoDelete_Vector<LF_None> s_objects;
    static Lorentz_Function *New();
    void Delete();
  };

  class LF_Pol: public Lorentz_Function {
  public:
    LF_Pol();
    int NofIndex() const;
    std::string String(int shortversion) const;
    Lorentz_Function *GetCopy() const;
    static ATOOLS::AutoDelete_Vector<LF_Pol> s_objects;
    static Lorentz_Function *New();
    void Delete();
  };

 /*! 
    \file
    \brief this file contains the classes MODEL::Lorentz_Function
           and MODEL::lf   
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
    defined in the class MODEL::lf, and a number of arguments. These arguments called 
    partarg are of the integer type. Per default they are set to '-1'. Attached to a 
    Lorentz_Function is a string expression that contains the information of the 
    Lorentz_Function type and its arguments, see Tab. ?? for a complete list of 
    Lorentz_Function types defined in MODEL++, there associated number of arguments 
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

#define DEFINE_LF_GETTER(CLASS,TAG,INFO)				\
  DECLARE_GETTER(CLASS,TAG,Lorentz_Function,LF_Key);			\
  Lorentz_Function *ATOOLS::Getter<Lorentz_Function,LF_Key,CLASS>::	\
  operator()(const LF_Key &key) const { return CLASS::New(); }		\
  void ATOOLS::Getter<Lorentz_Function,LF_Key,CLASS>::			\
  PrintInfo(std::ostream &str,const size_t width) const { str<<INFO; }

#endif

