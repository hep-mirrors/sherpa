#ifndef AMEGIC_Amplitude_Zfunc_H
#define AMEGIC_Amplitude_Zfunc_H

#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include "ATOOLS/Org/Smart_Pointer.H"

namespace AMEGIC {
  class Zfunc_Calc;

  class Pair {
  public:
    int pold,pnew;
    
    Pair(int _pold,int _pnew) : pold(_pold), pnew(_pnew) {}
  };

  struct SpinorDirection {
    int from,to;
    SpinorDirection* Next;
  };

  class CValue {
  public:
    int     * p_args;
    ATOOLS::Kabbala   m_value;
    CValue() { p_args=0; }
    CValue(int * _args, const ATOOLS::Kabbala _value) : p_args(_args) , m_value(_value) {} 
    inline bool Compare(int* args, int n){
      for(int i=0;i<n;++i) if (args[i]!=p_args[i]) return false;
      return true;
    }
  };
  typedef std::vector<CValue> Calc_List;
  typedef Calc_List::iterator CL_Iterator;

  class Zfunc {
  public:

    std::string m_type;
    int       m_narg,m_ncoupl;
    int*      p_arguments;
    Complex*  p_couplings;
    // propagators
    int       m_nprop;
    Argument* p_propagators;
    //Lorentzfunction sign
    int       m_sign;
    std::string m_str;
    Calc_List m_calclist;
    
    Zfunc*            p_equal;    
    SP(Zfunc_Calc)    p_calculator;

    Zfunc();
    Zfunc(const Zfunc& z);
    virtual ~Zfunc();
    	
    Zfunc& operator=(const Zfunc& z);
    //virtual part
    virtual void ClearCalcList();
    virtual void KillZList() {}
    virtual void ReplaceProp(std::vector<Pair>*);

    virtual Zfunc* operator[](const int i) {return this;} 

    virtual void Print();
    virtual int  GetSize() {return 1;}
    virtual char GetOp() {return 0;}

    virtual void SetSign(const int,const int) 
    { std::cerr<<"Called virtual function Zfunc::SetSign()!!!"<<std::endl; }
    virtual int  GetSign(const int i) {   
      std::cerr<<"Called virtual function Zfunc::GetSign()!!!"<<std::endl;
      return 1;
    }
    virtual int  GetSumIndex() {   
      std::cerr<<"Called virtual function Zfunc::GetSumIndex()!!!"<<std::endl;
      return 0;
    }
  };

  class Zfunc_Group: public Zfunc {
    int  m_sumindex;
    char m_op;
  public:
    std::vector<Zfunc*> m_zlist;
    std::vector<int>    m_zsigns;

    Zfunc_Group() : Zfunc() {}
    Zfunc_Group(const Zfunc& z);
    Zfunc_Group(Zfunc& z1,Zfunc& z2,int si,Pfunc_List* pl);
    
    Zfunc* operator[](const int i) { return m_zlist[i]; } 
    void Print();
    int  GetSize() { return m_zlist.size(); }
    char GetOp() { return m_op; } 
   
    void SetSign(const int pos,const int _sign) { m_zsigns[pos] = _sign; }
    int  GetSign(const int i)                   { return m_zsigns[i]; }
    int  GetSumIndex() { if (m_op=='*') return m_sumindex; else return -1; }


    void ClearCalcList();
    void KillZList();
    void ReplaceProp(std::vector<Pair>*);  
  };

  typedef std::vector<Zfunc*> Zfunc_List;
  typedef Zfunc_List::iterator Zfunc_Iterator;

  /*!
    \class CValue
    \brief Element for Zfunc::m_calclist
    
    This class contains a list of arguments and a ATOOLS::Kabbala result of a Z-function, 
    to keep track for what sets of arguments a Z function was already calculated.
  */
  /*!
    \fn CValue::Compare(int* args, int n)
    Returns true, if args contains the same arguments as the member p_args, otherwise false.
  */
  /*!
    \class Zfunc
    \brief Container for a single Z function
    
    This class contains everything neccessary for the calculation of a single Z function and
    keeps track over already calculated values.

    A Z function is part of an amplitude, usually a vertex or a group of connected vertices
    including the external legs in terms of spinors or polarization vectors. It is described
    by a type and a list of arguments, usually referring to the legs.
  */
  /*!
    \var int Zfunc::m_type
    Type of the building block to calculate this Z function.
    For a list of all available types see class zl.
  */
  /*!
    \var int Zfunc::m_ncoupl
    The number of coupling constants.
  */
  /*!
    \var int Zfunc::m_narg
    The number of arguments.
  */
  /*!
    \var int* Zfunc::p_arguments
    Pointer to an array of arguments. 
    
    Possible arguments:
    - external particles (0-9)
    - a propagator in the numbering scheme of the corresponding amplitude:
      -# fermion propagators (101-199)
      -# boson propagators (>201)
    - 99 as a dummy to indicate explicit polarization vectors

    Arguments come always in pairs (except dummy arguments for scalars).
    One pair can contain:
    - two fermion arguments from a fermion line
    - for the old external boson treatment, the two spinor arguments to construct the polarization vector
    - one boson argument and 99, that indicates that there are explicit polarization for the boson
    - two times the same argument for a spin-2 particle

    A dummy argument for a scalar is set for Z-functions, describing a vertex with scalars,
    but where the building block to calculate the Z-function does 
    not depend on the scalars momentum. These arguments are always at the end of the array.
    With the dummys, the whole information about the connection between vertices of an 
    amplitude is contained in the argument structure of Z-functions. This allows a more 
    efficient algorithm to build Super_Amplitude's.
  */
  /*!
    \var Complex* Zfunc::p_couplings;
    Pointer to an array of complex coupling constants.
  */
  /*!
    \var int Zfunc::m_nprop
    Number of propagators in the array Zfunc::p_propagators.
  */
  /*!
    \var Argument* Zfunc::p_propagators
    A pointer to an array of propagators, needed for calculation of the Z function.
    (Usually internal propagators of the building block)
  */
  /*!
    \var int Zfunc::m_sign
    An extra sign for the Z function, occuringfrom permutations in the legs of it's vertices 
    (see Lorentz_Function::GetSign).
  */
  /*!
    \var std::string Zfunc::m_str
    A unique string for the Z function within all Z functions of the process.
    Used to identify equal Z function in different amplitudes to construct a Super_Amplitude.
  */
  /*!
    \var Calc_List Zfunc::m_calclist
    Vector of CValue, to keep track over all sets of arguments, 
    for which the Z function was calculated.
  */
  /*!
    \var Zfunc* Zfunc::p_equal
    A pointer to an identical Zfunc (possibly from a different amplitude).
    It may be set to this.

    p_equal allows to use Zfunc::m_calclist of identical Zfuncs together.
  */
  /*!
    \var Zfunc_Calc* Zfunc::p_calculator
    A pointer to the corresponding calculator class.
  */

  /*!
    \fn Zfunc::Zfunc()
    Constructs an empty object.
  */
  /*!
    \fn Zfunc::Zfunc(const Zfunc& z);
    Copy constructor.
  */
  /*!
     \fn Zfunc& Zfunc::operator=(const Zfunc& z);
     Copies a Zfunc including all elements of the Zfunc::p_arguments,
     Zfunc::p_couplings and Zfunc::p_propagators.
  */
  /*!
    \fn Zfunc::ClearCalcList();
    Deletes all entries in Zfunc::m_calclist.
  */
  /*!
    \fn Zfunc::ReplaceProp(std::vector<Pair>*);
    Changes the labels of propagators (Zfunc::p_propagators) and 
    arguments (Zfunc::p_arguments) from Pair::pold to Pair::new.
  */
  /*!
    \fn virtual Zfunc* Zfunc::operator[](const int i) 
    Returns the pointer this.
  */
  /*!
    \fn virtual void Zfunc::Print()
    Prints the Zfunc with type, arguments and coupling constants.
  */
  /*!
    \class Zfunc_Group
    \brief Contains a list a Zfunc, connected by an operation + or *.
    
    This class is derived from Zfunc and contains a list of Z functions.
    When an Single_Amplitude or Super_Amplitude is calculated in Single_Amplitude_Base,
    it's Zfunc structure is transformed into a recursive structure of Zfunc_Group's, 
    connected by the operations + and *.
    
    - '+' arrise when several amplitudes, containing identical parts are merged together to
    a Super_Amplitude. 
    The '+'-group can have two or more sub Z functions, that all have 
    the same argument list (up to the sequence).
    - '*' combines two Z functions (or groups) that are connected by a propagator. The 
    sub Z functions have only the argument referring to this propagator common. The 
    '*'-operation here goes along with the summation over spinor/antispinor modes for
    fermion propagators or polarization vectors for boson propagators.

    Zfunc_Group has no defined type and calculator. The calculation is done 
    in Single_Amplitude_Base::SingleZGoupValue by 
    performing the operator and the calculators of the sub Z functions.
  */
  /*!
    \var int  Zfunc_Group::m_sumindex
    For a Zfunc_Group with operation '*':
    Label of the connecting propagator to perform the sum over all propagating modes.
  */
  /*!
    \var char Zfunc_Group::m_op
    Operation to be performed to combine the sub Z functions.
    Possible values are '+' or '*'.
  */
  /*!
    \var std::vector<Zfunc*> Zfunc_Group::m_zlist;
    A vector of all the sub Z functions.
  */
  /*!
    \var std::vector<int>    Zfunc_Group::m_zsigns
    For a Zfunc_Group with operation '+':
    A sign for each sub Z function.
  */
  /*!
    \fn Zfunc_Group::Zfunc_Group()
    Constructs an empty object.
  */
  /*!
    \fn Zfunc_Group::Zfunc_Group(const Zfunc& z)
    Constructor for a group with operation '+'.

    Copies the argument list of the Zfunc z to the group. 
    This can be done with any Z function of the latter group, since they must
    have identical arguments.

    Zfunc_Group::m_zlist must be filled separatly.
  */
  /*!
    \fn Zfunc_Group::Zfunc_Group(Zfunc& z1,Zfunc& z2,int si,Pfunc_List* pl);
     Constructor for a group with operation '*'.
     
     The new array of arguments contains both arrays of the sub Z functions, except
     the label for the connecting propagator. This is now becoming an internal 
     argument, to be kept in Zfunc_Group::m_sumindex.

     Both propagator arrays are copied into the group array and are removed from 
     the sub Z functions.

     All coupling constants are copied to the group.
  */
  /*!
    \fn Zfunc* Zfunc_Group::operator[](const int i)
    Returns the Zfunc #i from Zfunc_Group::m_zlist.
  */
  /*!
    \fn void Zfunc_Group::Print()
    Prints the Zfunc_Group with its arguments and sub Z functions.
  */
  /*!
    \fn int  Zfunc_Group::GetSize()
    Returns the number of cointained sub Z functions.
  */
  /*!
    \fn char Zfunc_Group::GetOp()
    Returns the operator, + or *.
  */
  /*!
    \fn void Zfunc_Group::SetSign(const int pos,const int _sign)
    Set the sign for sub Z function at position pos. (only '+'-groups)
  */
  /*!
    \fn int  Zfunc_Group::GetSign(const int i)
    Returns the sign for the sub Z function at position i. (only '+'-groups)
  */
  /*!
    \fn int  Zfunc_Group::GetSumIndex()
    Returns the index of the connecting propagator of two Z functions for '*'-groups.

  */
  /*!
    \fn void Zfunc_Group::ClearCalcList()
    Deletes all entries in m_calclist including all lists of its sub Z functions.
    
   */
  /*!
    \fn void Zfunc_Group::ReplaceProp(std::vector<Pair>*)
    See Zfunc::ReplaceProp.
  */ 
}

#endif
