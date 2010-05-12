#ifndef AMEGIC_Amplitude_Single_Amplitude_Base_H
#define AMEGIC_Amplitude_Single_Amplitude_Base_H

#include "AMEGIC++/Amplitude/Amplitude_Base.H"
#include "ATOOLS/Math/Kabbala.H"

namespace AMEGIC {

  class Argument;

  class Single_Amplitude_Base: public virtual Amplitude_Base {
  protected:
    int amplnumber;

    Zfunc_List *zlist;
    Pfunc_List plist;
    int * b;
    int   N;

    String_Handler     * shand;
    Basic_Sfuncs       * BS;
    ATOOLS::Flavour * fl;

    virtual ATOOLS::Kabbala SingleZvalue(Zfunc*,std::vector<int>*,std::vector<int>*,int last=0);
    ATOOLS::Kabbala SingleZGroupvalue(Zfunc* z, std::vector<int>* iz, std::vector<int>* iargs,int last=0);
    ATOOLS::Kabbala SingleZvalueTensor(Zfunc* z,std::vector<int>* iz, std::vector<int>* iargs,int);
    int     FillArgs(Zfunc* z,int* args,std::vector<int>*,std::vector<int>*);
    void    SetLoopVar(std::vector<int>&,std::vector<std::vector<int> >&);

    ATOOLS::Kabbala SingleMassTerms(int iz,int iarg);
    ATOOLS::Flavour* GetPflav(int pn);
    ATOOLS::Kabbala GetProp(Zfunc* z);
    void    GroupZfuncs();
  public:
    int sign;
    int m_oqed,m_oqcd;
    Single_Amplitude_Base(int*,int,Basic_Sfuncs*,ATOOLS::Flavour* ,String_Handler*); 
    Single_Amplitude_Base(String_Handler*,int); 
    virtual ~Single_Amplitude_Base();

    Complex Zvalue(int);
    Complex Zvalue(int,int*);
    Complex Zvalue(String_Handler * sh,int ihel);

    Zfunc_List* GetZlist();
    Pfunc_List* GetPlist();
    int GetSign();
    void SetSign(int);

    void ClearCalcList();
    void KillZList();
    void PrintGraph(); 
    void SetNumber(int& i);
    Amplitude_Base* GetAmplitude(const int n);
    const int GetNumber() { return amplnumber; }
    
    void    SetOrderQED(int);
    void    SetOrderQCD(int);
    virtual int     GetOrderQED();
    virtual int     GetOrderQCD();

  };

  /*!
    \class Single_Amplitude_Base
    \brief Numerical calculation and string generation for one amplitude.

    The numerical calculation and string generation for a Single_Amplitude as well as
    for a Super_Amplitude is performed in this class.

     
  */
  /*!
    \var Zfunc_List Single_Amplitude_Base::zlist
    A list of Z functions in the amplitude; Calling 
    Single_Amplitude_Base::GroupZfuncs transforms it into a Zfunc_Group,
    that recursivley contains all original Z functions.
  */
  /*!
    \var Pfunc_List Single_Amplitude_Base::plist
    A list of all propagators in the amplitude.
  */
  /*!
    \fn ATOOLS::Kabbala Single_Amplitude_Base::SingleZvalue(Zfunc* z,std::vector<int>* iz,std::vector<int>* iargs,int last=0)
    Calculates the Z function z. 
    
    The vector iz contains a list of labels for propagators or external particles, 
    while iargs is a list of corresponding helicity/polarization states.

    When SingleZvalue is first called from Single_Amplitude_Base::Zvalue, iz is 
    just a list 0,1,2,...,n with n the number of external particles.

    This method may be also called recursively from Single_Amplitude_Base::SingleZGroupvalue
    or Single_Amplitude_Base::SingleZvalueTensor.

    - Single_Amplitude_Base::FillArgs generates the argument array for the 
    calculator of the Zfunc (Zfunc_Calc).
    - In the case, the calculated process has external spin-2 particles, 
    Single_Amplitude_Base::FillArgs returns their position in iz. The corresponding
    dummy polarization modes (see class mt) in iargs will be replaced
    in Single_Amplitude_Base::SingleZvalueTensor by sums of outer products
    of two polarization vectors.
    - The Zfunc::m_calclist is checked, if the Z function was already calculated
    for the current set of arguments. If yes, the result will be returned directly.
    - Calculation:
      - If the Zfunc is not elementary, but a Zfunc_Group, it is calculated
      in Single_Amplitude_Base::SingleZGroupvalue.
      - For an elementary Zfunc its calculator (Zfunc_Calc* Zfunc::p_calculator)
      is started.
    - The result is saved in Zfunc::m_calclist and returned.
  */
  /*!
    \fn ATOOLS::Kabbala Single_Amplitude_Base::SingleZGroupvalue(Zfunc* z, std::vector<int>* iz, std::vector<int>* iargs,int last=0)
    Calculates a Zfunc_Group. 

    - For a group with the operator '+' it's simply the sum of all sub Z functions, 
    calculated by Single_Amplitude_Base::SingleZvalue with the same set of arguments.

    - For a group with the operator '*' the product of the two Z functions is calculated.
    Depending on the propagator connecting these functions also a sum over the 
    helicity/polarization modes has to be performed in order to reproduce the 
    propagator's Dirac/Lorentz structure. The list of corresponding arguments is determined
    in Single_Amplitude_Base::SetLoopVar.    
  */
  /*!
    \fn ATOOLS::Kabbala Single_Amplitude_Base::SingleZvalueTensor(Zfunc* z,std::vector<int>* iz, std::vector<int>* iargs,int)
  */
  /*!
    \fn int Single_Amplitude_Base::FillArgs(Zfunc* z,int* args,std::vector<int>*,std::vector<int>*)
    Generates the argument array args for the calculator of the Zfunc (Zfunc_Calc).

    For every element of Zfunc::p_arguments args contains a pair of elements:
    - external particles: number of the particle; helicity or polarization
    - boson propagators: propagator label; polarization (0 for scalars)
    - fermion propagators:
      - outer particle version: number of an external particle; helicity
      - inner particle version: propagator label; helicity
      - the spinor type is encoded in an extra \f$-\f$ sign in the arguments 
  */
  /*!
    \fn void Single_Amplitude_Base::SetLoopVar(std::vector<int>& iz ,std::vector<std::vector<int> >& iarg)
    Determines a list of arguments for the decomposition of a propagator.
  */
  /*!
    \fn ATOOLS::Kabbala Single_Amplitude_Base::SingleMassTerms(int iz,int iarg)
    Calculates the mass term for the decomposition of fermionic propagators.
  */
  /*!
    \fn ATOOLS::Flavour* Single_Amplitude_Base::GetPflav(int pn)
    Returns the flavour of the propagator with label pn.
  */
  /*!
    \fn ATOOLS::Kabbala Single_Amplitude_Base::GetProp(Zfunc*)
    Returns the product of the propagators in the Zfunc.
    If the Zfunc is elementary only propagators with the flag Pfunc::on are
    returned. This avoids double countings.
   */
  /*!
   \fn void Single_Amplitude_Base::GroupZfuncs();
   Rearranges the the Z funtions of the amplitude from a linear list to
   a recursive structure of Zfunc_Group. 

   Within the helicity method an amplitude 
   (helicities/polarizations of external particles are fixed) is a product of
   Z functions, with a sum over spinor/helicity or polarization modes of the 
   propagator.

   The procedure in this method corresponds to a factoring out in order 
   to reduce the number of complex multiplications.
  */
  /*!
   \var int Single_Amplitude_Base::sign
   An overall sign for the amplitude.
  */
  /*!
    \fn Complex Single_Amplitude_Base::Zvalue(int,int*)
    Performs the calculation of an amplitude.

    - If Single_Amplitude_Base::zlist contains more than one Zfunc,
    Single_Amplitude_Base::GroupZfuncs is called.
    - The remaining Zfunc in the list has only arguments of external particles
    and is calculated using Single_Amplitude_Base::SingleZvalue.
    - The result is multiplied with all propagator factors, 
    that were factored out.
  */
  /*!
    \fn void Single_Amplitude_Base::ClearCalcList()
    Executes Zfunc::ClearCalcList() for all Z functions in the amplitude. 
  */
 }
#endif












