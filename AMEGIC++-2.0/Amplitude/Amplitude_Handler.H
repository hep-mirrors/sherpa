#ifndef Amplitude_Handler_h
#define Amplitude_Handler_h

#include "Interaction_Model_Base.H"
#include "Flavour.H"
#include "Polarisation.H"
#include "Basic_Sfuncs.H"
#include "String_Handler.H"
#include "MyComplex.H"
#include "Point.H"
#include "Amplitude_Generator.H"
#include "CFColor.H"
#include "Amplitude_Group.H"

namespace AMEGIC {
  //! This is the master class for the generation and calculation of an amplitude.
  class Amplitude_Handler : public Amplitude_Group {
    CFColor              * CFCol_Matrix;   
    Amplitude_Generator  * gen;
    int       ngraph,namplitude;
    int       sw_probabs;
    double  * probabs;
    Complex * Mi;
    int probs;
    String_Handler* shand;
    std::vector<Point*> pointlist;
    int  PropProject(Amplitude_Base*,int);
    int  CompareZfunc(Amplitude_Base*,Zfunc*,Amplitude_Base*,Zfunc*);
    void CheckEqual(Single_Amplitude*);
    void CheckEqualInGroup();
    void BuildGlobalString(Single_Amplitude*);
    void PreCluster(Single_Amplitude* firstgraph);
    void OptimizeProps(int,Single_Amplitude*);
    int  TOrder(Single_Amplitude* a);
  public:
    Amplitude_Handler(int,ATOOLS::Flavour *,int *,Polarisation *,
		      Interaction_Model_Base *,Topology *,int,int,
		      Basic_Sfuncs *,String_Handler *,std::string);
    ~Amplitude_Handler();
    void Kicker(int*,int,std::string);
    Complex Zvalue(int,int*);
    Complex Zvalue(String_Handler*,int);
    double Differential(int ihel,int* signlist=0)      { return (Zvalue(ihel,signlist)).real();}
    double Differential(String_Handler * sh,int ihel)  { return (Zvalue(sh,ihel)).real(); }

    inline int GetGraphNumber() {return ngraph;}  
    inline int GetRealGraphNumber() {return namplitude;}  
    Point* GetPointlist(int);
    void Reset_ProbAbs();
    double Get_Probab(int);
    int  FourVertex(int i) { return ExistFourVertex(GetPointlist(i)); } 
    bool  ExistFourVertex(Point *);
  };

    //! The alternative color matrix, for more details see CFColor.
    //! Generates the Feynman diagrams.
    //! The number of amplitudes
    //! Not used.
    //! Not used.
    //! The list of results for the single amplitudes.
    //! Not used.
    //! The string handler supervises everything with strings. 
    //! list of pointers to the pointlists
    /*!
      This method provides a mapping of propagator numbers on momentum numbers.
    */
    /*!
      Compares whether two Zfunctions are completeley identical :
      types, number of propagators, arguments, coupling constants, 
      propagator momenta and flavours. For the handling of arguments
      and propagators, PropProject is used.
    */
    /*!
      Here it is checked whether two amplitudes are completely identical.
      This is used to prevent eventual double counting of amplitudes.
      It maps amplitudes on their Zfunctions which are then compared
      via CompareZfunc.
    */
    //! The whole generation of the amplitudes takes place in this constructor.
    //! Destructing stuff.....
    //! Kicking amplitudes... 
    //! Calculates the differential cross section 
    //! Calculates the differential cross section 
    //! Calculates the differential cross section with strings
    //    double Differential(int);
    //! Returns the number of graphs.
    //! Returns the point list of a certain graph.
    //! Currently not used.
    //! Currently not used.
    /*!
      Checks whether in the full Pointlist of the argument a
      four vertex is already present.
    */
    /*!
      ExistFourVertex actually checks recursively, whether starting form
      the point there is a middle leg. If so a true is returned, if not, the
      left and the right leg are checked via calling ExistFourVertex(p->left/right).
    */

}
#endif







