#ifndef Amplitude_Handler_h
#define Amplitude_Handler_h

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
    //! The alternative color matrix, for more details see CFColor.
    CFColor              * CFCol_Matrix;   
    //! Generates the Feynman diagrams.
    Amplitude_Generator  * gen;
    //! The number of amplitudes
    int       ngraph,namplitude;
    //! Not used.
    int       sw_probabs;
    //! Not used.
    double  * probabs;
    //! The list of results for the single amplitudes.
    Complex * Mi;
    //! Not used.
    int probs;
    //! The string handler supervises everything with strings. 
    String_Handler* shand;
    //! list of pointers to the pointlists
    std::vector<Point*> pointlist;
    /*!
      This method provides a mapping of propagator numbers on momentum numbers.
    */
    int  PropProject(Amplitude_Base*,int);
    /*!
      Compares whether two Zfunctions are completeley identical :
      types, number of propagators, arguments, coupling constants, 
      propagator momenta and flavours. For the handling of arguments
      and propagators, PropProject is used.
    */
    int  CompareZfunc(Amplitude_Base*,Zfunc*,Amplitude_Base*,Zfunc*);
    /*!
      Here it is checked whether two amplitudes are completely identical.
      This is used to prevent eventual double counting of amplitudes.
      It maps amplitudes on their Zfunctions which are then compared
      via CompareZfunc.
    */
    void CheckEqual(Single_Amplitude*);
    void CheckEqualInGroup();
    void BuildGlobalString(Single_Amplitude*);
    void PreCluster(Single_Amplitude* firstgraph);
    void OptimizeProps(int,Single_Amplitude*);
    int  TOrder(Single_Amplitude* a);
  public:
    //! The whole generation of the amplitudes takes place in this constructor.
    Amplitude_Handler(int,APHYTOOLS::Flavour*,int*,Polarisation*,Topology*,
		      Basic_Sfuncs*,String_Handler*,std::string);
    //! Destructing stuff.....
    ~Amplitude_Handler();
    //! Kicking amplitudes... 
    void Kicker(int*,int,std::string);
    //! Calculates the differential cross section 
    Complex Zvalue(int,int*);
    Complex Zvalue(String_Handler*,int);
    //! Calculates the differential cross section 
    double Differential(int ihel,int* signlist=0)      { return (Zvalue(ihel,signlist)).real();}
    //! Calculates the differential cross section with strings
    //    double Differential(int);
    double Differential(String_Handler * sh,int ihel)  { return (Zvalue(sh,ihel)).real(); }

    //! Returns the number of graphs.
    inline int GetGraphNumber() {return ngraph;}  
    inline int GetRealGraphNumber() {return namplitude;}  
    //! Returns the point list of a certain graph.
    Point* GetPointlist(int);
    //! Currently not used.
    void Reset_ProbAbs();
    //! Currently not used.
    double Get_Probab(int);
    /*!
      Checks whether in the full Pointlist of the argument a
      four vertex is already present.
    */
    int  FourVertex(int i) { return ExistFourVertex(GetPointlist(i)); } 
    /*!
      ExistFourVertex actually checks recursively, whether starting form
      the point there is a middle leg. If so a true is returned, if not, the
      left and the right leg are checked via calling ExistFourVertex(p->left/right).
    */
    bool  ExistFourVertex(Point *);
  };
}
#endif







