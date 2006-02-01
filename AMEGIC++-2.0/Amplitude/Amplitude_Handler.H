#ifndef Amplitude_Handler_h
#define Amplitude_Handler_h

#include "Interaction_Model_Base.H"
#include "Flavour.H"
#include "Polarisation.H"
#include "Basic_Sfuncs.H"
#include "String_Handler.H"
#include "MyComplex.H"
#include "Point.H"
#include "CFColor.H"
#include "Amplitude_Group.H"
#include "Helicity.H"
#include "Spin_Correlation_Tensor.H"
#include "AMEGIC_SCT.H"

namespace AMEGIC {
  //! This is the master class for the generation and calculation of an amplitude.
  class Process_Info;
  class Amplitude_Handler : public Amplitude_Group {
    String_Handler* shand;
    CFColor              * CFCol_Matrix;   
    int       ngraph,namplitude,ntotal;
    int       sw_probabs;
    double  * probabs;
    Complex * Mi;
    int probs;
    Single_Amplitude* firstgraph;
    std::vector<Point*> pointlist;
    bool m_print_graph;
    AMEGIC_SCT* p_SCT;
    int  PropProject(Amplitude_Base*,int);
    int  CompareZfunc(Amplitude_Base*,Zfunc*,Amplitude_Base*,Zfunc*);
    void CheckEqual(Single_Amplitude*);
    void CheckEqualInGroup();
    void BuildGlobalString(Single_Amplitude*);
    void PreCluster(Single_Amplitude* firstgraph);
    void OptimizeProps(int,Single_Amplitude*);
    int  TOrder(Single_Amplitude* a);
    int  SingleCompare(Point*,Point*,double &);
  public:
    Amplitude_Handler(int,ATOOLS::Flavour *,int *,Process_Info*,
		      Interaction_Model_Base *,Topology *,int &,int &,
		      Basic_Sfuncs *,String_Handler *, bool print_graph=false);
    ~Amplitude_Handler();
    void ConstructSignalAmplitudes(int N,ATOOLS::Flavour* fl,int* b,
				   Process_Info* pinfo,Single_Amplitude** sglist,
				   Basic_Sfuncs* BS);
    void CompleteAmplitudes(int,ATOOLS::Flavour*,int*,Polarisation*,
			    Topology*,Basic_Sfuncs*,std::string);
    void Kicker(int*,int,std::string);
    Complex Zvalue(int,int*);
    Complex Zvalue(int);
    Complex Zvalue(String_Handler*,int);
    double Zvalue(Helicity* hel);
    double Differential(int ihel,int* signlist)      { return (Zvalue(ihel,signlist)).real();}
    double Differential(int ihel)                    { return (Zvalue(ihel)).real();}
    double Differential(String_Handler * sh,int ihel)  { return (Zvalue(sh,ihel)).real(); }
    //! Returns an SCT; hel=ptr to the Helicity class; nIn=Number of incoming particles.
    ATOOLS::Spin_Correlation_Tensor* GetSpinCorrelations(Helicity* hel, size_t nIn);
    inline int GetGraphNumber() {return ngraph;}  
    inline int GetRealGraphNumber() {return namplitude;}  
    inline int GetTotalGraphNumber() {return ntotal;}

    Single_Amplitude* GetFirstGraph() {return firstgraph;}
    Point* GetPointlist(int);
    void   FillPointlist();
    void   Reset_ProbAbs();
    double Get_Probab(int);
    int    FourVertex(int i) { return ExistFourVertex(GetPointlist(i)); } 
    bool   ExistFourVertex(Point *);
    int    CompareAmplitudes(Amplitude_Handler*,double &);
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







