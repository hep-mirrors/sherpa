#ifndef AMEGIC_Amplitude_Amplitude_Handler_H
#define AMEGIC_Amplitude_Amplitude_Handler_H

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "AMEGIC++/Main/Polarisation.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Handler.H"
#include "ATOOLS/Math/MyComplex.H"
#include "AMEGIC++/Main/Point.H"
#include "AMEGIC++/Amplitude/CFColor.H"
#include "AMEGIC++/Amplitude/Amplitude_Group.H"
#include "AMEGIC++/Main/Helicity.H"
// #include "ATOOLS/Phys/Spin_Correlation_Tensor.H"
// #include "AMEGIC++/Main/AMEGIC_SCT.H"

namespace HELICITIES { class Amplitude_Tensor; }

namespace AMEGIC {

  typedef std::map<int,CFColor*> CFC_map;
  typedef CFC_map::iterator CFC_iterator;

  //! This is the master class for the generation and calculation of an amplitude.
  class Process_Tags;
  class Amplitude_Handler : public Amplitude_Group {
    String_Handler* shand;
    CFColor              * CFCol_Matrix;   
    CFC_map                CFCol_MMatrixMap;   
    int       ngraph,namplitude,ntotal;
    int       sw_probabs;
    Complex * Mi;
    int probs;
    Single_Amplitude* firstgraph;
    std::vector<Point*> pointlist;
    bool m_print_graph;
    std::vector<Amplitude_Base*> m_ramplist;

    std::map<ATOOLS::Flavour,ATOOLS::Flavour> m_flavourmap;

//     AMEGIC_SCT* p_SCT;
    int  PropProject(Amplitude_Base*,int);
    int  CompareZfunc(Amplitude_Base*,Zfunc*,Amplitude_Base*,Zfunc*);
    void CheckEqual(Single_Amplitude*);
    void CheckEqualInGroup();
    void BuildGlobalString(Single_Amplitude*);
    void PreCluster(Single_Amplitude* firstgraph);
    void OptimizeProps(int,Single_Amplitude*);
    int  TOrder(Single_Amplitude* a);
    int  SingleCompare(Point*,Point*,double &, std::map<std::string,Complex> &);
//     void Kicker(int*,int,std::string);
    void RestoreAmplitudes(std::string path);
  public:
    Amplitude_Handler(int,ATOOLS::Flavour *,int *,Process_Tags*,
		      MODEL::Model_Base *,Topology *,int &,int &,
		      Basic_Sfuncs *,String_Handler *, bool print_graph=false, bool create_4V=true);
    ~Amplitude_Handler();
    void ConstructSignalAmplitudes(int N,ATOOLS::Flavour* fl,int* b,
				   Process_Tags* pinfo,Single_Amplitude** sglist,
				   Basic_Sfuncs* BS);
    void CompleteAmplitudes(int,ATOOLS::Flavour*,int*,Polarisation*,
			    Topology*,Basic_Sfuncs*,std::string,char=0,char=0);
    void CompleteLibAmplitudes(int N,std::string,std::string,char=0,char=0,ATOOLS::Flavour* =0);
    void StoreAmplitudeConfiguration(std::string lib);

    Complex Zvalue(int,int*);
    Complex Zvalue(int);
    Complex Zvalue(String_Handler*,int);
    Complex Zvalue(int,int,int);
    double Zvalue(Helicity* hel);
    double Differential(int ihel,int* signlist)      { return (Zvalue(ihel,signlist)).real();}
    double Differential(int ihel)                    { return (Zvalue(ihel)).real();}
    double Differential(String_Handler * sh,int ihel)  { return (Zvalue(sh,ihel)).real(); }
    double Differential(int ihel,int i,int k)        { return (Zvalue(ihel,i,k)).real();}
    //! Returns an SCT; hel=ptr to the Helicity class; nIn=Number of incoming particles.
//     ATOOLS::Spin_Correlation_Tensor* GetSpinCorrelations(Helicity* hel, size_t nIn);
    void   FillAmplitudes(HELICITIES::Amplitude_Tensor*,Helicity*,double);
    inline int GetGraphNumber() {return ngraph;}  
    inline int GetRealGraphNumber() {return namplitude;}  
    inline int GetTotalGraphNumber() {return ntotal;}

    Single_Amplitude* GetFirstGraph() {return firstgraph;}
    Point* GetPointlist(int);
    void   FillPointlist();
    int    FourVertex(int i) { return ExistFourVertex(GetPointlist(i)); } 
    bool   ExistFourVertex(Point *);
    int    CompareAmplitudes(Amplitude_Handler*,double &, std::map<std::string,Complex> &);
    std::map<ATOOLS::Flavour,ATOOLS::Flavour>& GetFlavourmap() { return m_flavourmap; }
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







