#ifndef Amplitude_Generator_h
#define Amplitude_Generator_h

#include "Flavour.H"

#include "Interaction_Model_Base.H"
#include "Vertex.H"
#include "Topology.H"
#include "Single_Amplitude.H"
#include "Basic_Sfuncs.H"
#include "String_Handler.H"

#include <map>
#include <vector>

//#define _USE_MPI_
//#ifdef _USE_MPI_
//#include <mpi++.h>
//#endif

namespace AMEGIC {
  class Pre_Amplitude {
  public:
    Point* p;
    int on;
    int top;
    int perm;
    Pre_Amplitude():p(0),on(1),top(0),perm(0) {}
    Pre_Amplitude(Point* _p, int _top, int _perm):p(_p),on(1),top(_top),perm(_perm) {}
  };


  typedef  std::vector<Single_Vertex *>  Vertex_List;
  typedef  std::map<ATOOLS::Flavour, Vertex_List> Vertex_Table;
  typedef  std::vector<Pre_Amplitude> Pre_Ampl_List;


  class Mpi_Tag {
  public:
    enum code {  complete_buffer = 11,
		 partial_buffer  = 12,
                 empty_buffer    = 13,
                 on_switches     = 15,
		 cfc             = 22,
                 dummy           = 42};
  };

  struct Prea_Buffer {
    MPI_Point    * buff;
    int          * counters;
    int            nbuff;
    Mpi_Tag::code  tag;
#ifdef _USE_MPI_
    MPI::Request   req;
#endif
  };

  typedef  std::vector<Prea_Buffer> Prea_Buffers;
  
  class Amplitude_Generator {
    static int NMAX;
    ATOOLS::Flavour     * fl;
    int                    * b;
    Interaction_Model_Base * p_model;
    Topology               * top;
    int                      N,nEW,nQCD;

    Single_Topology        * single_top;
    int                      prenum;
    Pre_Amplitude          * prea;
    Basic_Sfuncs           * BS;
    String_Handler         * shand;
    
    
    
    Vertex_Table     v_table;
    //exhibits the points, the topology and the permutation number
    Pre_Ampl_List    prea_table;
    //  Counter_List     counter_table;
    
    void Set_End(Point*,int* &,int&);
    void Next_P(Point*,Point* &);
    void Print_P(Point*);
    int  Match_Vertex(Single_Vertex*,ATOOLS::Flavour*,Complex*);
    int  Check_End(Point*,ATOOLS::Flavour);  
    void Set_Props(Point*,int,Single_Amplitude* &,int*, int, int);
    int  Kill_Off(Single_Amplitude* &);
    int  Single_Compare(Point*,Point*);
    int  FindQEDOrder(Point*,int&);
    int  FindQCDOrder(Point*,int&);
    void KillHigherOrders(Single_Amplitude* &);
    void Compare(Single_Amplitude* &);
    void Unite(Point*,Point*);
    int  CompareColors(Color_Function*,std::vector<int>,Color_Function*,std::vector<int>);

    void CreateSingleAmplitudes(Single_Amplitude* &);
    void CheckFor4Vertices(Single_Amplitude* &);
    int  EvalPointlist(Point*&, Point*&,Point*&,std::vector<Point*>&);
    int  ShrinkProps(Point*&,Point*&,Point*&,Point*&,std::vector<Point*>&);
    Point* FindNext(Point*);
    

    MPI_Point *              s_buffer;
    Prea_Buffers             r_buffers;
    int *                    s_counter_buffer;
    void CollectPreAmplitudes();
    void DistributePreAmplitudes();
    void CommitMPITypes() ;
    void ReplaceVertex(Point *);
  
    int Count4G(Point * p);
    int CountRealAmplitudes(Single_Amplitude* first);
  public:
    static int * DivideComparisons(int cpu_size, int nampl, int rank);
    Amplitude_Generator(int,ATOOLS::Flavour *,int *,Interaction_Model_Base *,Topology *,
			int,int,Basic_Sfuncs *,String_Handler *);
    ~Amplitude_Generator(); 
    Single_Amplitude* Matching();
  };
}
#endif









