#ifndef AMEGIC_Amplitude_Amplitude_Generator_H
#define AMEGIC_Amplitude_Amplitude_Generator_H

#include "ATOOLS/Phys/Flavour.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Vertex.H"
#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Amplitude/Single_Amplitude.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Handler.H"

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
    MODEL::Model_Base      * p_model;
    Topology               * top;
    int                      N,nEW,nQCD;

    Single_Topology        * single_top;
    std::vector<Pre_Amplitude> prea;
    Basic_Sfuncs           * BS;
    String_Handler         * shand;

    bool                     m_create_4V;
    
    MODEL::Vertex_Table     v_table;
    //exhibits the points, the topology and the permutation number
    Pre_Ampl_List    prea_table;
    //  Counter_List     counter_table;
    
    void Set_End(Point*,int* &,int&);
    void Next_P(Point*,Point* &);
    void Print_P(Point*);
    int  MatchVertex(MODEL::Single_Vertex*,ATOOLS::Flavour*,std::vector<Complex>&);
    int  CheckEnd(Point*,ATOOLS::Flavour);  
    void SetProps(Point*,int,Single_Amplitude* &,int*, int, int);
    int  Kill_Off(Single_Amplitude* &);
    int  SingleCompare(Point*,Point*);
    int  FindQEDOrder(Point*,int&);
    int  FindQCDOrder(Point*,int&);
    void CountOrders(Single_Amplitude* &);
    bool CheckOrders(Point * p);
    void Compare(Single_Amplitude* &);
    void Unite(Point*,Point*);
    int  CompareColors(Point*,Point*);
    int  Compare5Vertex(Point*,Point*);

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
  
    void Kill5VertexArtefacts(Single_Amplitude* first);
    int  Is5VertexArtefact(Point* p, int &tcnt);

    int Count4G(Point * p);
    int CountRealAmplitudes(Single_Amplitude* first);
  public:
    static int * DivideComparisons(int cpu_size, int nampl, int rank);
    Amplitude_Generator(int,ATOOLS::Flavour *,int *,MODEL::Model_Base *,Topology *,
			int,int,Basic_Sfuncs *,String_Handler *,bool=true);
    ~Amplitude_Generator(); 
    Single_Amplitude* Matching();
    void GetOrders(int& qed,int& qcd) {qed=nEW;qcd=nQCD;}
  };
}
#endif









