#ifndef AMEGIC_Amplitude_Amplitude_Generator_H
#define AMEGIC_Amplitude_Amplitude_Generator_H

#include "ATOOLS/Phys/Flavour.H"

#include "AMEGIC++/Amplitude/Vertex.H"
#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Amplitude/Single_Amplitude.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Handler.H"

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


  class Amplitude_Generator {
    static int NMAX;
    ATOOLS::Flavour     * fl;
    int                    * b;
    Amegic_Model        * p_model;
    Topology               * top;
    int                      N,ntchan_min;
    std::vector<int> order;

    Single_Topology        * single_top;
    std::vector<Pre_Amplitude> prea;
    Basic_Sfuncs           * BS;
    String_Handler         * shand;

    bool                     m_create_4V;
    
    Vertex_Table     v_table;
    //exhibits the points, the topology and the permutation number
    Pre_Ampl_List    prea_table;
    
    void Set_End(Point*,int* &,int&);
    void Next_P(Point*,Point* &);
    void Print_P(Point*);
    int  MatchVertex(AMEGIC::Single_Vertex*,ATOOLS::Flavour*,std::vector<Complex>&);
    int  CheckEnd(Point*,ATOOLS::Flavour);  
    void SetProps(Point*,int,Single_Amplitude* &,int*, int, int);
    int  Kill_Off(Single_Amplitude* &);
    int  SingleCompare(Point*,Point*);
    void CountOrders(Single_Amplitude* &);
    bool CheckOrders(Point * p);
    bool CheckTChannels(Point * p);
    void Compare(Single_Amplitude* &);
    void Unite(Point*,Point*);
    int  CompareColors(Point*,Point*);
    int  Compare5Vertex(Point*,Point*);

    void CreateSingleAmplitudes(Single_Amplitude* &);
    void CheckFor4Vertices(Single_Amplitude* &);
    int  EvalPointlist(Point*&, Point*&,Point*&,std::vector<Point*>&);
    int  ShrinkProps(Point*&,Point*&,Point*&,Point*&,std::vector<Point*>&);
    Point* FindNext(Point*);
    
    void Kill5VertexArtefacts(Single_Amplitude* first);
    int  Is5VertexArtefact(Point* p, int &tcnt);

    int Count4G(Point * p);
    int CountRealAmplitudes(Single_Amplitude* first);
  public:
    Amplitude_Generator(int,ATOOLS::Flavour *,int *,Amegic_Model *,Topology *,
			std::vector<int>,int,Basic_Sfuncs *,String_Handler *,bool=true);
    ~Amplitude_Generator(); 
    Single_Amplitude* Matching();
  };
}
#endif









