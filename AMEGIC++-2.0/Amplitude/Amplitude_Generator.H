#ifndef Amplitude_Generator_h
#define Amplitude_Generator_h

#include "Flavour.H"
#include "Vertex.H"
#include "Topology.H"
#include "Single_Amplitude.H"
#include "Basic_Sfuncs.H"
#include "String_Handler.H"

namespace AMEGIC {
struct Pre_Amplitude {
  Point* p;
  int on;
};

class Amplitude_Generator {
  int N;
  APHYTOOLS::Flavour* fl;
  int* b;
  Topology* top;
  Single_Topology* single_top;
  int prenum;
  Pre_Amplitude* prea;
  Basic_Sfuncs* BS;
  String_Handler* shand;


  void Set_End(Point*,int* &,int&);
  void Next_P(Point*,Point* &);
  void Print_P(Point*);
  int  Match_Vertex(Single_Vertex*,APHYTOOLS::Flavour*,Complex*);
  int  Check_End(Point*,APHYTOOLS::Flavour,Vertex*); 
  void Set_Props(Point*,int,Single_Amplitude* &,int*);
  void Kill_Off(Single_Amplitude* &);
  int  Single_Compare(Point*,Point*);
  void Compare(Single_Amplitude* &);
  void Unite(Point*,Point*);
  void CheckFor4Vertices(Single_Amplitude* &);
  int  EvalPointlist(Point*&, Point*&,Point*&,std::vector<Point*>&);
  int  ShrinkProps(Point*&,Point*&,Point*&,Point*&,std::vector<Point*>&);
  Point* FindNext(Point*);
public:
  Amplitude_Generator(int,APHYTOOLS::Flavour*,int*,Topology*,Basic_Sfuncs*,
		      String_Handler*);
  ~Amplitude_Generator(); 
  Single_Amplitude* Matching();
};
}
#endif









