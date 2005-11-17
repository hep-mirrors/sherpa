#ifndef Mom_H
#define Mom_H
/*!
  \file Color.H
  \brief Declares the classes Color_Term,
  CNumber, Fundamental, Adjoint and Expression.
*/
#include <iostream.h>
#include "Node.H"
#include "Vector.H"
#include "MyComplex.H"
#include "MHV_Calculator.H"
#include "Pfunc.H"

using namespace ATOOLS;

namespace AMEGIC {

  class MHV_Calculator;

  
  class Mom: public Vec4D {
  private:
 
    int hel;
    int part;
    

  public: 

    // constructor
    Mom(Vec4D&,int,int) ; 
    Mom(Vec4D&) ; 

    // destructor
    //~Mom();

    int GetHelicity(); 
    int GetPartType();
    void Print();

  };

  class MomentumList: public std::vector<Mom*> {
  private:

    int* hlist; 
    int* plist;
    size_t m_size;   
    double etha[2];

    bool Get(const char* file);
    void MakeHlist(); 
    void MakePlist();
   
  public:

    // constructor
    MomentumList(const char*  file); 

    //destructor
    ~MomentumList();

    int* GetHList(); 
    int* GetPList();
    inline size_t Size() {return m_size;}

    Complex S0(int,int);
    Complex S1(int,int);
    int GetMomNumber(Pfunc*);
    Vec4D Momentum(size_t);

    void Print(); 

  };



  class Fullamplitude_MHV {
  private:
    
    int print;
    MomentumList* momentumlist;
    int *hlist, *plist, *perm;
    Complex amp;
    MHV_Calculator* calc;
    size_t* ma;
    void Permutation_pureg(int,int**);
    void Permutation_quark(int,int**);


  public:
    
    // constructor
    Fullamplitude_MHV(MomentumList*,int* hl,int* pl); 

    //destructor
    //~Fullamplitude_MHV();


    Complex Calculate(int print);
    
  };

}// end of namespace AMEGIC

#endif
