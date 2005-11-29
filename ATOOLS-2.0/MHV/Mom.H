#ifndef Mom_H
#define Mom_H

#include <iostream.h>
#include <vector>
#include "Vector.H"
#include "MyComplex.H"
#include "Pfunc.H"
#include "MHV_Calculator.H"
#include "PermStore.H"

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



  class Fullamplitude_MHV_test {
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
    Fullamplitude_MHV_test(MomentumList*,int* hl,int* pl); 

    //destructor
    //~Fullamplitude_MHV_test();


    Complex Calculate(int print);
    
  };



 
  class Fullamplitude_MHV {
  private:
    
    PermStore* permstore;
    MHV_Calculator* calc;

    int n_part; 
    Complex ampnc;
    Complex ampsq;

    int *plist, *qlist;
    int *perm, *permconj, *permgl, *permconjgl;
    std::vector<int> *permtb;
  
    int print;
    bool colorstore;
    int m_qpos, m_qposconj;

    void PermutationStoreColor_pureg(int,int**);
    void PermutationStoreAmp_pureg(int,int**);
    void Permutation_pureg(int,int**,bool);

    void PermutationStoreColor_quark2(int,std::vector<int>**);
    void PermutationStoreAmp_quark2(int,std::vector<int>**);
    void Permutation_quark2(int,int**,bool);

    //void PermutationStoreColor_quark4(int,int**);
    //void PermutationStoreAmp_quark4(int,int**);
    void Permutation_quark4(int,int**,bool);


  public:
    
    // constructor
    Fullamplitude_MHV(int n_part,int* plist,int print); 
    
    //destructor
    ~Fullamplitude_MHV();


    Complex MSquare(MomentumList*,int* hlist);
    
  };

  
}// end of namespace AMEGIC

#endif
