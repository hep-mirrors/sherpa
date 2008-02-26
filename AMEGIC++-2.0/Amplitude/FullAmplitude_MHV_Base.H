#ifndef FullAmplitude_MHV_Base_H
#define FullAmplitude_MHV_Base_H

#include <vector>
#include "Vector.H"
#include "MyComplex.H"
#include "Pfunc.H"
#include "MHVCalculator.H"
#include "PermStore.H"
#include "Permutation.H"
#include "Flavour.H"
#include "Model_Base.H"

using namespace ATOOLS;

namespace AMEGIC {

  class MHVCalculator;

    
  class FullAmplitude_MHV_Base {
  protected:

    MODEL::Interaction_Model_Base *p_model;

    PermStore* p_permstore;	
    Permutation* p_permutation;
    MHVCalculator *p_calc;
    MomentumList* p_BS;

    Complex** m_colorstore;
    Complex* m_ampstore;
	
    double m_cpl, p_norm; 
    bool colorflag;
        
    int n_part, maxn; 
	
    int *m_plist, *m_hlist;
    std::vector<Flavour*> m_flist;
    int *m_perm, *m_permgl;

    virtual void InitAmplitude();
    virtual bool AmpStore();
    virtual double Result(); 


  public:
	
    // constructor
    FullAmplitude_MHV_Base(MODEL::Model_Base *model,int part,int* plist,MomentumList* BS); 

    //destructor
    virtual ~FullAmplitude_MHV_Base(); 

    virtual double MSquare(int* hlist);   
    double MSquareHel();    
    inline double ParticlesNorm() {return p_norm;}
  };
    




  class FullAmplitude_MHV_PureG: public FullAmplitude_MHV_Base {
  private:

    void InitAmplitude();
    void PermutationStoreColor(int,int**);
    void ColorStore();
    bool AmpStore();   
	

  public:
	
    // constructor
    FullAmplitude_MHV_PureG(MODEL::Model_Base *model,int part,int* plist,MomentumList *BS); 

    //destructor
    ~FullAmplitude_MHV_PureG();

  };





  class FullAmplitude_MHV_Q2: public FullAmplitude_MHV_Base {
  private:
    
    void InitAmplitude();
    void PermutationStoreColor(int,int** p_adr);
    void ColorStore();
    bool AmpStore();		
    
    
  public:
    
    // constructor
    FullAmplitude_MHV_Q2(MODEL::Model_Base *model,int part,int* plist,MomentumList* BS); 

    //destructor
    ~FullAmplitude_MHV_Q2();
      
  };
    



  
  class FullAmplitude_MHV_Q4: public FullAmplitude_MHV_Base {
  private:
    
    MHVCalculator* p_calc_partner;

    void InitAmplitude();	
    void ColorStore();
    bool AmpStore();
    double Result();


  public:
    
    // constructor
    FullAmplitude_MHV_Q4(MODEL::Model_Base *model,int part,int* plist,MomentumList* BS); 

    //destructor
    ~FullAmplitude_MHV_Q4();
	
  };
    




  class FullAmplitude_MHV_Q2L2: public FullAmplitude_MHV_Base {
  private:
       
    int *m_qlist, *m_llist;       

    void InitAmplitude();
    void ColorStore();
    bool AmpStore();


    public:
	
    // constructor
    FullAmplitude_MHV_Q2L2(MODEL::Model_Base *model,int part,int* plist,MomentumList* BS); 

    //destructor
    ~FullAmplitude_MHV_Q2L2();
  	
  };




  FullAmplitude_MHV_Base* FullAmplitude_MHV_Handler(MODEL::Model_Base *model,int part,int* plist,MomentumList* BS); 
  

}// end of namespace AMEGIC

#endif
