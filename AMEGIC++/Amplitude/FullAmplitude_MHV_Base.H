#ifndef AMEGIC_Amplitude_FullAmplitude_MHV_Base_H
#define AMEGIC_Amplitude_FullAmplitude_MHV_Base_H

#include <vector>
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include "AMEGIC++/Amplitude/Zfunctions/MHVCalculator.H"
#include "ATOOLS/Math/PermStore.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Model_Base.H"

using namespace ATOOLS;

namespace AMEGIC {

  typedef std::map<int,Complex**> ColorStoreMap;

  class MHVCalculator;

    
  class FullAmplitude_MHV_Base {
  protected:

    MODEL::Model_Base *p_model;
    MODEL::Coupling_Map * m_cpls; 
    
    PermStore* p_permstore;	
    Permutation* p_permutation;
    MHVCalculator *p_calc;

    ColorStoreMap m_colormap;
    Complex** m_colorstore;
    Complex* m_ampstore;
    Complex* m_ampstore2;
	
    double m_cpl, p_norm; 
    bool colorflag;
        
    int n_part, maxn; 
	
    int *m_plist, *m_hlist;
    std::vector<Flavour*> m_flist;
    int *m_perm, *m_permgl;
    int m_emit,m_spect;
    bool m_dptgluon;
    double m_A,m_conv;
    Complex m_phase2;

    virtual void InitAmplitude();
    virtual bool AmpStore(MomentumList*);
    virtual bool AmpStoreDPT(MomentumList*);
    virtual double Result(Complex**); 
    virtual double ResultDPT(); 
    MODEL::Coupling_Data *p_aqcd, *p_aqed;
    double m_oqcd, m_oqed;
    
  public:
	
    // constructor
    FullAmplitude_MHV_Base(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist); 

    //destructor
    virtual ~FullAmplitude_MHV_Base(); 

    double MSquare(int* hlist,MomentumList* BS);   
    double MSquare(int,int);   
    double MSquareHel(MomentumList* BS);    
    void   CalculateAmps(int* hlist,MomentumList* BS);
    inline double ParticlesNorm() {return p_norm;}
    void SetSqMatrix(double A,Complex phase)  {m_A=A; m_phase2=phase*phase;}
  };
    




  class FullAmplitude_MHV_PureG: public FullAmplitude_MHV_Base {
  private:

    void InitAmplitude();
    void PermutationStoreColor(int,int**);
    void ColorStore();
    void ColorStoreDPT(int, int,Complex**);
    bool AmpStore(MomentumList*);   
    bool AmpStoreDPT(MomentumList*);   
	

  public:
	
    // constructor
    FullAmplitude_MHV_PureG(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist); 
    FullAmplitude_MHV_PureG(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist,int,int); 

    //destructor
    ~FullAmplitude_MHV_PureG();

  };





  class FullAmplitude_MHV_Q2: public FullAmplitude_MHV_Base {
  private:
    ATOOLS::Flavour m_efl,m_sfl;
    
    void InitAmplitude();
    void PermutationStoreColor(int,int** p_adr);
    void ColorStore();
    void ColorStoreDPT(int,int,Complex**);
    bool AmpStore(MomentumList *BS);		
    bool AmpStoreDPT(MomentumList *BS);		
    
    
  public:
    
    // constructor
    FullAmplitude_MHV_Q2(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist); 
    FullAmplitude_MHV_Q2(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist,int,int); 

    //destructor
    ~FullAmplitude_MHV_Q2();
      
  };
    



  
  class FullAmplitude_MHV_Q4: public FullAmplitude_MHV_Base {
  private:
    ATOOLS::Flavour m_efl,m_sfl;
    
    MHVCalculator* p_calc_partner;

    void InitAmplitude();	
    void ColorStore();
    void ColorStore(int,int,Complex**);
    bool AmpStore(MomentumList* BS);
    bool AmpStoreDPT(MomentumList* BS);
    double Result(Complex**);
    double ResultDPT();


  public:
    
    // constructor
    FullAmplitude_MHV_Q4(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist); 
    FullAmplitude_MHV_Q4(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist,int,int); 

    //destructor
    ~FullAmplitude_MHV_Q4();
	
  };
    




  class FullAmplitude_MHV_Q2L2: public FullAmplitude_MHV_Base {
  private:
       
    int *m_qlist, *m_llist;       

    void InitAmplitude();
    void ColorStore();
    bool AmpStore(MomentumList* BS);


    public:
	
    // constructor
    FullAmplitude_MHV_Q2L2(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist); 

    //destructor
    ~FullAmplitude_MHV_Q2L2();
  	
  };




//   FullAmplitude_MHV_Base* FullAmplitude_MHV_Handler(MODEL::Model_Base *model,int part,int* plist,MomentumList* BS); 
  FullAmplitude_MHV_Base* FullAmplitude_MHV_Handler(MODEL::Model_Base *model,MODEL::Coupling_Map *const cpls,int part,int* plist,
						    MomentumList* BS,bool &newamp,int emit=-1,int spec=-1); 
  

}// end of namespace AMEGIC

#endif
