#ifndef AMEGIC_Amplitude_Zfunctions_Mom_H
#define AMEGIC_Amplitude_Zfunctions_Mom_H
#include "AMEGIC++/Amplitude/Zfunctions/MHVCalculator.H"

#ifndef Basic_Sfuncs_In_MHV
#include <vector>
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include "ATOOLS/Math/PermStore.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace ATOOLS;


namespace AMEGIC {

    class MHVCalculator;



    class Mom: public Vec4D {
    private:
	
	int hel, part;
	Complex l0[2], l1[2];
	double *etha;

	void Compute_l();
	

    public: 
	
	// constructors
	Mom(Vec4D&,int,int); 
	Mom(Vec4D&,double *et=0); 

	
	void Print();

	friend class MomentumList;


    };




    class MomentumList: public std::vector<Mom*> {
    private:
	
	int *hlist, *plist;
	size_t m_size, m_in, m_out;
	size_t m_propsize;
	double etha[2];    
	Complex **s0, **s1;

	bool Get(const char* file);
	void MakeHlist(); 
	void MakePlist();
	void BuildMomList();
	void Compute_s();
   
    public:

	// constructor
	MomentumList(const char*  file);  
	MomentumList(size_t in, size_t out); 
	
	//destructor
	~MomentumList();
	

	void PutMomenta(const ATOOLS::Vec4D*);

	Vec4D Momentum(int);
	int GetMomNumber(Pfunc*);
	int * GetHList(); 
	int * GetPList();

	void Print();
		
	inline Complex S0(int i,int j) { if (i>j) return s0[i][j]; else return -s0[j][i]; }
	inline Complex S1(int i,int j) { if (i>j) return s1[i][j]; else return -s1[j][i]; }
	inline size_t Size() {return m_size;}

    };



    class Fullamplitude_MHV_test {
    private:
	
	MomentumList* momentumlist;
	int *hlist, *plist, *perm;
	Complex amp;
	MHVCalculator* calc;
	size_t* ma;
	void Permutation_pureg(int,int**);
	void Permutation_quark(int,int**);
	

    public:
	
	// constructor
	Fullamplitude_MHV_test(MomentumList*,int* hl,int* pl); 
	
	//destructor
	//~Fullamplitude_MHV_test();
	
	
	Complex Calculate();
	
    };



 
    class Fullamplitude_MHV_old {
    private:
	
	PermStore* permstore;
	MHVCalculator* calc;
	
	int n_part; 
	Complex ampnc;
	Complex ampsq;
	double m_cpl, p_norm; 
	
        const int *qlist;
	int *plist, *p_hlist;
	int *perm, *permconj, *permgl, *permconjgl;
	std::vector<int> *permtb;
	
	bool colorstore;
	int m_qpos, m_qposconj;
	
	size_t *mact1, *mact2;
	size_t qindex[2];
	
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
	Fullamplitude_MHV_old(int n_part,int* plist); 
	
	//destructor
	~Fullamplitude_MHV_old();
	
	
	double MSquare(MomentumList*,int* hlist);   
	Complex MSquareHel(MomentumList*);

	inline double ParticlesNorm() {return p_norm;}
	
    };
    
 
  
}// end of namespace AMEGIC

#endif
#endif
