#ifndef PHASIC_Main_Color_Integrator_H
#define PHASIC_Main_Color_Integrator_H

#include <vector>
#include <iostream>
#include <map>

namespace ATOOLS { class Cluster_Amplitude; }

namespace PHASIC {

  struct cls {

    enum scheme {
      unknown   = 0,
      sum       = 1,
      sample    = 2
    };

    static std::map<std::string,std::string> ColorSchemeTags();

  };// end of struct cls

  std::ostream &operator<<(std::ostream &str,const cls::scheme &s);

  class Representation {
  protected:

    size_t m_id, m_i, m_j;
    int    m_type, m_act;

  public:

    // constructor
    inline Representation(const size_t &id,
			  const int &type,const int &act): 
      m_id(id), m_i(0), m_j(0), m_type(type), m_act(act) {}

    // inline functions
    inline void SetI(const size_t &i) { m_i=i; }
    inline void SetJ(const size_t &j) { m_j=j; }

    inline size_t I() const { return m_i; }
    inline size_t J() const { return m_j; }

    inline size_t Id() const { return m_id; }

    inline int Type() const { return m_type; }
    inline int Act() const  { return m_act;  }

  };// end of class Representation

  std::ostream &operator<<(std::ostream &ostr,const Representation &v);
  
  typedef std::vector<int>    Int_Vector;
  typedef std::vector<double> Double_Vector;

  typedef short unsigned int      Idx_Type;
  typedef std::vector<Idx_Type>   Idx_Vector;
  typedef std::vector<Idx_Vector> Idx_Matrix;

  class Color_Integrator {
  private:

    std::vector<Representation*> m_ids;

    Idx_Matrix m_orders;
    Idx_Vector m_lastids, m_lastperm, m_nexti;

    Double_Vector m_weights, m_cweights, m_alpha;

    size_t m_pairs, m_lastconf, m_alphamode, m_csmode;
    bool   m_check, m_on, m_otfcc, m_fincc, m_valid, m_nogen;

    double m_weight, m_n, m_nv, m_max, m_over, m_cmax, m_mean, m_cweight;

    double Factorial(const double &n) const;

    size_t GenerateIndex();

    bool DiceColours();
    int  Dice();
    bool LookUp();
    
    int  ConstructConfigurations(Idx_Vector ids,Idx_Vector perm,
				 bool sing,double weight,
				 Idx_Vector &nexti,bool one,size_t depth);
    void InitConstruction(Idx_Vector &ids,Idx_Vector &perm,
			  Idx_Vector &nexti);
    int  ConstructConfigurations();

    bool TrivialCheck();
    bool CheckPermutation(const Idx_Vector &perm);

    bool AddConfiguration(const size_t &i);
    bool RInitialize();

  public:

    // constructor
    Color_Integrator();

    // destructor
    ~Color_Integrator();

    // member functions
    bool ConstructRepresentations(const Idx_Vector &ids,
				  const Int_Vector &types,
				  const Int_Vector &acts);

    bool NextOrder();
    bool GenerateOrders();
    bool GeneratePoint();
    bool GenerateType(const size_t &type,const bool orders=false);
    bool Initialize();

    size_t IdentifyType(const Idx_Vector &perm) const;

    void SetI(const Int_Vector &i);
    void SetJ(const Int_Vector &j);

    Int_Vector I() const;
    Int_Vector J() const;

    bool SetConfiguration(const size_t &id);

    void SetAlpha(const Double_Vector &alpha);

    void SetPoint(const ATOOLS::Cluster_Amplitude *const ampl);

    // inline functions
    inline void SetNoGen(const bool nogen)     { m_nogen=nogen;     }
    inline void SetCheck(const bool check)     { m_check=check;     }
    inline void SetOTFCC(const bool otf)       { m_otfcc=otf;       }

    inline void SetOn(const bool on)           { m_on=on;           }
    inline void SetAlphaMode(const size_t am)  { m_alphamode=am;    }

    inline bool NoGen() const   { return m_nogen;   }
    inline bool Check() const   { return m_check;   }
    inline bool OTFCC() const   { return m_otfcc;   }

    inline bool On() const { return m_on; }

    inline bool ValidPoint() const   { return m_valid; }

    inline void RestartOrders() { m_fincc=true; }

    inline double GlobalWeight() const { return m_cweight; }

    inline const Idx_Matrix    &Orders() const  { return m_orders; }
    inline const Double_Vector &Weights() const { return m_weights; }

  };// end of class Color_Integrator

}// end of namespace PHASIC

#endif
