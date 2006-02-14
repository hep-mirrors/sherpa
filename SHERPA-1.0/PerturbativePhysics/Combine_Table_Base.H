#ifndef Combine_Table_Base_H
#define Combine_Table_Base_H

#include <vector>
#include <map>
#include "Point.H" 
#include "Vector.H"

namespace ATOOLS { class Jet_Finder; }

namespace SHERPA {

  class Leg {
  private:
    AMEGIC::Point *p_point;
    int m_anti, m_nqcd, m_nqed, m_pqcd, m_pqed;
    AMEGIC::lf::code m_type;
    friend std::ostream &operator<<(std::ostream &ostr,const Leg &leg);
  public:
    // constructors
    inline Leg(AMEGIC::Point *const point=NULL,const int anti=1):  
      p_point(point), m_anti(anti) {}
    inline Leg(const Leg &leg): 
      p_point(leg.p_point), m_anti(leg.m_anti) {}
    // member functions
    void DetermineCouplings(const int type);
    // inline functions
    inline void SetAnti(int a)                       { m_anti=a;      }
    inline void SetPoint(AMEGIC::Point *const point) { p_point=point; }
    inline int                  Anti() const  { return m_anti;  }
    inline const AMEGIC::Point *Point() const { return p_point; }
    inline ATOOLS::Flavour Flav() 
    { return m_anti==-1?p_point->fl.Bar():p_point->fl; }
    inline const int OrderQCD() const { return m_nqcd; }
    inline const int OrderQED() const { return m_nqed; }
    inline const int NQCD() const { return m_pqcd; }
    inline const int NQED() const { return m_pqed; }
    inline const AMEGIC::lf::code Type() const { return m_type; }
  };// end of class Leg

  std::ostream &operator<<(std::ostream &ostr,const Leg &leg);

  class Combine_Data;
  class Combine_Table_Base;

  class Combine_Key {
  public:
    int m_i, m_j;
    ATOOLS::Flavour m_flav;
  public:
    // constructors
    Combine_Key();
    Combine_Key(const int i,const int j, 
		const ATOOLS::Flavour &flav=ATOOLS::Flavour(ATOOLS::kf::none));
  };// end of class Combine_Key

  bool operator<(const Combine_Key &a,const Combine_Key &b);

  std::ostream& operator<<(std::ostream &ostr,const Combine_Key &ck);

  class Combine_Data {
  public:
    double m_pt2ij, m_sij, m_prop, m_coupling, m_weight;
    int    m_strong;
    Combine_Table_Base *p_down;
    std::vector<int>    m_graphs;
    friend std::ostream &operator<< (std::ostream&,const Combine_Table_Base &);
  public:
    // constructors
    Combine_Data();
    Combine_Data(const Combine_Data &cd);
    Combine_Data(const double pt2ij,const int ngraph);
    // destructor
    ~Combine_Data();
  };

  std::ostream &operator<<(std::ostream &ostr,const Combine_Data &cd);

  typedef std::map<Combine_Key,Combine_Data> CD_List;

  class Combine_Table_Base {
  protected:
    static int s_all;
    int        m_no, m_nstrong;
    int        m_isr1on, m_isr2on, m_isrshoweron;
    int        m_nlegs, m_nampl, m_graph_winner;
    double     m_x1, m_x2;
    double     m_kt2min, m_kt2QCD, m_kt2QED;

    CD_List               m_combinations;
    CD_List::iterator     m_cdata_winner;

    Combine_Table_Base  * p_up;
    Leg                ** p_legs;
    ATOOLS::Jet_Finder  * p_jf;
    ATOOLS::Vec4D       * p_moms;

    Leg **p_hard;

    Leg  CombinedLeg(Leg *legs,const int i,const int j);
    Leg *CombineLegs(Leg *legs,const int i,const int j,
		     const int nlegs,const double pt2ij);
    void CombineMoms(ATOOLS::Vec4D *moms,const int i,const int j,
		     const int maxl);
    void CombineMoms(ATOOLS::Vec4D *moms,const int i,const int j,
		     const int maxl,ATOOLS::Vec4D *&omoms);
    double      Sprime() const;
    friend std::ostream& operator<< (std::ostream&,const Combine_Table_Base &);
    friend std::ostream &operator<<(std::ostream&,const Combine_Data &);
  public:
    // constructors
    Combine_Table_Base(ATOOLS::Jet_Finder * jf, ATOOLS::Vec4D * moms, 
		       Combine_Table_Base * up, int isrmode, int isrshoweron);
    // destructor
    virtual ~Combine_Table_Base();
    // member functions
    virtual CD_List::iterator    CalcPropagator(CD_List::iterator &cit) = 0;
    virtual Combine_Table_Base * CalcJet(int nl,const double x1,const double x2,
					 ATOOLS::Vec4D *moms=NULL) = 0;
    virtual void                 FillTable(Leg **legs,const int nlegs,const int nampl) = 0;
    virtual double               GetWinner(int &i, int &j) = 0;

    bool Combinable(const Leg &a,const Leg &b) const ;

    void   IdentifyHardProcess();
    int    AddCouplings(int &nqed,int &nqcd) const;
    double MinKt2() const;
    double MinKt2QCD() const;
    double MinKt2QED() const;
    double Kt2() const;
    double Kt2QCD() const;
    double Kt2QED() const;
    // inline functions
    inline const int                  NLegs()                              const { return m_nlegs; }
    inline const int                  OrderStrong()                        const { return m_nstrong; }
    inline const int                  NAmplitudes()                        const { return m_nampl; }
    inline Combine_Table_Base * const Up()                                 const { return p_up; }
    inline ATOOLS::Vec4D              Momentum(const int i)                const { return p_moms[i]; }
    inline ATOOLS::Flavour            Flav(const int i)                    const { return p_legs[m_graph_winner][i].Flav(); }
    inline const Leg                & GetLeg(const int ngraph,const int l) const { return p_legs[ngraph][l]; }
    inline const Leg                & GetLeg(const int l)                  const { return p_legs[m_graph_winner][l]; }
    inline ATOOLS::Vec4D      * const Momenta()                            const { return p_moms; }
    inline Leg               ** const GetLegs()                            const { return p_legs; }
    inline Leg               ** const GetHardLegs()                        const { return p_hard; }
    inline void                       GetX1X2(double & x1, double & x2)    const { x1=m_x1; x2=m_x2; }
    inline static void                ResetCounter()                             { s_all=0; }
  };// end of class Combine_Table_Base

}// end of namespace SHERPA

#endif



