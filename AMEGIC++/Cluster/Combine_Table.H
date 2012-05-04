#ifndef AMEGIC_Cluster_Combine_Table_H
#define AMEGIC_Cluster_Combine_Table_H

#include <vector>
#include <map>
#include "AMEGIC++/Main/Point.H" 
#include "AMEGIC++/Main/Process_Base.H" 
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PDF/Main/Cluster_Definitions_Base.H"

namespace PDF { class Cluster_Definitions_Base; }

namespace AMEGIC {

  std::ostream &operator<<(std::ostream &str,const std::vector<int> &info);

  class Leg {
  private:
    AMEGIC::Point *p_point;
    int m_anti, m_nqcd, m_nqed, m_pqcd, m_pqed, m_ext, m_qcdjets, m_kin;
    long unsigned int m_id;
    double *p_qmin;
    std::string m_type;
    ATOOLS::Flavour m_mapfl;
    friend std::ostream &operator<<(std::ostream &ostr,const Leg &leg);
  public:
    // constructors
    Leg(AMEGIC::Point *const point=NULL,const int anti=1);
    Leg(const Leg &leg);
    // member functions
    void DetermineCouplings(const int type);
    // inline functions
    inline void SetAnti(int a)                       { m_anti=a;      }
    inline void SetExternal(int ext)                 { m_ext=ext;     }
    inline void SetPoint(AMEGIC::Point *const point) { p_point=point; }
    inline int                  Anti() const      { return m_anti;  }
    inline int                  External() const  { return m_ext;   }
    inline const AMEGIC::Point *Point() const     { return p_point; }
    inline void SetQCDJets(const int qcdjets) { m_qcdjets=qcdjets; }
    inline int  QCDJets() const               { return m_qcdjets;  }
    inline void   SetID(const size_t &id) { m_id=id;     }
    inline size_t ID() const              { return m_id; }
    inline void   SetKin(const int kin) { m_kin=kin; }
    inline int    Kin() const       { return m_kin;       }
    inline void    SetPQMin(double *const qmin) { p_qmin=qmin;   }
    inline double *PQMin() const                { return p_qmin; }
    inline void   SetQMin(const double qmin) { *p_qmin=qmin; }
    inline double QMin() const { return p_qmin!=NULL?*p_qmin:-1.0; }
    inline ATOOLS::Flavour Flav() 
    { return m_anti==-1?m_mapfl.Bar():m_mapfl; }
    inline int OrderQCD() const { return m_nqcd; }
    inline int OrderQED() const { return m_nqed; }
    inline int NQCD() const { return m_pqcd; }
    inline int NQED() const { return m_pqed; }
    inline const std::string &Type() const { return m_type; }
    inline void SetMapFlavour(const ATOOLS::Flavour &fl) { m_mapfl=fl; }
    inline ATOOLS::Flavour MapFlavour() const { return m_mapfl; }
  };// end of class Leg

  std::ostream &operator<<(std::ostream &ostr,const Leg &leg);

  class Combine_Data;
  class Combine_Table;

  class Combine_Key {
  public:
    int m_i, m_j, m_k;
    ATOOLS::Flavour m_flav;
  public:
    // constructors
    Combine_Key();
    Combine_Key(const int i,const int j,const int k,
		const ATOOLS::Flavour &flav);
  };// end of class Combine_Key

  bool operator<(const Combine_Key &a,const Combine_Key &b);

  std::ostream& operator<<(std::ostream &ostr,const Combine_Key &ck);

  class Combine_Data {
  public:
    PDF::CParam m_pt2ij;
    int    m_strong, m_calc;
    Combine_Table *p_down;
    ATOOLS::Flavour m_mo;
    std::vector<int>    m_graphs;
    friend std::ostream &operator<< (std::ostream&,const Combine_Table &);
  public:
    // constructors
    Combine_Data();
    Combine_Data(const double pt2ij,const int ngraph);
    // destructor
    ~Combine_Data();
  };

  std::ostream &operator<<(std::ostream &ostr,const Combine_Data &cd);

  typedef std::map<Combine_Key,Combine_Data> CD_List;

  class Combine_Table {
  protected:
    static int s_all;
  
    ATOOLS::Mass_Selector *p_ms;

    int        m_no, m_nstrong;
    int        m_nlegs, m_nampl, m_graph_winner;

    CD_List               m_combinations, m_rejected;
    CD_List::iterator     m_cdata_winner;

    Combine_Table  * p_up;
    Leg                ** p_legs;
    PDF::Cluster_Definitions_Base *p_clus;
    ATOOLS::Vec4D       * p_moms;

    Leg **p_hard;
    int **p_hardc, *p_channel;
    double *p_scale;

    ATOOLS::Vec4D *p_save_moms;
    size_t         m_nl;

    AMEGIC::Process_Base *p_proc;
    
    ATOOLS::DecayInfo_Vector m_decids, *p_decids;

    Leg  CombinedLeg(Leg *legs,const int i,const int j);
    Leg *CombineLegs(Leg *legs,const int i,const int j,
		     const int nlegs,const int kin);
    bool CombineMoms(ATOOLS::Vec4D *moms,const int i,const int j,
		     const int maxl);
    bool CombineMoms(ATOOLS::Vec4D *moms,const int i,const int j,
		     const int maxl,ATOOLS::Vec4D *&omoms);
    void AddPossibility(const int i,const int j,const int k,
			const int ngraph);
    bool SelectWinner(const size_t &mode);
    Combine_Table *CreateNext();
    Combine_Table *NextTable(Combine_Table *tab,
			     const int mode,const double &kt2);

    double Sprime() const;

    ATOOLS::Flavour MatchFlavour(const Leg &a,const Leg &b,const Leg &c,int mode) const;
    friend std::ostream& operator<< (std::ostream&,const Combine_Table &);
    friend std::ostream &operator<<(std::ostream&,const Combine_Data &);
  public:
    // constructors
    Combine_Table(AMEGIC::Process_Base *const proc,
		  ATOOLS::Mass_Selector *const ms,
		  PDF::Cluster_Definitions_Base *clus, 
		  ATOOLS::Vec4D * moms,Combine_Table * up,
		  ATOOLS::DecayInfo_Vector *const decids);
    // destructor
    ~Combine_Table();
    // member functions
    CD_List::iterator    CalcPropagator(CD_List::iterator &cit,int mode);
    Combine_Table * CalcJet(int nl,ATOOLS::Vec4D *moms,const size_t mode,const double &kt2);
    void                 FillTable(Leg **legs,const int nlegs,const int nampl);
    double               GetWinner(int &i, int &j,int &k,double &mu2);
    bool                 IdentifyHardProcess();
    int                  IdentifyHardPropagator(double &scale) const;

    bool Combinable(const Leg &a,const Leg &b,const int i,const int j) const ;

    int    AddCouplings(int &nqed,int &nqcd) const;
    // inline functions
    inline int                  Number()                             const { return m_no; }
    inline int                  NLegs()                              const { return m_nlegs; }
    inline int                  OrderStrong()                        const { return m_nstrong; }
    inline int                  NAmplitudes()                        const { return m_nampl; }
    inline Combine_Table * Up()                                 const { return p_up; }
    inline Combine_Table * Down()                               const { return m_combinations.empty()?NULL:m_cdata_winner->second.p_down; }
    inline ATOOLS::Vec4D              Momentum(const int i)                const { return p_moms[i]; }
    inline ATOOLS::Flavour            Flav(const int i)                    const { return p_legs[m_graph_winner][i].Flav(); }
    inline const Leg                & GetLeg(const int ngraph,const int l) const { return p_legs[ngraph][l]; }
    inline Leg                & GetLeg(const int l)                  const { return p_legs[m_graph_winner][l]; }
    inline ATOOLS::Vec4D      * Momenta()                            const { return p_moms; }
    inline Leg               ** GetLegs()                            const { return p_legs; }
    inline Leg               ** GetHardLegs()                        const { return p_hard; }
    inline int               ** GetHardInfo()                        const { return p_hardc; }
    inline static void                ResetCounter()                             { s_all=0; }
    inline const ATOOLS::DecayInfo_Vector &Decays() const { return m_decids; }

  };// end of class Combine_Table

}// end of namespace AMEGIC

#endif



