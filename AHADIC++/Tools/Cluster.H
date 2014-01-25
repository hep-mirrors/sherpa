#ifndef AHADIC_Tools_Cluster_H
#define AHADIC_Tools_Cluster_H

#include "AHADIC++/Tools/Proto_Particle.H"

namespace AHADIC {
  static ATOOLS::Vec4D s_axis = ATOOLS::Vec4D(1.,1.,0.,0.);

  class Cluster;
  class Cluster_List : public std::list<Cluster * > {
  private:
    static std::list<Cluster_List * > s_actives;
  public:
    Cluster_List() {
      s_actives.push_back(this);
    }
    
    ~Cluster_List() {
      s_actives.remove(this);
    }
    
    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Cluster_Lists -----"<<std::endl;
      for (std::list<Cluster_List * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->size()<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  typedef Cluster_List::iterator       Cluster_Iterator;
  typedef Cluster_List::const_iterator Cluster_Const_Iterator;

  std::ostream & operator<<(std::ostream & s, const Cluster_List & pl);


  class Cluster {
  private:
    bool                m_active;
    Proto_Particle    * p_trip, * p_anti;
    ATOOLS::Vec4D       m_momentum;
    ATOOLS::Flavour     m_flav;
    double              m_kt, m_z;

    ATOOLS::Poincare    m_boost,m_rotate;
    bool                m_hasboost,m_hasrotate;

    Cluster_List        m_clusters;
    Cluster           * p_left, * p_right;
    Cluster           * p_prev;
    Cluster           * p_nbtrip, * p_nbanti;
    long int            m_number;

    std::vector<ATOOLS::Flavour> m_decayproducts;

    static long int     s_cluster_count, s_cluster_number;
    static std::list<Cluster *> s_actives;

  public:
    Cluster(ATOOLS::Vec4D=ATOOLS::Vec4D(0.,0.,0.,0.),
	    ATOOLS::Flavour=ATOOLS::Flavour(kf_cluster),
	    bool active=true);
    Cluster(Proto_Particle *,Proto_Particle *);
    ~Cluster();

    void               Update();
    bool               CheckConsistency(std::ostream & s,std::string method);
    bool               EnsureMomentum();
    ATOOLS::Blob     * ConstructDecayBlob();

    ATOOLS::Particle * GetSelf() const;

    void               RescaleMomentum(ATOOLS::Vec4D);
    void               BoostInCMSAndRotateOnZ();
    void               RotateAndBoostBack();
    void               BoostInCMS();
    void               BoostBack();
    void               Boost(ATOOLS::Poincare &);
    void               BoostBack(ATOOLS::Poincare &);
    void               Rotate(ATOOLS::Poincare &);
    void               RotateBack(ATOOLS::Poincare &);
    void               BoostBack(ATOOLS::Vec4D &);
    void               RotateAndBoostBack(ATOOLS::Vec4D &);

    bool               Active() const;
    ATOOLS::Flavour    GetFlav() const;
    Cluster_List *     GetClusters();
    Cluster *          GetLeft()  const;
    Cluster *          GetRight() const;
    Cluster *          GetPrev()  const;
    Cluster *          GetNBTrip()  const;
    Cluster *          GetNBAnti() const;
    void               SetActive(bool);
    void               SetFlav(ATOOLS::Flavour);
    void               SetLeft(Cluster *);
    void               SetRight(Cluster *);
    void               SetPrev(Cluster *);
    void               SetNBTrip(Cluster *);
    void               SetNBAnti(Cluster *);


    Proto_Particle *   GetTrip() const;
    Proto_Particle *   GetAnti() const;
    void               SetTrip(Proto_Particle *);
    void               SetAnti(Proto_Particle *);
    ATOOLS::Vec4D      Momentum() const;
    void               SetMomentum(const ATOOLS::Vec4D);
    double             Mass2() const;
    double             Mass() const;
    double             KT() const;
    double             Z() const;
    void               SetKT(const double);
    void               SetZ(const double);

    void               push_back(ATOOLS::Flavour);
    void               push_back(Cluster * cluster);
    void               clear();
    size_t             size() const;
    ATOOLS::Flavour    operator[](size_t pos);


    void               Print();
    void               Delete();

    friend std::ostream& operator<<(std::ostream &,const Cluster &);

    inline long int    Number()             { return m_number; }
    static long int    RemainingClusters()  { return s_cluster_count; }
    static void        ResetClusterNumber() { s_cluster_number = 0; }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Clusters -----"<<std::endl;
      for (std::list<Cluster * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->m_number<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };


  
  inline bool Cluster::Active() const                  { return m_active;  }
  inline Cluster * Cluster::GetLeft() const            { return p_left; }
  inline Cluster * Cluster::GetRight() const           { return p_right;}
  inline Cluster * Cluster::GetPrev() const            { return p_prev; }
  inline Cluster * Cluster::GetNBTrip() const          { return p_nbtrip;  }
  inline Cluster * Cluster::GetNBAnti() const          { return p_nbanti;  }
  inline ATOOLS::Flavour Cluster::GetFlav() const      { return m_flav; }
  inline Cluster_List * Cluster::GetClusters()         { return &m_clusters; } 
  inline Proto_Particle * Cluster::GetTrip() const     { return p_trip; }
  inline Proto_Particle * Cluster::GetAnti() const     { return p_anti; }
  inline ATOOLS::Vec4D Cluster::Momentum() const       { return m_momentum; }
  inline double Cluster::Mass2() const { return m_momentum.Abs2(); }
  inline double Cluster::Mass() const  { return sqrt(m_momentum.Abs2()); }
  inline void Cluster::SetActive(bool act)             { m_active = act;}
  inline void Cluster::SetFlav(ATOOLS::Flavour fl)     { m_flav= fl; }
  inline void Cluster::SetLeft(Cluster * left)         { p_left= left;  }
  inline void Cluster::SetRight(Cluster * right)       { p_right  = right; }
  inline void Cluster::SetPrev(Cluster * prev)         { p_prev= &*prev;}
  inline void Cluster::SetNBTrip(Cluster * clu)        { p_nbtrip = clu;}
  inline void Cluster::SetNBAnti(Cluster * clu)        { p_nbanti = clu;}
  inline void Cluster::SetTrip(Proto_Particle * trip) { p_trip = trip; }
  inline void Cluster::SetAnti(Proto_Particle * anti) { p_anti = anti; }
  inline void Cluster::SetMomentum(const ATOOLS::Vec4D p) { m_momentum = p; }
  inline void Cluster::push_back(ATOOLS::Flavour flav) 
  { m_decayproducts.push_back(flav); }
  inline void  Cluster::push_back(Cluster * cluster) 
  { m_clusters.push_back(cluster); }
  inline void   Cluster::clear()   { m_decayproducts.clear(); }
  inline size_t Cluster::size() const { return m_decayproducts.size(); }
  inline ATOOLS::Flavour Cluster::
  operator[](size_t pos)  { return m_decayproducts[pos]; }
  inline double Cluster::KT() const{ return m_kt; }
  inline double Cluster::Z() const { return m_z; }
  inline void  Cluster::SetKT(const double kt)  { m_kt = kt; }
  inline void  Cluster::SetZ(const double z) { m_z = z; }
}
#endif
