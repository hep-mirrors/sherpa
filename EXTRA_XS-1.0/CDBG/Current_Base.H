#ifndef Current_Base_H
#define Current_Base_H

#include "Vector.H"
#include "Flavour.H"
#include "Getter_Function.H"
#include "Node.H"

#include <vector>

namespace EXTRAXS {

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;

  typedef std::vector<std::string>   String_Vector;
  typedef std::vector<String_Vector> String_Matrix;

  typedef ATOOLS::Node<std::string> Graph_Node;

  class Vertex;

  typedef std::vector<Vertex*> Vertex_Vector;

  char ParticleType(const ATOOLS::Flavour &fl);

  size_t     MakeId(const Int_Vector &ids);
  Int_Vector MakeId(const size_t &id,const size_t &n);

  class Current_Base {
  protected:

    ATOOLS::Flavour m_fl;

    Vertex_Vector m_in, m_out;
    Int_Vector    m_id, m_fid;

    ATOOLS::Vec4D m_p;

    size_t m_key;
    char   m_type;

    virtual std::string CLabel() const;

    void CollectGraphs(Graph_Node *graph,
		       const std::string &lastv) const;

  public:

    Current_Base(const ATOOLS::Flavour &fl);

    virtual ~Current_Base();

    // member functions
    void FindPermutations();

    template <typename CType> void AddJ(const CType &j);

    template <typename CType> const std::vector<CType> &J() const;
    
    virtual void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			    const int cr,const int ca) = 0;
    virtual void SetGauge(const ATOOLS::Vec4D &k) = 0;

    virtual void Evaluate() = 0;
    virtual void AddPropagator() = 0;

    virtual Complex Contract(const Current_Base &c,
			     const long unsigned int &hm,
			     const long unsigned int &hp) const = 0;

    virtual void Print() const = 0;
    virtual char Type() const = 0;

    size_t CId() const;

    void CollectGraphs(Graph_Node *graph) const;

    // inline functions
    inline void SetP(const ATOOLS::Vec4D &p) { m_p=p; }

    inline void SetId(const Int_Vector &id)   { m_id=id;   }
    inline void SetFId(const Int_Vector &fid) { m_fid=fid; }
    inline void SetKey(const size_t &key)     { m_key=key; }

    inline ATOOLS::Flavour Flav() const { return m_fl; }

    inline ATOOLS::Vec4D P() const { return m_p; }

    inline const Int_Vector &Id() const   { return m_id;   }
    inline const Int_Vector &FId() const  { return m_fid;  }
    inline size_t            Key() const  { return m_key;  }

    inline void AttachIn(Vertex *const v)  { m_in.push_back(v);  }
    inline void AttachOut(Vertex *const v) { m_out.push_back(v); }

  };// end of class Current_Base
  
  struct Current_Key {
    ATOOLS::Flavour m_fl;
    std::string     m_model;
    inline Current_Key(const ATOOLS::Flavour &fl,const std::string &model): 
      m_fl(fl), m_model(model) {}
    std::string Type() const;
  };// end of struct Current_Key

  typedef ATOOLS::Getter_Function<Current_Base,Current_Key> Current_Getter;

  std::ostream &operator<<(std::ostream &str,const Current_Base &c);

}// end of namespace EXTRAXS

#endif
