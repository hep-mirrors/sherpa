#ifndef COMIX_Amplitude_Current_Base_H
#define COMIX_Amplitude_Current_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "COMIX/Main/Model.H"
#include "ATOOLS/Org/Node.H"

#include <vector>

namespace COMIX {

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;

  typedef std::map<long unsigned int,size_t>       LongUInt_Map;
  typedef std::map<long unsigned int,LongUInt_Map> LongUIntMap_Matrix;

  typedef std::vector<std::string>   String_Vector;
  typedef std::vector<String_Vector> String_Matrix;

  typedef ATOOLS::Node<std::string> Graph_Node;

  class Vertex_Base;

  typedef std::vector<Vertex_Base*> Vertex_Vector;

  char ParticleType(const ATOOLS::Flavour &fl);

  struct Current_Key {
    ATOOLS::Flavour m_fl;
    Model *p_model;
    inline Current_Key(const ATOOLS::Flavour &fl,
		       Model *const model): 
      m_fl(fl), p_model(model) {}
    std::string Type() const;
  };// end of struct Current_Key

  class Current_Base {
  protected:

    ATOOLS::Flavour m_fl;

    Vertex_Vector m_in, m_out;
    Int_Vector    m_id, m_fid;

    ATOOLS::Vec4D m_p;

    size_t m_key, m_oew, m_oqcd;
    char   m_type;

    double m_mass, m_width;
    bool   m_msv, m_zero;
    int    m_dir, m_cut, m_osd;

    mutable std::string m_psinfo;

    virtual std::string CLabel() const;

    void CollectGraphs(Graph_Node *graph,
		       const std::string &lastv) const;

  public:

    Current_Base(const Current_Key &key);

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

    template <class SType> inline std::complex<SType> 
    Contract(const Current_Base &c,const long unsigned int &hm,
	     const long unsigned int &hp) const;
    template <class SType> inline void 
    Contract(const Current_Base &c,const LongUIntMap_Matrix &hmp,
	     std::vector<std::complex<SType> > &ress) const;

    virtual void Print() const = 0;
    virtual char Type() const = 0;

    size_t CId() const;

    void CollectGraphs(Graph_Node *graph) const;

    void DetachOut(Vertex_Base *const v);

    std::string PSInfo() const;

    virtual void ResetZero() = 0;

    // inline functions
    inline void SetP(const ATOOLS::Vec4D &p) { m_p=p; }

    inline void SetId(const Int_Vector &id)   { m_id=id;   }
    inline void SetFId(const Int_Vector &fid) { m_fid=fid; }
    inline void SetKey(const size_t &key)     { m_key=key; }

    inline void SetOrderEW(const size_t &oew)   { m_oew=oew;   }
    inline void SetOrderQCD(const size_t &oqcd) { m_oqcd=oqcd; }

    inline void SetDirection(const int &dir) { m_dir=dir; }
    inline void SetCut(const int &cut)       { m_cut=cut; }
    inline void SetOnShell(const int &osd)   { m_osd=osd; }

    inline ATOOLS::Flavour Flav() const { return m_fl; }

    inline ATOOLS::Flavour RFlav() const 
    { return m_fl.IsAnti()?m_fl.Bar():m_fl; }
    inline ATOOLS::Flavour AFlav() const 
    { return m_fl.IsAnti()?m_fl:m_fl.Bar(); }

    inline ATOOLS::Vec4D P() const { return m_p; }

    inline const Int_Vector &Id() const   { return m_id;   }
    inline const Int_Vector &FId() const  { return m_fid;  }
    inline size_t            Key() const  { return m_key;  }

    inline size_t OrderEW() const  { return m_oew;  }
    inline size_t OrderQCD() const { return m_oqcd; }

    inline int Direction() const { return m_dir; }
    inline int Cut() const       { return m_cut; }
    inline int OnShell() const   { return m_osd; }

    inline void AttachIn(Vertex_Base *const v)  { m_in.push_back(v);  }
    inline void AttachOut(Vertex_Base *const v) { m_out.push_back(v); }

    inline size_t NIn() const  { return m_in.size();  }
    inline size_t NOut() const { return m_out.size(); }

    inline const Vertex_Vector &In() const  { return m_in;  }
    inline const Vertex_Vector &Out() const { return m_out; }

    inline bool Dangling() const { return m_out.empty(); }

    inline double Mass() const  { return m_mass;  }
    inline double Width() const { return m_width; }

    inline bool Zero() const { return m_zero; }

  };// end of class Current_Base
  
  template <typename SType>
  class Current_Contractor_Base {
  public:

    virtual ~Current_Contractor_Base() {}

    // member functions
    virtual std::complex<SType> 
    SContract(const Current_Base &c,const long unsigned int &hm,
	      const long unsigned int &hp) const = 0;

    virtual void 
    SContract(const Current_Base &c,const LongUIntMap_Matrix &hmp,
	      std::vector<std::complex<SType> > &ress) const = 0;

  };// end of class Current_Base

  template <class SType> inline std::complex<SType> 
  Current_Base::Contract(const Current_Base &c,const long unsigned int &hm,
			 const long unsigned int &hp) const
  { return dynamic_cast<const Current_Contractor_Base<SType>*>
      (this)->SContract(c,hm,hp); }
  template <class SType> inline void 
  Current_Base::Contract(const Current_Base &c,const LongUIntMap_Matrix &hmp,
			 std::vector<std::complex<SType> > &ress) const
  { dynamic_cast<const Current_Contractor_Base<SType>*>
      (this)->SContract(c,hmp,ress); }

  typedef ATOOLS::Getter_Function<Current_Base,Current_Key,
				  std::less<std::string> > Current_Getter;

  std::ostream &operator<<(std::ostream &str,const Current_Base &c);

}// end of namespace COMIX

#endif
