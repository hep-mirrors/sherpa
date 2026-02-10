#ifndef METOOLS_HadronCurrents_Vector_Decays_H
#define METOOLS_HadronCurrents_Vector_Decays_H

#include "METOOLS/HadronCurrents/Resonance_Base.H"
#include "METOOLS/HadronCurrents/Propagator.H"

namespace METOOLS {
  ////////////////////////////////////////////////////////////////////
  //
  // Vector -> Pseudoscalar + Photon
  //
  // Important: assume REAL photon is always 2nd particle (with mass 0)
  //
  ////////////////////////////////////////////////////////////////////
  class V_PGamma: public Partial_Width_Base {
    double m_lambda2;
    const double ME2(const double & s) const;
    const double PS_2(const double & s) const;

    const double Calculate(const double & s);
  public:
    V_PGamma(const ATOOLS::Flavour & inflav,
	     const std::vector<ATOOLS::Flavour> & outflavs,
	     const double & BR);
    
    void Init(const std::vector<Propagator_Base *> & props,
	      const ATOOLS::axis & Qrange) {};
  };

  ////////////////////////////////////////////////////////////////////
  //
  // Vector -> Pseudoscalar + Pseudoscalar
  //
  ////////////////////////////////////////////////////////////////////
  class V_PP: public Partial_Width_Base {
    double m_lambda2;
    const double ME2(const double & s) const;
    const double PS_2(const double & s) const;

    const double Calculate(const double & s);
  public:
    V_PP(const ATOOLS::Flavour & inflav,
	 const std::vector<ATOOLS::Flavour> & outflavs,
	 const double & BR);
    
    void Init(const std::vector<Propagator_Base *> & props,
	      const ATOOLS::axis & Qrange) {};
  };

  ////////////////////////////////////////////////////////////////////
  //
  // Vector -> Off-shell (Pseudo-)Scalar + Pseudoscalar 
  //
  // Important:
  // - assume on-shell pseudoscalar is always last flavour
  // - assume off-shell particle decays into the first n-1 particle
  //   flavours
  //
  ////////////////////////////////////////////////////////////////////
  class V_PoffP : public Partial_Width_Base {
    ATOOLS::OneDim_Table * p_table;
    double       m_mPmin2, m_mPS, m_mPS2;
    const double Calculate(const double & Q2);
  public:
    V_PoffP(const ATOOLS::Flavour & inflav,
	    const std::vector<ATOOLS::Flavour> & outflavs,
	    const double & BR);
    ~V_PoffP();
    void Init(const std::vector<Propagator_Base * > & props,
	      const ATOOLS::axis & Qrange);
  };  

  ////////////////////////////////////////////////////////////////////
  //
  // Vector -> Off-shell Vector + Pseudoscalar 
  //
  // Note:      I assume this could also be used for
  //            Axial-Vector -> Off-shell Axial-Vector + Pseudoscalar
  // Important: assume off-shell particle is always 1st particle 
  //
  ////////////////////////////////////////////////////////////////////
  class V_VoffP : public Partial_Width_Base {
    ATOOLS::OneDim_Table * p_table;
    double       m_mVmin2, m_mPS, m_mPS2;
    const double Calculate(const double & Q2);
  public:
    V_VoffP(const ATOOLS::Flavour & inflav,
	    const std::vector<ATOOLS::Flavour> & outflavs,
	    const double & BR);
    ~V_VoffP();
    void Init(const std::vector<Propagator_Base * > & props,
	      const ATOOLS::axis & Qrange);
  };  

  ////////////////////////////////////////////////////////////////////
  //
  // Vector -> Off-shell Axial-Vector + Pseudoscalar 
  //
  // Note:      I assume this could also be used for
  //            Axial-Vector -> Off-shell Vector + Pseudoscalar
  // Important: assume off-shell particle is always 1st particle 
  //
  ////////////////////////////////////////////////////////////////////
  class V_AoffP : public Partial_Width_Base {
    ATOOLS::OneDim_Table * p_table;
    double       m_mAmin2, m_mPS, m_mPS2;
    const double Calculate(const double & Q2);
  public:
    V_AoffP(const ATOOLS::Flavour & inflav,
	    const std::vector<ATOOLS::Flavour> & outflavs,
	    const double & BR);
    ~V_AoffP();
    void Init(const std::vector<Propagator_Base * > & props,
	      const ATOOLS::axis & Qrange);
  };  
  

  
  ////////////////////////////////////////////////////////////////////
  //
  // (Axial-)Vector -> Pseudoscalar + Pseudoscalar + Pseudoscalar
  //
  // We assume that this decay is realised through an intermediate
  // vector cecaying into two pseudoscalars, i.e. structures like
  // V -> V^*+P,  V^* -> PP
  // where we include all interference effects.
  //
  ////////////////////////////////////////////////////////////////////
  class V_PPP_Arg_Base {
  protected:
    double m_mi, m_mj, m_mk, m_mi2, m_mj2, m_mk2, m_sum2;
    double m_sij_min, m_skj_min, m_sik_min;
    double m_sum2_ij, m_delta2_ij, m_sum2_kj, m_delta2_kj;
    size_t m_sij_steps, m_skj_steps;

    virtual double dg(const double & Q2,
		      const double & sij,const double & skj) const = 0;
  public:
    V_PPP_Arg_Base(const std::vector<ATOOLS::Flavour> & outflavs);
    virtual void SetProps(const std::vector<Propagator_Base *> & props) = 0;
    virtual void FillTable(ATOOLS::OneDim_Table * table);
  };
  // This helper encodes the structures relevant for a1(1260)
  class V1plus_PPP_Arg : public V_PPP_Arg_Base {
    Propagator_Base * p_propij, * p_propkj;
    double dg(const double & Q2,
	      const double & sij,const double & skj) const;
  public:
    V1plus_PPP_Arg(const std::vector<ATOOLS::Flavour> & outflavs) :
      V_PPP_Arg_Base(outflavs) {}
    void SetProps(const std::vector<Propagator_Base *> & props) ;
  };  

  // This helper encodes the structures relevant for omega(782)
  class V1minus_PPP_Arg : public V_PPP_Arg_Base {
    Propagator_Base * p_propij, * p_propkj, * p_propik;
    double dg(const double & Q2,
	      const double & sij,const double & skj) const;
  public:
    V1minus_PPP_Arg(const std::vector<ATOOLS::Flavour> & outflavs) :
      V_PPP_Arg_Base(outflavs) {}
    void SetProps(const std::vector<Propagator_Base *> & props) ;
  };  

  
  class V_PPP : public Partial_Width_Base {
    ATOOLS::OneDim_Table * p_table;
    V_PPP_Arg_Base       * p_argument;
    void         FillTable() { p_argument->FillTable(p_table); }
    const double Calculate(const double & Q2);
  public:
    V_PPP(const ATOOLS::Flavour & inflav,
	  const std::vector<ATOOLS::Flavour> & outflavs,
	  const double & BR);
    ~V_PPP();
    void Init(const std::vector<Propagator_Base *> & props,
	      const ATOOLS::axis & Qrange);
  };

  ////////////////////////////////////////////////////////////////////
  //
  // (Axial-)Vector -> (Pseudo-)Scalar + Off-shell (Pseudo-)Scalar 
  //
  // We assume that the off-shell particle decays into two scalars,
  // and we may capture interference effects between different
  // resulting amplitudes.
  // 
  ////////////////////////////////////////////////////////////////////
  class V_PPP_scalar : public Partial_Width_Base {
    ATOOLS::OneDim_Table * p_table;
    Propagator_Base      * p_props[3];
    double       m_sum2, m_sum2_ij, m_sum2_jk, m_sum2_ik;
    double       m_sikmin, m_sikmax;
    void         FillTable();
    const double Calculate(const double & Q2);
    const double dg(const double & Q2,const double & sij,const double & sjk) const;
 public:
    V_PPP_scalar(const ATOOLS::Flavour & inflav,
		 const std::vector<ATOOLS::Flavour> & outflavs,
		 const double & BR);
    ~V_PPP_scalar();
    void Init(const std::vector<Propagator_Base *> & props,
	      const ATOOLS::axis & Qrange);
  };  

}

#endif
