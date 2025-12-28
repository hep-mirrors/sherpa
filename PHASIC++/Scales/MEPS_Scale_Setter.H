#ifndef PHASIC__Scales__MEPS_Scale_Setter_H
#define PHASIC__Scales__MEPS_Scale_Setter_H

#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Scales/Color_Setter.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PHASIC++/Scales/Cluster_Definitions.H"
#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {

  class MEPS_Scale_Setter: public Scale_Setter_Base {
  private:

    PDF::Cluster_Definitions_Base *p_clu, *p_qdc;
    Core_Scale_Setter *p_core, *p_uoscale;
    Color_Setter *p_cs;
    ATOOLS::Mass_Selector *p_ms;
    Single_Process *p_sproc;

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

    std::shared_ptr<Color_Integrator> p_ci;

    double m_rsf, m_fsf;
    int    m_cmode, m_nmin;
    int    m_rproc, m_sproc, m_rsproc, m_vproc, m_nproc;

    static int s_nfgsplit, s_allowuo, s_nlocpl;

    int Select(const PDF::ClusterInfo_Vector &ccs,
	       const Int_Vector &on,const int mode=0) const;

    bool CoreCandidate(ATOOLS::Cluster_Amplitude *const ampl) const;
    bool CheckOrdering(ATOOLS::Cluster_Amplitude *const ampl,
		       const int ord=1) const;
    bool CheckSplitting(const PDF::Cluster_Info &cp,
			const int ord) const;
    bool CheckSubEvents(const PDF::Cluster_Config &cc) const;
    void Combine(ATOOLS::Cluster_Amplitude &ampl,
		 const PDF::Cluster_Info &ci) const;

    double SetScales(ATOOLS::Cluster_Amplitude *ampl);
    double Differential(ATOOLS::Cluster_Amplitude *const ampl,
			const int mode=0) const;

    bool ClusterStep(ATOOLS::Cluster_Amplitude *ampl,
		     ATOOLS::ClusterAmplitude_Vector &ampls,
		     const PDF::Cluster_Info &ci,const int ord) const;
    void Cluster(ATOOLS::Cluster_Amplitude *ampl,
		 ATOOLS::ClusterAmplitude_Vector &ampls,
		 const int ord) const;
    void SetCoreScale(ATOOLS::Cluster_Amplitude *const ampl) const;

    double UnorderedScale(ATOOLS::Cluster_Amplitude *const ampl) const;

  public:

    MEPS_Scale_Setter(const Scale_Setter_Arguments &args,
		      const int mode=1);

    ~MEPS_Scale_Setter();

    double Calculate(const ATOOLS::Vec4D_Vector &p,const size_t &mode);

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class MEPS_Scale_Setter

}// end of namespace PHASIC

#endif
