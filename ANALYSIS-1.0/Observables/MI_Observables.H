#ifndef MI_Observables_H
#define MI_Observables_H

#include "CXXFLAGS.H"

#ifdef USING__Amisic

#include "Primitive_Observable_Base.H"
#include "../../AMISIC++-1.0.7/Tools/Amisic_Histogram.H"

namespace ANALYSIS {
  class MI_Statistics: public Primitive_Observable_Base {  
  private:

    std::vector<ATOOLS::Histogram*> m_scales;

  public:

    MI_Statistics(const size_t scales,const std::string &listname,int mode=0);
    ~MI_Statistics();

    void Evaluate(const ATOOLS::Blob_List &blobs,double weight=1.,int ncount=1);
    void EndEvaluation(double scale);
    void Output(const std::string & pname);

    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

  };// end of class MI_Statistics

  class Forward_Backward_Eta_Correlation: public Primitive_Observable_Base {  
  private:

    AMISIC::Amisic_Histogram<double> m_etafw, m_etafwsq, m_etafwbw;
  public:

    Forward_Backward_Eta_Correlation(const int type,
				     const double detamin,const double detamax,
				     const int nbins,const std::string &jetlist,
				     const std::string &listname);

    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);

    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Forward_Backward_Eta_Correlation

  class Multiplicity_vs_JetPT: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;
    AMISIC::Amisic_Histogram<double> m_histogram;
  public:

    Multiplicity_vs_JetPT(const int type,
			  const double ptmin,const double ptmax,
			  const int nbins,const std::string &jetlist,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Multiplicity_vs_JetPT

  class Scalar_PT_Sum_vs_JetPT: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_PT_Sum_vs_JetPT(const int type,
			   const double ptmin,const double ptmax,
			   const int nbins,const std::string &jetlist,
			   const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_PT_Sum_vs_JetPT

  class Scalar_PT_Sum_vs_JetET: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_PT_Sum_vs_JetET(const int type,
			   const double ptmin,const double ptmax,
			   const int nbins,const std::string &jetlist,
			   const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_PT_Sum_vs_JetET

  class Multiplicity_vs_PT: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Multiplicity_vs_PT(const int type,
		       const double ptmin,const double ptmax,
		       const int nbins,const std::string &jetlist,
		       const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Multiplicity_vs_PT

  class Multiplicity_vs_DPhi: public Primitive_Observable_Base {  
  private:

    std::vector<std::string> m_jetlists;

    double m_offset;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Multiplicity_vs_DPhi(const int type,
			 const double dphimin,const double dphimax,
			 const int nbins,const double offset,
			 const std::vector<std::string> &jetlists,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Multiplicity_vs_DPhi

  class Scalar_PT_Sum_vs_DPhi: public Primitive_Observable_Base {  
  private:

    std::vector<std::string> m_jetlists;

    double m_offset;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_PT_Sum_vs_DPhi(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_PT_Sum_vs_DPhi

  class Scalar_ET_Sum_vs_DPhi: public Primitive_Observable_Base {  
  private:

    std::vector<std::string> m_jetlists;

    double m_offset;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_ET_Sum_vs_DPhi(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_ET_Sum_vs_DPhi

  class Multiplicity_vs_DEta: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;
    double      m_offset;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Multiplicity_vs_DEta(const int type,
			 const double dphimin,const double dphimax,
			 const int nbins,const double offset,
			 const std::vector<std::string> &jetlists,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Multiplicity_vs_DEta

  class Scalar_PT_Sum_vs_DEta: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;
    double      m_offset;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_PT_Sum_vs_DEta(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_PT_Sum_vs_DEta

  class Scalar_ET_Sum_vs_DEta: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;
    double      m_offset;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_ET_Sum_vs_DEta(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_ET_Sum_vs_DEta

  class Multiplicity_vs_Eta: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Multiplicity_vs_Eta(const int type,
			const double ptmin,const double ptmax,
			const int nbins,const std::string &jetlist,
			const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Multiplicity_vs_Eta

  class Scalar_PT_Sum_vs_Eta: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_PT_Sum_vs_Eta(const int type,
			 const double ptmin,const double ptmax,
			 const int nbins,const std::string &jetlist,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_PT_Sum_vs_Eta

  class Scalar_ET_Sum_vs_Eta: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    Scalar_ET_Sum_vs_Eta(const int type,
			 const double ptmin,const double ptmax,
			 const int nbins,const std::string &jetlist,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Scalar_ET_Sum_vs_Eta

  class MIScale_vs_JetPT: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;

  public:

    MIScale_vs_JetPT(const int type,
			  const double ptmin,const double ptmax,
			  const int nbins,const std::string &jetlist,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Blob_List & blobs,double value,int ncount);
    
    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class MIScale_vs_JetPT

}// end of namespace ANALYSIS

#endif
#endif
