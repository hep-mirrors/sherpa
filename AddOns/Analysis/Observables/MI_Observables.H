#ifndef Analysis_Observables_MI_Observables_H
#define Analysis_Observables_MI_Observables_H

#include "AddOns/Analysis/Observables/Normalized_Observable.H"

namespace ANALYSIS {

  class MI_Statistics: public Primitive_Observable_Base {  
  private:

    std::vector<ATOOLS::Histogram*> m_scales;

  public:

    MI_Statistics(const size_t scales,const std::string &listname,int mode=0);
    ~MI_Statistics();

    void Evaluate(const ATOOLS::Blob_List &blobs,double weight=1.,double ncount=1);
    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    void Output(const std::string & pname);

    Primitive_Observable_Base &operator+=(const Primitive_Observable_Base &obs);
    Primitive_Observable_Base *Copy() const;

  };// end of class MI_Statistics

  class Multiplicity_vs_JetPT: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    Multiplicity_vs_JetPT(const int type,
			  const double ptmin,const double ptmax,
			  const int nbins,const std::string &jetlist,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Multiplicity_vs_JetPT

  class Scalar_PT_Sum_vs_JetPT: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    Scalar_PT_Sum_vs_JetPT(const int type,
			   const double ptmin,const double ptmax,
			   const int nbins,const std::string &jetlist,
			   const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_PT_Sum_vs_JetPT

  class Scalar_PT_Sum_vs_JetET: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    Scalar_PT_Sum_vs_JetET(const int type,
			   const double ptmin,const double ptmax,
			   const int nbins,const std::string &jetlist,
			   const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_PT_Sum_vs_JetET

  class Multiplicity_vs_PT: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    Multiplicity_vs_PT(const int type,
		       const double ptmin,const double ptmax,
		       const int nbins,const std::string &jetlist,
		       const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Multiplicity_vs_PT

  class Multiplicity_vs_DPhi: public Normalized_Observable {  
  private:

    std::vector<std::string> m_jetlists;

    double m_offset;

  public:

    Multiplicity_vs_DPhi(const int type,
			 const double dphimin,const double dphimax,
			 const int nbins,const double offset,
			 const std::vector<std::string> &jetlists,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Multiplicity_vs_DPhi

  class Scalar_PT_Sum_vs_DPhi: public Normalized_Observable {  
  private:

    std::vector<std::string> m_jetlists;

    double m_offset;

  public:

    Scalar_PT_Sum_vs_DPhi(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_PT_Sum_vs_DPhi

  class Scalar_ET_Sum_vs_DPhi: public Normalized_Observable {  
  private:

    std::vector<std::string> m_jetlists;

    double m_offset;

  public:

    Scalar_ET_Sum_vs_DPhi(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_ET_Sum_vs_DPhi

  class Multiplicity_vs_DEta: public Normalized_Observable {  
  private:

    std::string m_jetlist;
    double      m_offset;

  public:

    Multiplicity_vs_DEta(const int type,
			 const double dphimin,const double dphimax,
			 const int nbins,const double offset,
			 const std::vector<std::string> &jetlists,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Multiplicity_vs_DEta

  class Scalar_PT_Sum_vs_DEta: public Normalized_Observable {  
  private:

    std::string m_jetlist;
    double      m_offset;

  public:

    Scalar_PT_Sum_vs_DEta(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_PT_Sum_vs_DEta

  class Scalar_ET_Sum_vs_DEta: public Normalized_Observable {  
  private:

    std::string m_jetlist;
    double      m_offset;

  public:

    Scalar_ET_Sum_vs_DEta(const int type,
			  const double dphimin,const double dphimax,
			  const int nbins,const double offset,
			  const std::vector<std::string> &jetlists,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_ET_Sum_vs_DEta

  class Multiplicity_vs_Eta: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    Multiplicity_vs_Eta(const int type,
			const double ptmin,const double ptmax,
			const int nbins,const std::string &jetlist,
			const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Multiplicity_vs_Eta

  class Scalar_PT_Sum_vs_Eta: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    Scalar_PT_Sum_vs_Eta(const int type,
			 const double ptmin,const double ptmax,
			 const int nbins,const std::string &jetlist,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_PT_Sum_vs_Eta

  class Scalar_ET_Sum_vs_Eta: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    Scalar_ET_Sum_vs_Eta(const int type,
			 const double ptmin,const double ptmax,
			 const int nbins,const std::string &jetlist,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Scalar_ET_Sum_vs_Eta

  class MIScale_vs_JetPT: public Normalized_Observable {  
  private:

    std::string m_jetlist;

  public:

    MIScale_vs_JetPT(const int type,
			  const double ptmin,const double ptmax,
			  const int nbins,const std::string &jetlist,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Blob_List & blobs,double value,double ncount);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class MIScale_vs_JetPT

}// end of namespace ANALYSIS

#endif
