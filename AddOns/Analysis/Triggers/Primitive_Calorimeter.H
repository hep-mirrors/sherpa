#ifndef Analysis_Triggers_Primitive_Calorimeter_H
#define Analysis_Triggers_Primitive_Calorimeter_H

#include "AddOns/Analysis/Triggers/Primitive_Detector.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"

namespace ATOOLS { class Particle_Qualifier_Base; }

namespace ANALYSIS {
  class Primitive_Calorimeter : public Primitive_Detector_Element {
  private:
    double m_mineta, m_maxeta, m_delta_eta, m_delta_phi;
    double *  p_costheta, *  p_sintheta, *  p_cosphi, *  p_sinphi;
    ATOOLS::Particle_Qualifier_Base *p_qualifier;

    void   SmearEnergy(const ATOOLS::Flavour &, double &);
    void   Reset();
  public:
    Primitive_Calorimeter(int,int,std::string) { abort(); }
    Primitive_Calorimeter(const double mineta=-2.5,const double maxeta=2.5,
			  const int neta=50, const int nphi=40,
			  const std::string qualifier="NotLepton");
    ~Primitive_Calorimeter();
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Print(std::ostream & =std::cout);

    double PseudoRapidityNAzimuthalAngle(const ATOOLS::Vec4D & p, double & phi);
    void   GetDimensions(int & neta,int & nphi,
			 double & mineta, double & maxeta); 
    void   GetCosSinTheta(const int i,double & costheta,double & sintheta);  
    void   GetCosSinPhi(const int i,double & cosphi,double & sinphi);  
    Primitive_Detector_Element * Copy() const ;
  };


  inline void Primitive_Calorimeter::GetDimensions(int & neta,int & nphi,
						   double & mineta, double & maxeta) { 
    neta     = m_nx;       nphi     = m_ny;
    mineta   = m_mineta;   maxeta   = m_maxeta;
  }
  inline void Primitive_Calorimeter::GetCosSinTheta(const int i,
						    double & costheta,double & sintheta) {
    if (i>-1&&i<m_nx) { costheta = p_costheta[i]; sintheta = p_sintheta[i]; return; }
    msg_Error()<<"Error in Primitive_Calorimeter::GetCosSinTheta("
		       <<i<<") : "<<std::endl
		       <<"   Out of bounds ("<<m_nx<<"), will return zeroes."<<std::endl;
    costheta = sintheta = 0.;
  }  
  inline void Primitive_Calorimeter::GetCosSinPhi(const int i,
						  double & cosphi,double & sinphi) {
    if (i>-1&&i<m_ny) { cosphi = p_cosphi[i]; sinphi = p_sinphi[i]; return; }
    msg_Error()<<"Error in Primitive_Calorimeter::GetCosSinPhi("
		       <<i<<") : "<<std::endl
		       <<"   Out of bounds ("<<m_ny<<"), will return zeroes."<<std::endl;
    cosphi = sinphi = 0.;
  }
}

#endif
