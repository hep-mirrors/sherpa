#ifndef Primitive_Detector_H
#define Primitive_Detector_H

#include "Primitive_Observable_Base.H"
#include "Vector.H"

namespace ATOOLS {
  struct Jet_Data {
    int i,j;
    Vec4D mom;
    double et;
    Jet_Data(int _i, int _j, const Vec4D & _mom, double _et) :
      i(_i),j(_j),mom(_mom),et(_et) {}
  };

  class Primitive_Detector : public Primitive_Observable_Base {
    int   m_mode;
    double m_dR, m_etjet_cut;

    double ** p_et_cal;
    int    ** p_jetno;
    double *  p_costheta;
    double *  p_sintheta;
    double *  p_cosphi;
    double *  p_sinphi;

    double m_res, m_res_ew, m_res_had;
    
    double m_ejet_min, m_etcal_min;

    int    m_neta, m_nphi;
    double m_mineta, m_maxeta, m_delta_eta, m_delta_phi;

    std::vector<Jet_Data> m_jets;

    void SmearEnergy(const Flavour &, double &);
    double PseudoRapidityNAzimuthalAngle(const Vec4D & p, double & phi);

    void Init();
    void FillCalorimeter(const Particle_List &);
    void CalcJets();
    void ClearAll();
    void Print(std::ostream & =std::cout);
  public:
    Primitive_Detector();
    Primitive_Detector(const double,const double,const int,const std::string &);
    ~Primitive_Detector();

    void Evaluate(const Particle_List & pl,double value, int ncount);
    void Evaluate(const Blob_List &,double value, int ncount);
    Primitive_Observable_Base * Copy() const;    
  };

}

#endif
