#ifndef Analysis_Observables_Event_Shapes_EE_H
#define Analysis_Observables_Event_Shapes_EE_H

#include "AddOns/Analysis/Triggers/Final_Selector.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"

namespace ANALYSIS {

  struct Event_Shape_EE_Data {
    friend std::ostream& operator<<( std::ostream&, const Event_Shape_EE_Data &);
    double        thrust, major, minor, oblateness;
    ATOOLS::Vec3D thrustaxis,majoraxis,minoraxis;
    Event_Shape_EE_Data(double thru,double m1,double m2,double oblate,
			ATOOLS::Vec3D ta,ATOOLS::Vec3D m1a,ATOOLS::Vec3D m2a);
  };

  class Event_Shapes_EE : public Final_Selector {
  private:
    unsigned int                      m_startaxes,m_maxidentaxes;
    double                            m_accuracy;
    std::string                       m_key;
    std::vector<ATOOLS::Vec3D>        m_vectors, m_vectors_save;

    void          Select(const ATOOLS::Particle_List &,double=1.,double ncount=1);
    void          CalculateLinears();

    void          RotateMoms(std::vector<ATOOLS::Vec3D> &,const ATOOLS::Vec3D &);
    ATOOLS::Vec3D NewAxis(const std::vector<ATOOLS::Vec3D> &,const ATOOLS::Vec3D &);
    double        CalculateThrust(const std::vector<ATOOLS::Vec3D> &,const ATOOLS::Vec3D &);
    double        SumP(const std::vector<ATOOLS::Vec3D> &);
    double        SumNP(const std::vector<ATOOLS::Vec3D> &,const ATOOLS::Vec3D &);
    unsigned int  ipow(int,int);
  protected:
    ATOOLS::Vec3D m_thrustaxis,m_majoraxis,m_minoraxis;
    double        m_thrust,m_major,m_minor,m_oblateness;
    double        m_cparameter, m_dparameter;
  public:
    Event_Shapes_EE(const std::string & inlistname,
		    const std::string & outlistname,
		    SP(ATOOLS::Particle_Qualifier_Base) const);
    void Evaluate(const ATOOLS::Blob_List &,double=1.,double ncount=1);
    void Evaluate(const ATOOLS::Particle_List &,double=1.,double ncount=1);
    Analysis_Object * GetCopy() const;

    double        Thrust()     { return m_thrust;     }
    double        Major()      { return m_major;      }
    double        Minor()      { return m_minor;      }
    double        Oblateness() { return m_oblateness; }
    ATOOLS::Vec3D ThrustAxis() { return m_thrustaxis; }
    ATOOLS::Vec3D MajorAxis()  { return m_majoraxis;  }
    ATOOLS::Vec3D MinorAxis()  { return m_minoraxis;  }
  };
}

#endif
