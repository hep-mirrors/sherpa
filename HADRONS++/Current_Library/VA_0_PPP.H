#ifndef HADRONS_Current_Library_VA_0_PPP_H
#define HADRONS_Current_Library_VA_0_PPP_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {
  class VA_0_PPP : public Current_Base {     
    class FormFactor {
    protected :
      virtual void FillResonances(const GeneralModel & model)=0;
    public :
      FormFactor(const GeneralModel & model) {}
      virtual ~FormFactor() {}
      virtual Complex F13(const double & s12, const double & s13, const double & s23,
			  const double & Q2) = 0;
      virtual Complex F23(const double & s12, const double & s13, const double & s23,
			  const double & Q2) = 0;
      virtual Complex FV(const double & s12, const double & s13, const double & s23,
			 const double & Q2) = 0;
      virtual Complex FS(const double & s12, const double & s13, const double & s23,
			 const double & Q2) = 0;
    };

  private:
    double       m_global;
    FormFactor * p_ff;
    FormFactor * SelectFormFactor(const ME_Parameters &parameters,const std::string& name);
    FormFactor * ThreePiFF(const ME_Parameters &parameters);
  public:
    VA_0_PPP(const ME_Parameters &parameters,const std::string & name);
    ~VA_0_PPP() { if (p_ff) { delete p_ff; p_ff = 0; } }
    void SetModelParameters( struct GeneralModel _md ) {}
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);

    /////////////////////////////////////////////////////////////////
    // Kuehn-Santamaria form factors
    /////////////////////////////////////////////////////////////////
    class FF_KS_pipipi : public FormFactor {
      int    m_version;
      double m_a1s_norm, m_rhos_norm;
      double m_mpi2, m_ma12;
      std::list<std::pair<ResonanceFlavour, double> > m_a1s, m_rhos;
      void    FillResonances(const GeneralModel & model);
      Complex F(const double & spipi, const double & Q2);
    public:
      FF_KS_pipipi(const GeneralModel & model,const int & version);
      Complex F13(const double & s12, const double & s13, const double & s23,
		  const double & Q2) {
	return F(s13,Q2);
      }
      Complex F23(const double & s12, const double & s13, const double & s23,
		  const double & Q2) {
	return F(s23,Q2);
      }
      Complex FV(const double & s12, const double & s13, const double & s23,
		 const double & Q2) {
	return Complex(0.,0.);
      }
      Complex FS(const double & s12, const double & s13, const double & s23,
		 const double & Q2);
    };
  };
}
#endif
