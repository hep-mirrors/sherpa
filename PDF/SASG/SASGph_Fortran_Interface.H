#ifndef PDF_SASG_SASGph_Fortran_Interface_H
#define PDF_SASG_SASGph_Fortran_Interface_H

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Phys/Flavour.H"

namespace PDF {
class SASGph_Fortran_Interface: public PDF_Base {
private:
	// PDF parameters
	int iset = 1;
	float P2 = 0;
	int IP2 = 0;
	//std::string m_set,m_path;
	float m_d, m_u, m_s, m_c, m_b, m_g, m_t;
	int m_nf;
public:
	SASGph_Fortran_Interface(const ATOOLS::Flavour);
	PDF_Base* GetCopy();

	void CalculateSpec(const double&, const double&);
	double GetXPDF(const ATOOLS::Flavour&);
	double GetXPDF(const kf_code&, bool);

	double GetXMin() {
		return 1.e-5;
	}
	double GetXMax() {
		return 1.;
	}
	double GetQ2Min() {
		return 0.25;
	}
	double GetQ2Max() {
		return 1.e6;
	}
};

/*!
 \class SASGph_Fortran_Interface
 \brief An interface to the SASG photon distribution.
 */
/*
 * See https://arxiv.org/pdf/hep-ph/9503384.pdf
 */
}

#endif
