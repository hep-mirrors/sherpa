#include "PDF/SASG/SASGph_Fortran_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Message.H"
#include <unistd.h> 

#include <iostream>

using namespace PDF;
using namespace ATOOLS;

extern "C" {
// CALL SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
void sasgam_(int&, float&, float&, float&, int&, float&, float*);
}

SASGph_Fortran_Interface::SASGph_Fortran_Interface(
		const ATOOLS::Flavour _bunch) {
	m_xmin = 1.e-5;
	m_xmax = 1.;
	m_q2min = .25;
	m_q2max = 1.e6;
	m_nf = 6;

	m_set = "SASG";
	m_bunch = _bunch;
	m_d = m_u = m_s = m_c = m_b = m_g = m_t = 0.;

	for (int i = 1; i < 6; i++) {
		m_partons.insert(Flavour((kf_code) (i)));
		m_partons.insert(Flavour((kf_code) (i)).Bar());
	}
	m_partons.insert(Flavour(kf_gluon));
	m_partons.insert(Flavour(kf_jet));
	m_partons.insert(Flavour(kf_quark));
	m_partons.insert(Flavour(kf_quark).Bar());
}

PDF_Base* SASGph_Fortran_Interface::GetCopy() {
	return new SASGph_Fortran_Interface(m_bunch);
}

void SASGph_Fortran_Interface::CalculateSpec(const double &_x,
		const double &_Q2) {
	float x = _x / m_rescale, Q2 = _Q2;

	float f2photon = 0;
	float f[13];

	// CALL SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
	sasgam_(iset, x, Q2, P2, IP2, f2photon, f);
	m_g = f[6];
	m_d = f[7];
	m_u = f[8];
	m_s = f[9];
	m_c = f[10];
	m_b = f[11];
	m_t = f[12];
}

double SASGph_Fortran_Interface::GetXPDF(const ATOOLS::Flavour &infl) {
	double value = 0.;

	if (infl.Kfcode() == kf_gluon)
		value = m_g;
	else if (infl.Kfcode() == kf_d)
		value = m_d;
	else if (infl.Kfcode() == kf_u)
		value = m_u;
	else if (infl.Kfcode() == kf_s)
		value = m_s;
	else if (infl.Kfcode() == kf_c)
		value = m_c;
	else if (infl.Kfcode() == kf_b)
		value = m_b;
	else if (infl.Kfcode() == kf_t)
		value = m_t;

	return m_rescale * value;
}

double SASGph_Fortran_Interface::GetXPDF(const kf_code &kf, bool anti) {
	double value = 0.;

	if (kf == kf_gluon)
		value = m_g;
	else if (kf == kf_d)
		value = m_d;
	else if (kf == kf_u)
		value = m_u;
	else if (kf == kf_s)
		value = m_s;
	else if (kf == kf_c)
		value = m_c;
	else if (kf == kf_b)
		value = m_b;
	else if (kf == kf_t)
		value = m_t;

	return m_rescale * value;
}

DECLARE_PDF_GETTER(SASGph_Getter);

PDF_Base* SASGph_Getter::operator()(const Parameter_Type &args) const {
	if (!args.m_bunch.IsPhoton())
		return NULL;
	return new SASGph_Fortran_Interface(args.m_bunch);
}

void SASGph_Getter::PrintInfo(std::ostream &str, const size_t width) const {
	str << "SASG photon PDF, see PRD45(1992)3986 and PRD46(1992)1973";
}

SASGph_Getter *p_get_sasg;

extern "C" void InitPDFLib() {
	p_get_sasg = new SASGph_Getter("SASG");
}

extern "C" void ExitPDFLib() {
	delete p_get_sasg;
}
