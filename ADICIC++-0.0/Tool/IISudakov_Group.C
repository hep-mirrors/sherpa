//bof
//Version: 4 ADICIC++-0.0/2006/08/12

//Implementation of IISudakov_Group.H.



#include "IISudakov_Group.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "IISudakov_Group.tpt.cc"





//=============================================================================



template class IISudakov_Group<Dipole::iiqbarq>;
template class IISudakov_Group<Dipole::iiqbarg>;
template class IISudakov_Group<Dipole::iigq>;
template class IISudakov_Group<Dipole::iigg>;



//=============================================================================



template<> const short IISudakov<Dipole::iiqbarq,Radiation::gluon>::s_x1pow=2;
template<> const short IISudakov<Dipole::iiqbarq,Radiation::gluon>::s_x3pow=2;
template<>
const double IISudakov<Dipole::iiqbarq,Radiation::gluon>::s_colfac=1.5*M_PI;
template<>
const double IISudakov<Dipole::iiqbarq,Radiation::gluon>::s_iieffbas=2.0;
template<>
const double IISudakov<Dipole::iiqbarq,Radiation::gluon>::s_pdfapprox=4;
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarq,Radiation::gluon>::s_pdf1stat=
Sudakov_Stats(Dipole::iiqbarq,"PDFQQAnn",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarq,Radiation::gluon>::s_pdf2stat=
Sudakov_Stats(Dipole::iiqbarq,"redundantPDFQQAnn",false);


template<>
const double IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_colfac=4.0*M_PI;
template<>
const double IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_iieffbas=5.0;
template<>
const double IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_pdfapprox=40;
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_pdf1stat=
Sudakov_Stats(Dipole::iiqbarq,"PDFQQCom",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_pdf2stat=
Sudakov_Stats(Dipole::iiqbarq,"PDFGQCom",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);


template class IISudakov<Dipole::iiqbarq,Radiation::gluon>;
template class IISudakov<Dipole::iiqbarq,Radiation::igluon>;



//-----------------------------------------------------------------------------



template<> const short IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_x1pow=2;
template<> const short IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_x3pow=3;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_colfac=4.0*M_PI/3.;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_iieffbas=2.0;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_pdfapprox=10;
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_pdf1stat=
Sudakov_Stats(Dipole::iiqbarg,"PDFQQAnn",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_pdf2stat=
Sudakov_Stats(Dipole::iiqbarg,"PDFGGAnn",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);


template<>
const double IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_colfac=32*M_PI/9.;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_iieffbas=5.0;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_pdfapprox=100;
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_pdf1stat=
Sudakov_Stats(Dipole::iiqbarg,"PDFGGCom",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);
template<>
Sudakov_Stats IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_pdf2stat=
Sudakov_Stats(Dipole::iiqbarg,"PDFGQCom",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);


template<>
const double IISudakov<Dipole::iiqbarg,Radiation::quark>::s_colfac=8.0*M_PI;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::quark>::s_iieffbas=1.0;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::quark>::s_pdfapprox=16;


template class IISudakov<Dipole::iiqbarg,Radiation::gluon>;
template class IISudakov<Dipole::iiqbarg,Radiation::igluon>;
template class IISudakov<Dipole::iiqbarg,Radiation::quark>;



//-----------------------------------------------------------------------------



template<> const short IISudakov<Dipole::iigq,Radiation::gluon>::s_x1pow=3;
template<> const short IISudakov<Dipole::iigq,Radiation::gluon>::s_x3pow=2;
template<>
const double IISudakov<Dipole::iigq,Radiation::gluon>::s_colfac=4.0*M_PI/3.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::gluon>::s_iieffbas=2.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::gluon>::s_pdfapprox=10;
template<>
Sudakov_Stats IISudakov<Dipole::iigq,Radiation::gluon>::s_pdf1stat=
Sudakov_Stats(Dipole::iigq,"PDFQQAnn",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);
template<>
Sudakov_Stats IISudakov<Dipole::iigq,Radiation::gluon>::s_pdf2stat=
Sudakov_Stats(Dipole::iigq,"PDFGGAnn",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);


template<>
const double IISudakov<Dipole::iigq,Radiation::igluon>::s_colfac=32*M_PI/9.;
template<>
const double IISudakov<Dipole::iigq,Radiation::igluon>::s_iieffbas=5.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::igluon>::s_pdfapprox=100;
template<>
Sudakov_Stats IISudakov<Dipole::iigq,Radiation::igluon>::s_pdf1stat=
Sudakov_Stats(Dipole::iigq,"PDFGGCom",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);
template<>
Sudakov_Stats IISudakov<Dipole::iigq,Radiation::igluon>::s_pdf2stat=
Sudakov_Stats(Dipole::iigq,"PDFGQCom",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);


template<>
const double IISudakov<Dipole::iigq,Radiation::quark>::s_colfac=8.0*M_PI;
template<>
const double IISudakov<Dipole::iigq,Radiation::quark>::s_iieffbas=1.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::quark>::s_pdfapprox=16;


template class IISudakov<Dipole::iigq,Radiation::gluon>;
template class IISudakov<Dipole::iigq,Radiation::igluon>;
template class IISudakov<Dipole::iigq,Radiation::quark>;



//-----------------------------------------------------------------------------



template<> const short IISudakov<Dipole::iigg,Radiation::gluon>::s_x1pow=3;
template<> const short IISudakov<Dipole::iigg,Radiation::gluon>::s_x3pow=3;
template<>
const double IISudakov<Dipole::iigg,Radiation::gluon>::s_colfac=4.0*M_PI/3.;
template<>
const double IISudakov<Dipole::iigg,Radiation::gluon>::s_iieffbas=2.0;
template<>
const double IISudakov<Dipole::iigg,Radiation::gluon>::s_pdfapprox=25;
template<>
Sudakov_Stats IISudakov<Dipole::iigg,Radiation::gluon>::s_pdf1stat=
Sudakov_Stats(Dipole::iigg,"redundantPDFGGAnn",false);
template<>
Sudakov_Stats IISudakov<Dipole::iigg,Radiation::gluon>::s_pdf2stat=
Sudakov_Stats(Dipole::iigg,"PDFGGAnn",true,0.0,1.0,
	      Sudakov_Calculator::Pdfrule==on ? 0.0 : -1.0);


template class IISudakov<Dipole::iigg,Radiation::gluon>;



//=============================================================================





//eof
