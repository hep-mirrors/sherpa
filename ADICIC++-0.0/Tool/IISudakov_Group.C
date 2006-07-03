//bof
//Version: 4 ADICIC++-0.0/2006/05/28

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
const double IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_colfac=4.0*M_PI;
template<>
const double IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_iieffbas=5.0;
template<>
const double IISudakov<Dipole::iiqbarq,Radiation::igluon>::s_pdfapprox=40;

template class IISudakov<Dipole::iiqbarq,Radiation::gluon>;
template class IISudakov<Dipole::iiqbarq,Radiation::igluon>;




template<> const short IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_x1pow=2;
template<> const short IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_x3pow=3;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_colfac=4.0*M_PI/3.;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_iieffbas=2.0;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::gluon>::s_pdfapprox=10;

template<>
const double IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_colfac=32*M_PI/9.;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_iieffbas=5.0;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::igluon>::s_pdfapprox=100;

template<>
const double IISudakov<Dipole::iiqbarg,Radiation::quark>::s_colfac=8.0*M_PI;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::quark>::s_iieffbas=1.0;
template<>
const double IISudakov<Dipole::iiqbarg,Radiation::quark>::s_pdfapprox=16;

template class IISudakov<Dipole::iiqbarg,Radiation::gluon>;
template class IISudakov<Dipole::iiqbarg,Radiation::igluon>;
template class IISudakov<Dipole::iiqbarg,Radiation::quark>;




template<> const short IISudakov<Dipole::iigq,Radiation::gluon>::s_x1pow=3;
template<> const short IISudakov<Dipole::iigq,Radiation::gluon>::s_x3pow=2;
template<>
const double IISudakov<Dipole::iigq,Radiation::gluon>::s_colfac=4.0*M_PI/3.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::gluon>::s_iieffbas=2.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::gluon>::s_pdfapprox=10;

template<>
const double IISudakov<Dipole::iigq,Radiation::igluon>::s_colfac=32*M_PI/9.;
template<>
const double IISudakov<Dipole::iigq,Radiation::igluon>::s_iieffbas=5.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::igluon>::s_pdfapprox=100;

template<>
const double IISudakov<Dipole::iigq,Radiation::quark>::s_colfac=8.0*M_PI;
template<>
const double IISudakov<Dipole::iigq,Radiation::quark>::s_iieffbas=1.0;
template<>
const double IISudakov<Dipole::iigq,Radiation::quark>::s_pdfapprox=16;

template class IISudakov<Dipole::iigq,Radiation::gluon>;
template class IISudakov<Dipole::iigq,Radiation::igluon>;
template class IISudakov<Dipole::iigq,Radiation::quark>;




template<> const short IISudakov<Dipole::iigg,Radiation::gluon>::s_x1pow=3;
template<> const short IISudakov<Dipole::iigg,Radiation::gluon>::s_x3pow=3;
template<>
const double IISudakov<Dipole::iigg,Radiation::gluon>::s_colfac=4.0*M_PI/3.;
template<>
const double IISudakov<Dipole::iigg,Radiation::gluon>::s_iieffbas=2.0;
template<>
const double IISudakov<Dipole::iigg,Radiation::gluon>::s_pdfapprox=25;

template class IISudakov<Dipole::iigg,Radiation::gluon>;



//=============================================================================



IISudakov_Stats::IISudakov_Stats(Dipole::Type t, const double& inc,
				 const string& n)
  : out(false), incs(1),
    bviols(0), tviols(0),
    blim(inc), tlim(inc),
    min(inc), max(inc),
    type(t), name(n) {}



IISudakov_Stats::IISudakov_Stats(Dipole::Type t, const string& n, bool o,
				 const double& bot, const double& top)
  : out(o), incs(0),
    bviols(0), tviols(0),
    blim(bot), tlim(top),
    min(top), max(bot),
    type(t), name(n) {
  assert(blim>=0.0);
  assert(blim<tlim);
}



IISudakov_Stats::~IISudakov_Stats() {
  //#ifdef SUDAKOV_CALCULATOR_OUTPUT
  if(incs) {
    cout<<string(120,'=')<<"\n"<<__PRETTY_FUNCTION__<<" {\n  "
	<<type<<":"<<name<<"   "<<blim<<"  "<<tlim<<"   |   "
	<<incs<<"   |   "
	<<"min="<<min<<"   max="<<max<<"   "
	<<"viols(-)="<<bviols<<" ("<<1.0*bviols/incs<<")   "
	<<"viols(+)="<<tviols<<" ("<<1.0*tviols/incs<<")   "
	<<"\n}"<<string(119,'=')<<endl;
  }
  //#endif
}



const bool IISudakov_Stats::Include(const double& weight) {
  ++incs;
  bool flag=false;
  if(weight<min) { min=weight; flag=true;}
  if(weight>max) { max=weight; flag=true;}
  if(weight<blim) { ++bviols;}// flag=true;}
  if(weight>tlim) { ++tviols;}// flag=true;}
//#ifdef SUDAKOV_CALCULATOR_OUTPUT
  if(out && flag) cout<<"  IISudakov_Stats("<<type<<":"<<name<<") { "<<weight
		      <<" ==> min="<<min<<",   max="<<max<<";   "
		      <<"viols(-)="<<bviols<<",   viols(+)="<<tviols<<";}\n";
//#endif
  return flag;
}



//=============================================================================





//eof
