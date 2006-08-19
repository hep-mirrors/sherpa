//bof
//Version: 4 ADICIC++-0.0/2006/08/19

//Implementation of Sudakov_Stats.hpp.



#include <iomanip>
#include "Sudakov_Stats.hpp"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



const bool   Sudakov_Stats::s_chk=true;//false;//true;
const size_t Sudakov_Stats::s_mode=1;



namespace ADICIC {
  template<> double Sudakov_Stats::CalcEstimate<0>(const size_t, const double&);
  template<> double Sudakov_Stats::CalcEstimate<1>(const size_t, const double&);
  template<> double Sudakov_Stats::CalcEstimate<2>(const size_t, const double&);
  template<> double Sudakov_Stats::CalcEstimate<3>(const size_t, const double&);
}



//=============================================================================



Sudakov_Stats::Sudakov_Stats(Dipole::Type t, const double& inc,
			     const string& n)
  : out(false), incs(1),
    bviols(0), tviols(0),
    blim(inc), tlim(inc),
    min(inc), max(inc),
    type(t), name(n), entries(), wgtsum() {}



Sudakov_Stats::Sudakov_Stats(Dipole::Type t, const string& n, bool out,
			     const double& bot, const double& top, double exa)
  : out(out), incs(0),
    bviols(0), tviols(0),
    blim(bot), tlim(top),
    min(top), max(bot),
    type(t), name(n), entries(), wgtsum() {
  assert(blim>=0.0);
  assert(blim<tlim);
  if(exa>=0.0) {
    entries.resize(4,Multidouble(100,1.0));
    wgtsum.resize(4,Multidouble(100,exa));
    //pure mean, sigma^2, wmax, max per x.
  }
}



Sudakov_Stats::~Sudakov_Stats() {
  if(incs) {
    cout<<string(120,'=')<<"\n"<<__PRETTY_FUNCTION__<<" {\n  "
	<<type<<":"<<name<<"   "<<blim<<"  "<<tlim<<"   |   "
	<<incs<<"   |   "
	<<"min="<<min<<"   max="<<max<<"   "
	<<"viols(-)="<<bviols<<" ("<<1.0*bviols/incs<<")   "
	<<"viols(+)="<<tviols<<" ("<<1.0*tviols/incs<<")   "
	<<"\n}"<<string(117,'=')
	<<(entries.empty()?"==":"*=")<<endl;
    if(/*s_chk && */!entries.empty()) {
      cout<<setiosflags(ios::left)<<setw(4)<<"x"<<"    "
	    <<setw(9)<<"entries"<<setw(18)<<"mean"<<setw(18)<<"sigma"<<"    "
	    <<setw(9)<<"entries"<<setw(18)<<"wmax"<<"    "
	    <<setw(9)<<"entries"<<"max"<<resetiosflags(ios::left)<<endl;
      cout<<string(120,'-')<<endl;
      for(size_t i=0; i<100; ++i) {
	cout<<setiosflags(ios::left)<<setw(4)<<i/100.0<<"    "
	    <<setw(9)<<entries[0][i]
	    <<setw(18)<<wgtsum[0][i]/Max(entries[0][i],1.0)
	    <<setw(18)<<sqrt(wgtsum[1][i]/Max(entries[0][i]-1,1.0))<<"    "
	    <<setw(9)<<entries[2][i]
	    <<setw(18)<<wgtsum[2][i]/Max(entries[2][i],1.0)<<"    "
	    <<setw(9)<<entries[3][i]<<wgtsum[3][i]
	    <<resetiosflags(ios::left)<<endl;
      }
      cout<<string(120,'=')<<endl;
    }
  }
}



const bool Sudakov_Stats::Include(const double& weight) {
  ++incs;
  bool flag=false;
  if(weight<min) { min=weight; flag=true;}
  if(weight>max) { max=weight; flag=true;}
  if(weight<blim) { ++bviols;}// flag=true;}
  if(weight>tlim) { ++tviols;}// flag=true;}
  if(out && flag) cout<<"  IISudakov_Stats("<<type<<":"<<name<<") { "<<weight
		      <<" ==> min="<<min<<",   max="<<max<<";   "
		      <<"viols(-)="<<bviols<<",   viols(+)="<<tviols<<";}\n";
  return flag;
}



const bool Sudakov_Stats::Include(const double& weight, const double& estim,
				  const double& x) {
  static const bool gate[4]=
    {s_chk || s_mode==0 || s_mode==1, s_chk || s_mode==1,
     s_chk || s_mode==2, s_chk || s_mode==3};
  assert(weight>=0.0 && estim>0.0);
  bool flag=this->Include(weight/estim);
  if(entries.empty()) return flag;
  assert(x>=0.0 && x<=1.0);
  size_t dx=size_t(100.0*x);
  if(gate[0]) {    //mean
    ++entries[0][dx]; wgtsum[0][dx]+=weight;}
  if(gate[1]) {    //sigma^2 (no double bookkeeping of entries)
    wgtsum[1][dx]+=ATOOLS::sqr(weight-wgtsum[0][dx]/entries[0][dx]);}
  if(gate[2]) {    //wmax
    if(wgtsum[2][dx]/entries[2][dx]<weight) {
      ++entries[2][dx]; wgtsum[2][dx]+=weight;}}
  if(gate[3]) {    //max
    if(wgtsum[3][dx]<weight) { ++entries[3][dx]; wgtsum[3][dx]=weight;}}
  return flag;
}



double Sudakov_Stats::GiveEstimate(const double& x, const double& resort) {
  assert(!entries.empty());
  assert(x>=0.0 && x<=1.0);
  size_t dx=size_t(100.0*x);
  return CalcEstimate<s_mode>(dx,resort);
}



//=============================================================================



template<size_t n>
double Sudakov_Stats::CalcEstimate(const size_t dx, const double& resort) {
  cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
      <<"Warning: Method has not been implemented yet!\n"<<endl;
  return resort;
}


namespace ADICIC {
  template<>
  double Sudakov_Stats::CalcEstimate<0>(const size_t dx, const double& resort) {
    return (entries[0][dx]>21.0 ? 2*wgtsum[0][dx]/entries[0][dx] : resort);
  }
  template<>
  double Sudakov_Stats::CalcEstimate<1>(const size_t dx, const double& resort) {
    return (entries[0][dx]>21.0 ?
	    wgtsum[0][dx]/entries[0][dx] +
	    (name.find("GQ")!=string::npos ? 2.2 :
	     name.find("GG")!=string::npos ? 1.7 : 1.0) *
	    2.2*sqrt(wgtsum[1][dx]/(entries[0][dx]-1)) :
	    resort);
  }
  template<>
  double Sudakov_Stats::CalcEstimate<2>(const size_t dx, const double& resort) {
    return (entries[2][dx]>16.0 ? 1.04*wgtsum[2][dx]/entries[2][dx] : resort);
  }
  template<>
  double Sudakov_Stats::CalcEstimate<3>(const size_t dx, const double& resort) {
    //return (entries[3][dx]>3.0 ? wgtsum[3][dx] : resort);
    static double b=8; static double logb=std::log10(b);
    return (entries[3][dx]>3.0 ?
	    logb*wgtsum[3][dx]/std::log10(b+wgtsum[3][dx]) : resort);
  }
}



//=============================================================================





//eof
