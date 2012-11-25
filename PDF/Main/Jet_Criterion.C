#include "PDF/Main/Jet_Criterion.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace PDF;
using namespace ATOOLS;

double PDF::Qij2(const Vec4D &pi,const Vec4D &pj,const Vec4D &pk,
		 const Flavour &fi,const Flavour &fj,
		 const double &dparam,const int imode)
{
  int mode(imode|2);
  Vec4D npi(pi), npj(pj);
  Flavour nfi((mode&2)?Flavour(kf_gluon):fi);
  Flavour nfj((mode&2)?Flavour(kf_gluon):fj);
  if (npi[0]<0.0) {
    npi=-pi-pj;
    if (!(mode&1)) nfi=fi==fj.Bar()?Flavour(kf_gluon):(fi.IsVector()?fj.Bar():fi);
  }
  else if (npj[0]<0.0) {
    npj=-pj-pi;
    if (!(mode&1)) nfj=fj==fi.Bar()?Flavour(kf_gluon):(fj.IsVector()?fi.Bar():fj);
  }
  if ((fi.IsPhoton() && fj.IntCharge()==0) ||
      (fj.IsPhoton() && fi.IntCharge()==0)) return -1.0;
  double D(1.0);
  if (fi.IsPhoton() || fj.IsPhoton()) {
    if (dparam<0.0) D=-dparam;
    else {
    if (pi[0]<0.0) return pj.PPerp2();
    if (pj[0]<0.0) return pi.PPerp2();
    return Min(pi.PPerp2(), pj.PPerp2())*sqr(pi.DR(pj)/dparam);
    }
  }
  double pipj(dabs(npi*npj)), pipk(dabs(npi*pk)), pjpk(dabs(npj*pk));
  double mti(sqr(Flavour(nfi).Mass())), mtj(sqr(Flavour(nfj).Mass()));
  if (pipj==0.0) {
    if (mti!=0.0 || mtj!=0.0) THROW(fatal_error,"Ill-defined mass term");
  }
  else {
    mti/=2.0*pipj;
    mtj/=2.0*pipj;
  }
  double Cij(nfj.IsVector()?Max(0.0,pipk/(pipj+pjpk)-mti):0.5*pipk/(pipk+pjpk));
  double Cji(nfi.IsVector()?Max(0.0,pjpk/(pipj+pipk)-mtj):0.5*pjpk/(pipk+pjpk));
  return 2.0*dabs(pi*pj)/(Cij+Cji)/D;
}

namespace PDF {

  struct Qij2_Key {
  public:
    
    double m_qij2;
    size_t m_i, m_j, m_k;
    
  public:
    
    inline Qij2_Key(const double &qij2=std::numeric_limits<double>::max(),
		    const size_t &i=0,const size_t &j=0,const size_t &k=0):
      m_qij2(qij2), m_i(i), m_j(j), m_k(k) {}
    
  };// end of struct Qij2_Key

}

double PDF::Qij2Min(const ATOOLS::Vec4D_Vector &p,NLO_subevtlist *const subs)
{
  DEBUG_FUNC(subs->back()->m_pname);
  Qij2_Key qij2;
  const ATOOLS::Flavour *f(subs->back()->p_fl);
  for (size_t n(0);n<subs->size()-1;++n) {
    size_t i((*subs)[n]->m_i), j((*subs)[n]->m_j), k((*subs)[n]->m_k);
    Vec4D pi(i<2?-p[i]:p[i]), pj(j<2?-p[j]:p[j]), pk(k<2?-p[k]:p[k]);
    Flavour fi(i<2?f[i].Bar():f[i]), fj(j<2?f[j].Bar():f[j]);
    double pt2ij=Qij2(pi,pj,pk,fi,fj);
    msg_Debugging()<<"("<<i<<")["<<fi<<"] & ("<<j<<")["<<fj
		   <<"] <-> ("<<k<<"), ptjk = "<<sqrt(pt2ij)<<"\n";
    if (pt2ij<qij2.m_qij2) qij2=Qij2_Key(pt2ij,i,j,k);
  }
  msg_Debugging()<<"q_ij = "<<sqrt(qij2.m_qij2)<<"\n";
  return qij2.m_qij2;
}
