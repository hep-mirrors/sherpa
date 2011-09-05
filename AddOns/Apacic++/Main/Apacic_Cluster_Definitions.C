#include "AddOns/Apacic++/Main/Apacic_Cluster_Definitions.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace APACIC;
using namespace ATOOLS;

Apacic_Cluster_Definitions::Apacic_Cluster_Definitions() {}

CParam Apacic_Cluster_Definitions::KPerp2
(const Cluster_Amplitude &ampl,int i,int j,int k,const Flavour &mo,
 ATOOLS::Mass_Selector *const ms)
{
  if (i<0 && j<0) return CParam(std::numeric_limits<double>::max(),
				std::numeric_limits<double>::max());
  Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j));
  size_t idi(li->Id()), idj(lj->Id());
  Vec4D pi(li->Mom()), pj(lj->Mom());
  Flavour fi(li->Flav()), fj(lj->Flav());
  if (((idi+idj)&3)==3) return CParam(2.0*dabs(pi*pj),2.0*dabs(pi*pj));
  if (fi.IsQuark()&&fj.IsQuark()) 
    return CParam(2.0*dabs(pi*pj),2.0*dabs(pi*pj));
  double cij(1.0), cji(1.0);
  if ((idi&3)==0 && (idj&3)==0) {
    Vec4D pij(pi+pj), pijm(pij.PSpat(),(Vec3D)(-pij));
    double zij(pi[0]/pij[0]);
    if ((zij<0.0 && !IsZero(zij)) || 
	(zij>1.0 && !IsEqual(zij,1.0))) 
      THROW(fatal_error,"Invalid kinematics");
    cij=fi.IsQuark()?1.0:zij/(1.0-zij);
    cji=fj.IsQuark()?1.0:(1.0-zij)/zij;
    if (IsZero(zij) || IsEqual(zij,1.0)) cji=cij=0.0;
  }
  else {
    Vec4D pa;
    if (idi&3) pa=rpa->gen.PBeam(idi&1?0:1);
    else pa=rpa->gen.PBeam(idj&1?0:1);
    double xi, xj;
    if (pa[3]>0.0) {
      xi=pi.PPlus()/pa.PPlus();
      xj=pj.PPlus()/pa.PPlus();
    }
    else {
      xi=pi.PMinus()/pa.PMinus();
      xj=pj.PMinus()/pa.PMinus();
    }
    if ((xi>1.0 && !IsEqual(xi,1.0)) ||
	(xj>1.0 && !IsEqual(xj,1.0))) 
      return CParam(sqr(rpa->gen.Ecms()),sqr(rpa->gen.Ecms()));
    if ((xi<0.0 && !IsZero(xi)) || (xj<0.0 && !IsZero(xj)))
      return CParam(sqr(rpa->gen.Ecms()),sqr(rpa->gen.Ecms()));
    if (idi&3) xi-=xj;
    else if (idj&3) xj-=xi;
    else THROW(fatal_error,"Invalid parton indices");
    cij=fi.IsQuark()?1.0:xi/xj;
    cji=fj.IsQuark()?1.0:xj/xi;
    if (IsZero(xi) || IsZero(xj)) cji=cij=0.0;
  }
  return CParam(2.0*Min(cij,cji)*dabs(pi*pj),
		2.0*Min(cij,cji)*dabs(pi*pj));
}

Vec4D_Vector Apacic_Cluster_Definitions::Combine
(const Cluster_Amplitude &ampl,int i,int j,int k,const Flavour &mo,
 ATOOLS::Mass_Selector *const ms)
{
  if (i>j) std::swap<int>(i,j);
  Vec4D_Vector after(ampl.Legs().size()-1);
  for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
    if (m==(size_t)j) continue;
    after[l]=ampl.Leg(m)->Mom();
    if (m==(size_t)i) after[l]+=ampl.Leg(j)->Mom();
    ++l;
  }
#ifdef USING__Naive_Apacic_Reco_Scheme
  if (i>1) {
    Vec4D split, spec;
    for (size_t l(ampl.NIn());l<after.size();++l) {
      if (l!=i) spec+=after[l];
      else split=after[l];
    }
    Poincare cms(split+spec);
    cms.Boost(split);
    cms.Boost(spec);
    for (size_t l(ampl.NIn());l<after.size();++l) cms.Boost(after[l]);
    Poincare zaxis(split,Vec4D::ZVEC);
    zaxis.Rotate(split);
    zaxis.Rotate(spec);
    for (size_t l(ampl.NIn());l<after.size();++l) zaxis.Rotate(after[l]);
    double sijk((split+spec).Abs2()), Q(sqrt(sijk)/2.0);
    double sk(spec.Abs2()), sij(sqr(ms->Mass(mo)));
    double zhij(0.5+0.5*(sij-sk)/sijk), zhhij(zhij*zhij-sij/sijk);
    double zhk(0.5+0.5*(sk-sij)/sijk), zhhk(zhk*zhk-sk/sijk);
    if (zhhij<0.0 || zhhk<0.0) THROW(fatal_error,"Invalid splitting");
    double zij(zhij+sqrt(zhhij)), bij(sij/(zij*sijk));
    double zk(zhk+sqrt(zhhk)), bk(sk/(zk*sijk));
    Vec4D nsplit((zij+bij)*Q,0.0,0.0,(zij-bij)*Q);
    Vec4D nspec((zk+bk)*Q,0.0,0.0,-(zk-bk)*Q);
    if (after.size()-2>ampl.NIn()) {
      Poincare oldk(spec), newk(nspec);
      for (size_t l(ampl.NIn());l<after.size();++l)
	if (l==i) after[l]=nsplit;
	else {
	  oldk.Boost(after[l]);
	  newk.BoostBack(after[l]);
	}
    }
    else {
      for (size_t l(ampl.NIn());l<after.size();++l)
	if (l==i) after[l]=nsplit;
	else after[l]=nspec;
    }
    for (size_t l(ampl.NIn());l<after.size();++l) {
      zaxis.RotateBack(after[l]);
      cms.BoostBack(after[l]);
    }
  }
#endif
  if (i<2) {
#ifndef USING__Strict_Apacic_Scheme
    size_t a(0), b(1);
    if (ampl.Leg(0)->Id()&2) std::swap<size_t>(a,b);
    Vec4D pa(after[a]), pb(after[b]);
    Poincare newcms(Vec4D(-pa[0]-pb[0],0.0,0.0,-pa[3]-pb[3]));
    Poincare oldcms(-pa-pb);
    oldcms.Boost(pa);
    Poincare zrot(pa,-Vec4D::ZVEC);
    for (size_t i(0);i<after.size();++i) {
      oldcms.Boost(after[i]);
      zrot.Rotate(after[i]);
      newcms.BoostBack(after[i]);
    }
#else
    size_t a(0), b(1);
    if (ampl.Leg(0)->Id()&2) std::swap<size_t>(a,b);
    Vec4D pa(after[a]), pb(after[b]);
    Vec4D p32(-ampl.Leg(0)->Mom()-ampl.Leg(1)->Mom());
    double s32(p32.Abs2()), s12((-pa-pb).Abs2());
    double z(s12/s32), y12(p32.Y());
    if (i==(int)a) y12+=0.5*log(z);
    else y12-=0.5*log(z);
    Poincare oldcms(-pa-pb);
    oldcms.Boost(pa);
    Poincare zrot(pa,-Vec4D::ZVEC);
    Poincare newcms(Vec4D(cosh(y12),0.0,0.0,sinh(y12)));
    for (size_t i(0);i<after.size();++i) {
      oldcms.Boost(after[i]);
      zrot.Rotate(after[i]);
      newcms.BoostBack(after[i]);
    }
#endif
  }
  return after;
}
