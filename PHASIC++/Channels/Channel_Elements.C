#include "PHASIC++/Channels/Channel_Elements.H"

#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;

Channel_Elements PHASIC::CE;

double PHASIC::SqLam(double s,double s1,double s2)
{
  double arg(sqr(s-s1-s2)-4.*s1*s2);
  if (arg>0.) return sqrt(arg)/s;
  return 0.;
}

double PHASIC::PeakedDist(double a,double cn,double cxm,double cxp,int k,double ran)
{
  double ce(1.-cn);
  if (ce!=0.) return k*(pow(ran*pow(a+k*cxp,ce)+(1.-ran)*pow(a+k*cxm,ce),1/ce)-a);
  return k*((a+k*cxm)*pow((a+k*cxp)/(a+k*cxm),ran)-a);
}

double PHASIC::PeakedWeight(double a,double cn,double cxm,double cxp,double res,int k,double &ran)
{
  double ce(1.-cn), w;
  if (ce!=0.) {
    double amin=pow(a+k*cxm,ce);
    w=pow(a+k*cxp,ce)-amin;
    ran=(pow(a+k*res,ce)-amin)/w;
    w/=k*ce;
  }
  else {
    double amin=a+k*cxm;
    w=log((a+k*cxp)/amin);
    ran=log((a+k*res)/amin)/w;
    w/=k;
  }
  return w;
}

double Channel_Elements::MasslessPropWeight
(double sexp,double smin,double smax,const double s,double &ran)
{
  if (s<smin || s>smax) {
    msg_Error()<<METHOD<<"(): Value out of bounds: "
	       <<smin<<" .. " <<smax<<" vs. "<<s<< std::endl;
  }
  double w(PeakedWeight(0.,sexp,smin,smax,s,1,ran)/pow(s,-sexp));
  if (IsBad(w)) msg_Error()<<METHOD<<"(): Weight is "<<w<<std::endl;
  return 1./w;
}

double Channel_Elements::MasslessPropMomenta
(double sexp,double smin,double smax, double ran)
{
  double s(PeakedDist(0.,sexp,smin,smax,1,ran));
  if (IsBad(s)) msg_Error()<<METHOD<<"(): Value is "<<s<<std::endl;
  return s;
}

double Channel_Elements::MassivePropWeight
(double m,double g,double smin,double smax,double s,double &ran)
{
  if (s<smin || s>smax)
    msg_Error()<<METHOD<<"(): Value out of bounds: "
	       <<smin<<" .. " <<smax<<" vs. "<<s<< std::endl;
  double m2(m*m), mw(m*g);
  double ymax(atan((smin-m2)/mw)), ymin(atan((smax-m2)/mw));
  double y(atan((s-m2)/mw));
  ran=(y-ymin)/(ymax-ymin);
  double w(mw/((s-m2)*(s-m2)+mw*mw));
  w=(ymin-ymax)/w;
  if (IsBad(w)) msg_Error()<<METHOD<<"(): Weight is "<<w<<std::endl;
  return 1./w;
}

double Channel_Elements::MassivePropMomenta
(double m,double g,double smin,double smax,double ran)
{
  double m2(m*m), mw(m*g), s;
  double ymax(atan((smin-m2)/mw)), ymin(atan((smax-m2)/mw));
  s=m2+mw*tan(ymin+ran*(ymax-ymin));
  if (IsBad(s)) msg_Error()<<METHOD<<"(): Value is "<<s<<std::endl;
  return s;
}

double Channel_Elements::ThresholdWeight
(double sexp,double m,double smin,double smax,double s,double &ran)
{
  if (s<smin || s>smax)
    msg_Error()<<METHOD<<"(): Value out of bounds: "
	       <<smin<<" .. " <<smax<<" vs. "<<s<< std::endl;
  double m2(m*m), sg(sqrt(s*s+m2*m2));
  double sgmin(sqrt(smin*smin+m2*m2)), sgmax(sqrt(smax*smax+m2*m2));
  double w=PeakedWeight(0.,sexp,sgmin,sgmax,sg,1,ran)/(s*pow(sg,-sexp-1.));
  if (IsBad(w)) msg_Error()<<METHOD<<"(): Weight is "<<w<<std::endl;
  return 1./w;
}

double Channel_Elements::ThresholdMomenta
(double sexp,double m,double smin,double smax,double ran)
{
  double m2(m*m);
  double sgmin(sqrt(smin*smin+m2*m2)), sgmax(sqrt(smax*smax+m2*m2));
  double s(sqrt(sqr(PeakedDist(0.,sexp,sgmin,sgmax,1,ran))-m2*m2));
  if (IsBad(s)) msg_Error()<<METHOD<<"(): Value is "<<s<<std::endl;
  return s;
}

void Channel_Elements::Isotropic2Momenta
(Vec4D p,double s1,double s2,Vec4D &p1,Vec4D &p2,double ran1,
 double ran2,double ctmin,double ctmax,const Vec4D &_xref)
{
  double s(dabs(p.Abs2())), rs(sqrt(dabs(s)));
  double e1((s+s1-s2)/rs/2.), m1(sqrt(e1*e1-s1));
  double ct(ctmin+(ctmax-ctmin)*ran1), st(sqrt(1.-ct*ct));
  double phi(2.*M_PI*ran2);
  Vec4D xref(_xref[0]<0.0?-_xref:_xref);
  Vec4D pl(p.PSpat2()?p:Vec4D::ZVEC);
  Poincare cms(p), zax(pl,xref);
  Vec4D n_perp(zax.PT()), l_perp(LT(pl,xref,n_perp));
  l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
  p1=Vec4D(e1,m1*ct*Vec3D(pl)/pl.PSpat());
  p1+=m1*st*(cos(phi)*l_perp+sin(phi)*n_perp);
  cms.BoostBack(p1);
  p2=p-p1;
}

double Channel_Elements::Isotropic2Weight
(const Vec4D &p1, const Vec4D &p2,double &ran1, double &ran2,
 double ctmin, double ctmax,const Vec4D &_xref)
{
  Vec4D p(p1+p2), p1h(p1);
  Vec4D xref(_xref[0]<0.?-_xref:_xref);
  Vec4D pl(p.PSpat2()?p:Vec4D::ZVEC);
  Poincare cms(p), zax(pl,xref);
  Vec4D n_perp(zax.PT()), l_perp(LT(pl,xref,n_perp));
  l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
  cms.Boost(p1h);
  double ct(Vec3D(p1h)*Vec3D(pl)/sqrt(p1h.PSpat2()*pl.PSpat2()));
  double cp(-l_perp*p1), sp(-n_perp*p1), norm(sqrt(cp*cp+sp*sp));
  cp/=norm;
  sp/=norm;
  ran1=(ct-ctmin)/(ctmax-ctmin);
  ran2=atan2(sp,cp)/(2.*M_PI);
  if (ran2<0.) ran2+=1.;
  double w((ctmax-ctmin)/2.);
  w*=M_PI*SqLam(p.Abs2(),p1.Abs2(),p2.Abs2())/2.;
  if (IsBad(w)) msg_Error()<<METHOD<<"(): Weight is "<<w<<"."<<std::endl;
  return 1./w;
}

void Channel_Elements::TChannelMomenta
(Vec4D p1in,Vec4D p2in,Vec4D &p1out,Vec4D &p2out,
 double s1out,double s2out,double mt,double ctexp,
 double ctmax,double ctmin,double ran1,double ran2)
{
  Vec4D pin(p1in+p2in);
  double s(pin.Abs2()), rs(sqrt(dabs(s)));
  double s1in(p1in.Abs2()), s2in(p2in.Abs2());
  double e1in((s+s1in-s2in)/2./rs), m1in(sqrt(e1in*e1in-s1in));
  double e1out((s+s1out-s2out)/2./rs), m1out(sqrt(e1out*e1out-s1out));
  double a=(mt*mt-s1in-s1out+2.*e1in*e1out)/(2.*m1in*m1out);
  if (a<=1.0+1.0e-6) a=1.0+1.0e-6;
  double aminct(PeakedDist(0.,ctexp,a-ctmax,a-ctmin,1,ran1));
  double ct(a-aminct), st(sqrt(1.-ct*ct));
  double phi(2.*M_PI*ran2);
  p1out=Vec4D(e1out,m1out*Vec3D(st*cos(phi),st*sin(phi),ct)); 
  Poincare cms(pin);
  cms.Boost(p1in);
  Poincare zax(p1in,(p1in[3]<0?-1:1) * Vec4D::ZVEC);
  zax.RotateBack(p1out);
  cms.BoostBack(p1out);
  p2out=pin-p1out;
}

double Channel_Elements::TChannelWeight
(const Vec4D &p1in,const Vec4D &p2in,const Vec4D &p1out,const Vec4D &p2out,
 double mt,double ctexp,double ctmax,double ctmin,double &ran1,double &ran2)
{
  Vec4D pin(p1in+p2in), p1inh(p1in), p1outh(p1out);
  double s(pin.Abs2()), rs(sqrt(dabs(s)));
  double s1in(p1in.Abs2()), s2in(p2in.Abs2());
  double s1out(p1out.Abs2()), s2out(p2out.Abs2());
  double e1in((s+s1in-s2in)/2./rs), m1in(sqrt(e1in*e1in-s1in));
  double e1out((s+s1out-s2out)/2./rs), m1out(sqrt(e1out*e1out-s1out));
  double a=(mt*mt-s1in-s1out+2.*e1in*e1out)/(2.*m1in*m1out);
  if (a<=1.0+1.0e-6) a=1.0+1.0e-6;
  Poincare cms(pin);
  cms.Boost(p1inh);
  Poincare zax(p1inh,(p1inh[3]<0?-1:1) * Vec4D::ZVEC);
  cms.Boost(p1outh);
  zax.Rotate(p1outh);
  double pa1(pow(a-ctmax,1.-ctexp));
  double ct(p1outh[3]/p1outh.PSpat());
  if (ct<ctmin || ct>ctmax) {
    msg_Error()<<METHOD<<"(): Error in momentum mapping."<<std::endl;
    ran1=ran2=-1.;
    return 0.;
  }
  ran1=(pow(a-ct,1.-ctexp)-pa1);
  ran1/=(pow(a-ctmin,1.-ctexp)-pa1);
  ran2=asin(p1outh[2]/p1outh.PPerp())/(2.*M_PI);
  if (p1outh[1]<0.) ran2=.5-ran2;
  if (ran2<0.) ran2+=1.;
  double aminct(a-ct);
  double w(PeakedWeight(0.,ctexp,a-ctmax,a-ctmin,aminct,1,ran1));
  w*=m1out*M_PI/(2.*rs)/pow(aminct,-ctexp);
  if (IsBad(w)) msg_Error()<<METHOD<<"(): Weight is "<<w<<"."<<std::endl;
  return 1./w;
}

using namespace std;

double PHASIC::Tj1(double cn,double amcxm,double amcxp,double ran)
{
  double ce(1.-cn);
  if (ce!=0.) return pow(ran*pow(amcxm,ce)+(1.-ran)*pow(amcxp,ce),1./ce);
  if(amcxp>0.) return exp(ran*log(amcxm)+(1.-ran)*log(amcxp));
  return exp(ran*log(-amcxm)+(1.-ran)*log(-amcxp));
}

double PHASIC::Hj1(double cn,double amcxm,double amcxp)
{
  double ce(1.-cn);
  if (ce!=0.) return (pow(amcxp,ce)-pow(amcxm,ce))/ce;
  return log(amcxp/amcxm);
}

void Channel_Elements::CheckMasses(const double &s1, Vec4D &p1,
                                   const double &s2, Vec4D &p2) const {
  if (dabs((s1 - p1.Abs2()) / p1[0]) > 1.e-6) {
    msg_Error() << METHOD << "(): Strong deviation in masses\n"
                << "s1,p1: " << s1 << ";" << p1 << " -> " << p1.Abs2() << " : "
                << dabs(s1 - p1.Abs2()) << ", "
                << "rel = " << dabs((s1 - p1.Abs2()) / p1[0]) << "." << endl;
  }
  if (dabs((s2 - p2.Abs2()) / p2[0]) > 1.e-6) {
    msg_Error() << METHOD << "(): Strong deviation in masses\n"
                << "s2,p2: " << s2 << ";" << p2 << " -> " << p2.Abs2() << " : "
                << dabs(s2 - p2.Abs2()) << ", "
                << "rel = " << dabs((s2 - p2.Abs2()) / p2[0]) << "." << endl;
  }
}

double Channel_Elements::Anisotropic2Weight(const Vec4D& p1,const Vec4D& p2,
					    double& ran1,double& ran2,
					    double ctexp,
					    double ctmin,double ctmax,
					    const Vec4D &_xref)
{
  DEBUG_FUNC("");
  Vec4D p = p1 + p2, xref = _xref[0] < 0. ? -_xref : _xref;
  DEBUG_VAR(p);
  DEBUG_VAR(xref);
  DEBUG_VAR(p1);
  DEBUG_VAR(p2);
  double s = p.Abs2();
  double s1 = p1.Abs2();
  double s2 = p2.Abs2();
  double pabs = sqrt(dabs(s));
  Vec4D p1h;
  p1h[0] = (s + s1 - s2) / pabs / 2.;
  double p1mass = pabs * SqLam(s, s1, s2) / 2.;
  double pmass = sqrt(dabs(p[0] * p[0] - s));
  double a = p[0] * p1h[0] / pmass / p1mass;

  if ((1. >= a) && (a >= 0.))
    a = 1.0000000001;
  double ct = (pabs * p1[0] - p[0] * p1h[0]) / pmass / p1mass;
  if ((ct < ctmin) || (ct > ctmax))
    return 0.;

  double wt =
      1. / (M_PI * SqLam(s, s1, s2) / 4. * pow(a + ct, ctexp) *
            PeakedWeight(a, ctexp, ctmin, ctmax, ct, 1, ran1));
  p1h = p1;
  Vec4D pref(p[0], 0., 0., pmass);
  Poincare Rot(pref, p);
  Rot.RotateBack(p1h);
  Vec4D p1ref=p1h;
  Poincare Boo(pref);
  Boo.Boost(p1h);

  Vec4D zax(p), n_perp(0.0,cross(Vec3D(zax),Vec3D(xref)));
  if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
    msg_Debugging()<<"Set fixed n_perp\n";
    xref=Vec4D(0.,1.,0.,0.);
    n_perp=Vec4D(0.0,cross(Vec3D(zax),Vec3D(xref)));
    if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
      msg_Debugging()<<"Set fixed n_perp\n";
      zax=Vec4D(0.,0.,0.,1.);
      n_perp=Vec4D(0.0,cross(Vec3D(zax),Vec3D(xref)));
    }
  }
  n_perp *= 1.0 / n_perp.PSpat();
  Vec4D l_perp(0., cross(Vec3D(n_perp), Vec3D(zax)));
  l_perp *= 1.0 / l_perp.PSpat();
  DEBUG_VAR(l_perp << " " << p * l_perp << " " << l_perp.Abs2());
  DEBUG_VAR(n_perp << " " << p * n_perp << " " << n_perp.Abs2());

  double cp(-l_perp * p1), sp(-n_perp * p1), norm(sqrt(cp * cp + sp * sp));
  cp /= norm;
  sp /= norm;
  DEBUG_VAR(cp << " " << sp);

  ran2 = asin(sp) / (2. * M_PI);
  if (cp < 0.)
    ran2 = .5 - ran2;
  if (ran2 < 0.)
    ran2 += 1.;

  DEBUG_VAR(ran1 << " " << ran2);
  // ran2        = ::asin(p1h[1]/p1h.PPerp())/(2.*M_PI);
  // if(p1h[2]<0.) ran2=.5-ran2;
  // if (ran2<0.) ran2+=1.;
  if (!(wt > 0) && !(wt < 0))
    msg_Error() << "Anisotropic2Weight produces a nan!" << endl;

  return wt;
}

void Channel_Elements::Anisotropic2Momenta(Vec4D p, double s1, double s2,
                                           Vec4D &p1, Vec4D &p2, double ran1,
                                           double ran2, double ctexp,
                                           double ctmin, double ctmax,
                                           const Vec4D &_xref) {
  DEBUG_FUNC("");
  double s = p.Abs2();
  double pabs = sqrt(dabs(s));
  double e1 = (s + s1 - s2) / pabs / 2.;
  double p1m = pabs * SqLam(s, s1, s2) / 2.;
  double pmass = sqrt(dabs(p[0] * p[0] - s));
  double a = p[0] * e1 / pmass / p1m;
  if ((1. >= a) && (a >= 0.))
    a = 1.0000000001;
  double ct = PeakedDist(a, ctexp, ctmin, ctmax, 1, ran1);
  double st = sqrt(1. - sqr(ct));
  double phi = 2. * M_PI * ran2;

  Poincare cms(p);
  Vec4D xref(_xref[0] < 0.0 ? -_xref : _xref);
  DEBUG_VAR(p);
  DEBUG_VAR(xref);
  Vec4D zax(p), n_perp(0.0, cross(Vec3D(zax), Vec3D(xref)));
  if (n_perp.PSpat2() <= rpa->gen.SqrtAccu()) {
    msg_Debugging() << "Set fixed n_perp\n";
    xref = Vec4D(0., 1., 0., 0.);
    n_perp = Vec4D(0.0, cross(Vec3D(zax), Vec3D(xref)));
    if (n_perp.PSpat2() <= rpa->gen.SqrtAccu()) {
      msg_Debugging() << "Set fixed n_perp\n";
      zax = Vec4D(0., 0., 0., 1.);
      n_perp = Vec4D(0.0, cross(Vec3D(zax), Vec3D(xref)));
    }
  }
  n_perp *= 1.0 / n_perp.PSpat();
  Vec4D l_perp(0., cross(Vec3D(n_perp), Vec3D(zax)));
  l_perp *= 1.0 / l_perp.PSpat();
  p1 = Vec4D(e1, p1m * ct * Vec3D(zax) / zax.PSpat());
  p1 += p1m * st * (cos(phi) * l_perp + sin(phi) * n_perp);
  cms.BoostBack(p1);
  DEBUG_VAR(p << " "
              << p * (p1m * st * (cos(phi) * l_perp + sin(phi) * n_perp)));
  DEBUG_VAR(l_perp << " " << p * l_perp << " " << l_perp.Abs2());
  DEBUG_VAR(n_perp << " " << p * n_perp << " " << n_perp.Abs2());

  double cp(-l_perp * p1), sp(-n_perp * p1), norm(sqrt(cp * cp + sp * sp));
  cp /= norm;
  sp /= norm;
  DEBUG_VAR(cp << " " << sp << " " << cp / cos(phi) - 1. << " "
               << sp / sin(phi) - 1.);

  p2 = p + (-1.) * p1;
  DEBUG_VAR(cos(phi) << " " << sin(phi));
  DEBUG_VAR(ran1 << " " << ran2);
  DEBUG_VAR(p << p.Mass());
  DEBUG_VAR(p1 << p1.Mass());
  DEBUG_VAR(p2 << p2.Mass());

  CheckMasses(s1, p1, s2, p2);
}

double Channel_Elements::BremsstrahlungWeight(double ctexp, double ctmin,
                                              double ctmax, const Vec4D &q,
                                              const Vec4D &p1) {
  Vec4D p = q + p1;
  double sp = p.Abs2();
  double P = Vec3D(p).Abs();
  double sq = q.Abs2();
  double Q = Vec3D(q).Abs();
  double ct = Vec3D(p) * Vec3D(q) / (P * Q);
  if ((ct > ctmax) || (ct < ctmin))
    return 0.;
  double p1m = sqrt(p1.Abs2());
  double ctkin = (2. * p[0] * q[0] - sq - sp + p1m * p1m) / (2. * P * Q);
  if ((0. < ctkin) && (ctkin < 1.))
    ctkin = 1.;
  double amct = ctkin - ct;
  return 1. / (-2. * M_PI * pow(amct, ctexp) *
               Hj1(ctexp, ctkin - ctmin, ctkin - ctmax));
}

void Channel_Elements::BremsstrahlungMomenta(
    Vec4D &p, const double p1mass, const double Eq, const double sq,
    const double ctmin, const double ctmax, const double ctexp, Vec4D &q,
    Vec4D &p1, const double ran1, const double ran2) {
  /* Decay p -> q + p1, q is space-like with energy Eq given from outside
     cos(pq) is constriained by ctmin and ctmax. */
  double sp    = p.Abs2();
  double P     = Vec3D(p).Abs();
  Vec4D  pnorm = Vec4D(1.,0.,0.,1.);
  double Q     = Vec3D(q).Abs();
  double ctkin = (2.*p[0]*Eq-sq-sp+p1mass*p1mass)/(2.*P*Q); 
  if ((0.<ctkin) && (ctkin<1.)) ctkin = 1.;
  double cth = ctkin-Tj1(ctexp,ctkin-ctmin,ctkin-ctmax,ran1);
  double sth = sqrt(1.-cth*cth);
  double cph = cos(2.*M_PI*ran2);
  double sph = sqrt(1.-cph*cph);
  Vec4D qref = Vec4D(Eq,Q*Vec3D(sth*cph,sth*sph,cth)); 
  Poincare rot(pnorm,p);
  q=rot*qref;
  p1 = p+(-1.)*q;  
}

/* Propagators and other 1-dimensional Distributions */

double Channel_Elements::ExponentialMomenta(double sexp, double smin,
                                            double smax, double masses[],
                                            double ran) {
  double sMod = ExponentialDist(sexp, 0, smax - smin, ran);
  double s = sMod + smin;
  if (!(s > 0) && !(s < 0) && s != 0) {
    cout.precision(12);
    cout << "ExpMom : " << sexp << " " << smin << " " << smax << " " << s << " "
         << ran << endl;
    msg_Error() << "ExponentialMomenta produced a nan !" << endl;
  }
  return s;
}
double Channel_Elements::ExponentialWeight(double sexp, double smin,
                                           double smax, double masses[],
                                           const double s, double &ran) {
  if (s < smin || s > smax || smin == smax) {
    ran = -1.;
    return 0.;
  }
  double wt_inv = 0;
  double cbwt =
    PHASIC::ExponentialWeight(sexp, 0, smax - smin); // 1/integral
  wt_inv =
      1 / (exp(-sexp * (s - smin)) * cbwt); // weight = P(s)/I, this is inverse

  if (!(wt_inv > 0) && !(wt_inv < 0) && wt_inv != 0) {
    msg_Error() << "ExponentialWeight produces a nan: " << wt_inv << endl
                << "   smin,s,smax = " << smin << " < " << s << " < " << smax
                << "   sexp = " << sexp << endl;
  }
  return wt_inv;
}

double Channel_Elements::AntennaWeight(double amin, double amax, const double a,
                                       double &ran) {
  if (a < amin || a > amax || amin == amax) {
    ran = -1.;
    return 0.;
  }

  double wt = 1. / (a * (1. - a) *
                    BoundaryPeakedWeight(amin, amax, a, ran));
  if (!(wt > 0) && !(wt < 0) && wt != 0) {
    msg_Error() << "AntennaWeight produces a nan: " << wt << endl
                << "   amin,a,amax = " << amin << " < " << a << " < " << amax
                << endl;
  }
  return wt;
}

double Channel_Elements::AntennaMomenta(double amin, double amax, double ran) {
  double a = BoundaryPeakedDist(amin, amax, ran);
  if (!(a > 0) && !(a < 0) && a != 0)
    msg_Error() << "AntennaMomenta produced a nan !" << endl;
  return a;
}

double Channel_Elements::LLPropWeight(double sexp, double pole, double smin,
                                      double smax, double s, double &ran) {
  if (s < smin || s > smax || smin == smax) {
    ran = -1.;
    return 0.;
  }
  double wt =
      1. / (pow(pole - s, sexp) *
            PeakedWeight(pole, sexp, smin, smax, s, -1, ran));

  if (!(wt > 0) && !(wt < 0) && wt != 0) {
    msg_Error() << " In LL_Weight : " << smin << " < " << s << " < " << smax
                << " ^ " << sexp << ", " << pole << " wt = " << wt << endl
                << "LLPropWeight produces a nan: " << wt << endl;
  }
  return wt;
}

double Channel_Elements::LLPropMomenta(double sexp, double pole, double smin,
                                       double smax, double ran) {
  double s;
  if (smin == smax)
    s = smax;
  else
    s = PeakedDist(pole, sexp, smin, smax, -1, ran);
  if (!(s > 0) && !(s < 0) && s != 0)
    msg_Error() << "LLPropMomenta produced a nan !" << endl;
  if ((s < smin) || (s > smax))
    msg_Error() << "LLPropMomenta out of bounds !" << endl;
  return s;
}

///////////////////////////////////////////////////////////////////////////
double Channel_Elements::GenerateDMRapidityUniform(
    const double masses[], const Double_Container &spinfo,
    Double_Container &xinfo, const double cosXi, const double ran,
    const int mode) {
  double s = spinfo[3];
  double xmin, xmax, x;

  if (ATOOLS::IsEqual(cosXi, -1)) {
    xmin = xinfo[0] = 0.5 + (sqr(masses[0]) - sqr(masses[1])) / (2. * s);
    xmax = xinfo[1] = xinfo[0];
    x = xinfo[0];
    return x;
  } else {
    double test = ATOOLS::Max(masses[0] / sqrt(s), masses[1] / sqrt(s));
    xmin = xinfo[0] = ATOOLS::Max(0.5 - 0.5 * std::abs(cosXi), test);
    xmax = xinfo[1] = ATOOLS::Min(1 - xinfo[0], 1 - test);

    x = xmin + (xmax - xmin) * ran;
  }

  if (ATOOLS::IsZero(x))
    x = 0.;
  if (x < xmin || x > xmax) {
    msg_Error() << METHOD << spinfo << "," << xinfo << ","
                << "): "
                << " X out of bounds ! " << std::endl
                << "  x=" << x << endl;
    // If x is close to any bound, set it to this bound
    if (ATOOLS::IsEqual(x, xmin)) {
      msg_Error() << "Setting x to lower bound  xmin=" << xmin << endl;
      x = xmin;
    }
    if (ATOOLS::IsEqual(x, xmax)) {
      msg_Error() << "Setting x to upper bound xmax=" << xmax << endl;
      x = xmax;
    }
  }
  return x;
}

double Channel_Elements::GenerateDMAngleUniform(const double ran,
                                                const int mode) const {
  double cosxi_min = -1.;
  double cosxi_max = 1.;

  double cosxi = cosxi_min + (cosxi_max - cosxi_min) * ran;

  if (ATOOLS::IsZero(cosxi))
    cosxi = 0.;
  if (cosxi < cosxi_min || cosxi > cosxi_max) {
    msg_Error() << METHOD << " cosXi out of bounds ! " << std::endl
                << "  cosXi=" << cosxi << endl;
    // If x is close to any bound, set it to this bound
    if (ATOOLS::IsEqual(cosxi, cosxi_min)) {
      msg_Error() << "Setting cosXi to lower bound  cosXimin=" << cosxi_min
                  << endl;
      cosxi = cosxi_min;
    }
    if (ATOOLS::IsEqual(cosxi, cosxi_max)) {
      msg_Error() << "Setting cosXi to upper bound cosXimax=" << cosxi_max
                  << endl;
      cosxi = cosxi_max;
    }
  }
  return cosxi;
  // return -1; // TESTING
}

////////////////////////////////////////////////////////////////////////////

// treated from here

double Channel_Elements::GenerateYUniform(const double tau,
                                          const Double_Container &xinfo,
                                          const Double_Container &yinfo,
                                          const double ran,
                                          const int mode) const {
  /*!
    The boundaries for y are
    \begin{align}
    \frac{1}{2}\log\frac{x_{1, min}^2}{\tau} \le y \le
    \log\frac{1}{2}\frac{x_{1, max}^2}{\tau} \frac{1}{2}\log\frac{\tau}{x_{2,
    max}^2} \le y \le \log\frac{1}{2}\frac{\tau}{x_{2, min}^2} \end{align} where
    $x_{1/2, max}$ stem from the corresponding Base or the hard process
    respectively and x_{1, min} = xinfo[0] x_{1, max} = xinfo[2] x_{2, min} =
    xinfo[1] x_{2, max} = xinfo[3]
  */
  double logtau = 0.5 * log(tau);
  if (mode == 1)
    return logtau;
  if (mode == 2)
    return -logtau;
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  double y = ymin + (ymax - ymin) * ran;
  if (ATOOLS::IsZero(y))
    y = 0.;
  if (y < ymin || y > ymax) {
    msg_Error() << "Channel_Elements::GenerateYUniform(" << tau << "," << xinfo
                << "," << yinfo << "): "
                << " Y out of bounds !\n"
                << "   ymin, ymax vs. y : " << ymin << " " << ymax << " vs. "
                << y << "\n";
    // If y is close to any bound, set it to this bound
    if (ATOOLS::IsEqual(y, ymin)) {
      msg_Error() << "Setting y to lower bound  ymin=" << ymin << endl;
      y = ymin;
    }
    if (ATOOLS::IsEqual(y, ymax)) {
      msg_Error() << "Setting y to upper bound ymax=" << ymax << endl;
      y = ymax;
    }
  }
  return y;
}

double Channel_Elements::WeightYUniform(const double tau,
                                        const Double_Container &xinfo,
                                        const Double_Container &yinfo,
                                        double &ran, const int mode) const {
  /*
    See GenerateYUniform for details
  */
  if (mode != 3)
    return 1.;
  double logtau = 0.5 * log(tau);
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  if (yinfo[2] < ymin || yinfo[2] > ymax)
    return 0.0;
  ran = (yinfo[2] - ymin) / (ymax - ymin);
  return (ymax - ymin);
}

const double pre = 1.0;

double Channel_Elements::GenerateYCentral(const double tau,
                                          const Double_Container &xinfo,
                                          const Double_Container &yinfo,
                                          const double ran,
                                          const int mode) const {
  double logtau = 0.5 * log(tau);
  if (mode == 1)
    return logtau;
  if (mode == 2)
    return -logtau;
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  double y = pre * tan(ran * atan(ymax / pre) + (1. - ran) * atan(ymin / pre));
  if (ATOOLS::IsZero(y))
    y = 0.;
  if (y < ymin || y > ymax) {
    msg_Error() << "Channel_Elements::GenerateYCentral(" << tau << "," << xinfo
                << "," << yinfo << "): "
                << " Y out of bounds ! " << std::endl
                << "   ymin, ymax vs. y : " << ymin << " " << ymax << " vs. "
                << y << endl;
    if (ATOOLS::IsEqual(y, ymin)) {
      msg_Error() << "Setting y to lower bound  ymin=" << ymin << endl;
      y = ymin;
    }
    if (ATOOLS::IsEqual(y, ymax)) {
      msg_Error() << "Setting y to upper bound ymax=" << ymax << endl;
      y = ymax;
    }
  }
  return y;
}

double Channel_Elements::WeightYCentral(const double tau,
                                        const Double_Container &xinfo,
                                        const Double_Container &yinfo,
                                        double &ran, const int mode) const {
  if (mode != 3)
    return 1.;
  double logtau = 0.5 * log(tau);
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  if (yinfo[2] < ymin || yinfo[2] > ymax)
    return 0.0;
  double atey = atan(ymin / pre);
  double wt = atan(ymax / pre) - atey;
  ran = (atan(yinfo[2] / pre) - atey) / wt;
  return wt / pre * (pre * pre + yinfo[2] * yinfo[2]);
}

double Channel_Elements::GenerateYForward(
    const double yexponent, const double tau, const Double_Container &xinfo,
    const Double_Container &yinfo, const double ran, const int mode) const {
  double logtau = 0.5 * log(tau);
  if (mode == 1)
    return logtau;
  if (mode == 2)
    return -logtau;
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  double ypeak = ymax - xinfo[3];
  if (yexponent>=1. && ATOOLS::IsEqual(ypeak, ymax)) {
    if (ypeak > 0)
      ypeak *= 1.00000001;
    if (ypeak < 0)
      ypeak /= 1.00000001;
  }
  double y = PeakedDist(ypeak, yexponent, ymin, ymax, -1, ran);
  if (ATOOLS::IsZero(y))
    y = 0.;
  if (y < ymin || y > ymax) {
    msg_Error() << "Channel_Elements::GenerateYForward(" << tau << "," << xinfo
                << "," << yinfo << "): "
                << " Y out of bounds ! " << std::endl
                << "   ymin, ymax vs. y : " << ymin << " " << ymax << " vs. "
                << y << endl;
    if (ATOOLS::IsEqual(y, ymin)) {
      msg_Error() << "Setting y to lower bound  ymin=" << ymin << endl;
      y = ymin;
    }
    if (ATOOLS::IsEqual(y, ymax)) {
      msg_Error() << "Setting y to upper bound ymax=" << ymax << endl;
      y = ymax;
    }
  }
  return y;
}

double Channel_Elements::WeightYForward(const double yexponent,
                                        const double tau,
                                        const Double_Container &xinfo,
                                        const Double_Container &yinfo,
                                        double &ran, const int mode) const {
  if (mode != 3)
    return 1.;
  double logtau = 0.5 * log(tau);
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  if (yinfo[2] < ymin || yinfo[2] > ymax)
    return 0.0;
  double ypeak = ymax - xinfo[3];
  if (yexponent>=1. && ATOOLS::IsEqual(ypeak, ymax)) {
    if (ypeak > 0)
      ypeak *= 1.00000001;
    if (ypeak < 0)
      ypeak /= 1.00000001;
  }

  double wt = PeakedWeight(ypeak, yexponent, ymin, ymax,
                                           yinfo[2], -1, ran) *
              pow(ypeak - yinfo[2], yexponent);
  if (!(wt > 0) && !(wt < 0) && wt != 0) {
    msg_Error() << "WeightYForward produces a nan!" << endl
                << ymax << " " << ymin << " " << yexponent << " " << yinfo[2]
                << " " << xinfo[3] << endl;
  }
  return wt;
}

double Channel_Elements::GenerateYBackward(
    const double yexponent, const double tau, const Double_Container &xinfo,
    const Double_Container &yinfo, const double ran, const int mode) const {
  double logtau = 0.5 * log(tau);
  if (mode == 1)
    return logtau;
  if (mode == 2)
    return -logtau;
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  double ypeak = -ymin - xinfo[2];
  if (yexponent>=1. && ATOOLS::IsEqual(ypeak, -ymin)) {
    if (ypeak > 0)
      ypeak *= 1.00000001;
    if (ypeak < 0)
      ypeak /= 1.00000001;
  }

  double y =
      -PeakedDist(ypeak, yexponent, -ymax, -ymin, -1, ran);
  if (ATOOLS::IsZero(y))
    y = 0.;
  if (y < ymin || y > ymax) {
    msg_Error() << "Channel_Elements::GenerateYBackward(" << logtau << ","
                << xinfo << "," << yinfo << "): ymin, ymax vs. y : " << ymin
                << " " << ymax << " vs. " << y << " (" << (ymin > y) << ", "
                << (ymax < y) << ")\n";
    if (ATOOLS::IsEqual(y, ymin)) {
      msg_Error() << "Setting y to lower bound  ymin=" << ymin << endl;
      y = ymin;
    }
    if (ATOOLS::IsEqual(y, ymax)) {
      msg_Error() << "Setting y to upper bound ymax=" << ymax << endl;
      y = ymax;
    }
  }
  return y;
}

double Channel_Elements::WeightYBackward(const double yexponent,
                                         const double tau,
                                         const Double_Container &xinfo,
                                         const Double_Container &yinfo,
                                         double &ran, const int mode) const {
  if (mode != 3)
    return 1.;
  double logtau = 0.5 * log(tau);
  double ymin = ATOOLS::Max(xinfo[0] - logtau, logtau - xinfo[3]);
  double ymax = ATOOLS::Min(xinfo[2] - logtau, logtau - xinfo[1]);
  ymin = ATOOLS::Max(yinfo[0], ymin);
  ymax = ATOOLS::Min(yinfo[1], ymax);
  if (yinfo[2] < ymin || yinfo[2] > ymax)
    return 0.0;
  double ypeak = -ymin - xinfo[2];
  if (yexponent>=1. && ATOOLS::IsEqual(ypeak, -ymin)) {
    if (ypeak > 0)
      ypeak *= 1.00000001;
    if (ypeak < 0)
      ypeak /= 1.00000001;
  }

  double wt = PeakedWeight(ypeak, yexponent, -ymax, -ymin,
                                           -yinfo[2], -1, ran) *
              pow(ypeak + yinfo[2], yexponent);
  if (!(wt > 0) && !(wt < 0) && wt != 0) {
    msg_Error() << "WeightYBackward produces a nan!" << endl
                << ymax << " " << ymin << " " << yexponent << " " << yinfo[2]
                << " " << xinfo[3] << endl;
  }
  return wt;
}

double PHASIC::ExponentialDist(double ca,double cxm,double cxp,double ran)
{
  double res = 0.;
  if (!IsZero(ca))  res = -log(ran*exp(-ca*cxp) + (1-ran)*exp(-ca*cxm))/ca;
  else              msg_Error()<<"Flat distribution specified, expected exponential"<<endl;
  return res;
}
double PHASIC::ExponentialWeight(double ca,double cxm,double cxp)
{
  double wt = 0;
  if (!IsZero(ca))  wt = ca/(exp(-ca*cxm) - exp(-ca*cxp)); // 1/integral 
  else              msg_Error()<<"Flat distribution specified, expected exponential"<<endl;
  return wt;
}

double PHASIC::BoundaryPeakedDist(double cxm,double cxp,double ran)
  //  1/(x(1-x))
{
  double fxp=1./cxp-1.;
  double fxm=1./cxm-1.;
  double pw = pow(fxm/fxp,ran);
  return pw/(fxm+pw);
}

double PHASIC::BoundaryPeakedWeight(double cxm,double cxp,double res,double &ran)
  //  1/(x(1-x))
{
  double fxp=1./cxp-1.;
  double fxm=1./cxm-1.;
  double wt=log(fxm/fxp);
  ran = log(fxm/(1./res-1.))/wt;
  return wt;
}
