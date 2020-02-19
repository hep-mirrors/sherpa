#include "CSSHOWER++/Showers/Kinematics_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

Kinematics_Base::Kinematics_Base(): m_evolscheme(0), p_ms(NULL)
{
  Settings& s = Settings::GetMainSettings();
  std::string dipole_string = s["DIPOLES"]["CASE"].Get<std::string>();
  if      (dipole_string == "CS")    m_dipole_case = EXTAMP::DipoleCase::CS;
  else if (dipole_string == "IDa")   m_dipole_case = EXTAMP::DipoleCase::IDa;  // ee > bbWW with mapping a
  else if (dipole_string == "IDb")   m_dipole_case = EXTAMP::DipoleCase::IDb;  // ee > bbWW with mapping b
  else if (dipole_string == "IDin")  m_dipole_case = EXTAMP::DipoleCase::IDin; // pp > bbWW
  else if (dipole_string == "RES")   m_dipole_case = EXTAMP::DipoleCase::RES;  // ee > guu
  else if (dipole_string == "ID")    m_dipole_case = EXTAMP::DipoleCase::ID;   // ee > guu
  else                               m_dipole_case = EXTAMP::DipoleCase::CS;
}


void Kinematics_Base::SetFixVec(Parton *const p,Vec4D mom) const
{
  if (p->GetNext()) SetFixVec(p->GetNext(),mom);
  if (p->FixSpec()==Vec4D()) return;
  Vec4D oldp(p->OldMomentum()), ref(p->FixSpec());
  if (IsEqual(oldp,mom,rpa->gen.SqrtAccu())) {
    p->SetFixSpec(ref);
    p->SetOldMomentum(oldp);
    return;
  }
  Vec4D np(0.0,cross(Vec3D(oldp),Vec3D(mom)));
  if (np.PSpat2()<=1.0e-6) {
    msg_Debugging()<<"Set fixed n_perp\n";
    np=Vec4D(0.0,cross(Vec3D(oldp),Vec3D(1.0,1.0,0.0)));
  }
  np*=1.0/np.PSpat();
  Vec4D lp(0.0,cross(Vec3D(oldp),Vec3D(np)));
  lp*=1.0/lp.PSpat();
  Vec4D pl(0.0,(Vec3D(ref)*Vec3D(lp))*lp);
  Vec4D pn(0.0,(Vec3D(ref)*Vec3D(np))*np);
  Poincare rot(oldp,mom);
  if (oldp.Abs2()>1.0e-3 && mom.Abs2()>1.0e-3) {
    Poincare oldcms(oldp), newcms(mom);
    oldcms.Boost(ref);
    rot.Rotate(ref);
    newcms.BoostBack(ref);
  }
  else {
    double D(oldp*ref), Q(mom[0]), P(mom.PSpat());
    double s((ref-pl-pn).Abs2()), S(mom.Abs2());
    double E(P*sqrt(D*D-s*S)), sgn(Sign((Q*D+E)/(S*ref[0])));
    ref=Vec4D(Q*D+sgn*E,Vec3D(mom)*(P*D+sgn*Q*sqrt(D*D-s*S))/P)/S;
    ref+=rot*(pn+pl);
    if (!IsEqual((ref+mom).Abs2(),(p->FixSpec()+oldp).Abs2(),
		 rpa->gen.SqrtAccu())) {
      Vec4D ref2=Vec4D(Q*D-sgn*E,Vec3D(mom)*(P*D-sgn*Q*sqrt(D*D-s*S))/P)/S;
      ref2+=rot*(pn+pl);
      if (dabs((ref2+mom).Abs2()-(p->FixSpec()+oldp).Abs2())<
	  dabs((ref+mom).Abs2()-(p->FixSpec()+oldp).Abs2())) std::swap<Vec4D>(ref,ref2);
      if (!IsEqual((ref+mom).Abs2(),(p->FixSpec()+oldp).Abs2(),
		   rpa->gen.SqrtAccu())) {
	msg_Error()<<METHOD<<"(): Reconstruction failed {\n"
		   <<"  p_old  = "<<oldp<<" "<<oldp.Abs2()<<"\n"
		   <<"  p_new  = "<<mom<<" "<<mom.Abs2()<<"\n"
		   <<"  k_old  = "<<p->FixSpec()<<" "<<p->FixSpec().Abs2()<<"\n"
		   <<"  k_new1 = "<<ref<<" "<<ref.Abs2()<<"\n"
		   <<"  k_new2 = "<<ref2<<" "<<ref2.Abs2()<<"\n"
		   <<"  s_old  = "<<(p->FixSpec()+oldp).Abs2()
		   <<", s_new1 = "<<(ref+mom).Abs2()
		   <<", s_new2 = "<<(ref2+mom).Abs2()<<"\n}"<<std::endl;
      }
    }
  }
  p->SetFixSpec(ref);
  p->SetOldMomentum(mom);
}

double Kinematics_FF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &mi2,const double &mj2,const double &mk2,
			     const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
                             const ATOOLS::Vec4D pa, const ATOOLS::Vec4D pb,
                             const ATOOLS::Vec4D pi) const
{
  if(m_dipole_case==EXTAMP::IDa)
  if(m_evolscheme==1) return (pa*pi)*(pb*pi)/(pa*pb);
  if(m_evolscheme==2) return pa*pi;

  double pipj=(Q2-mi2-mj2-mk2)*y;
  if (m_evolscheme==0 || m_evolscheme==2) {
    double kt2=pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+mj2;
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    double kt2=pipj*z*(1.0-z);
    if (fla.IsFermion()) kt2=pipj*(flc.IsVector()?(1.0-z):z);
    else if (flc.IsFermion()) kt2=pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+mj2;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

double Kinematics_FF::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &mi2,const double &mj2,const double &mk2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2)) return -1.0;
  double kt2=_kt2;
  if (m_evolscheme==2 || m_evolscheme==3) kt2=kt2-mi2-mj2;
  if (m_evolscheme==0 || m_evolscheme==2) {
    return (kt2/(z*(1.0-z))+(1.0-z)/z*mi2+z/(1.0-z)*mj2)/(Q2-mi2-mj2-mk2);
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    if (fla.IsFermion()) {
      if (flc.IsFermion()) return kt2/z/(Q2-mi2-mj2-mk2);
      return kt2/(1.0-z)/(Q2-mi2-mj2-mk2);
    }
    if (flc.IsFermion()) return kt2/(Q2-mi2-mj2-mk2);
    return kt2/(z*(1.0-z))/(Q2-mi2-mj2-mk2);
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

void V_Calculator::calculate_helpers(const double &Q2,const double &k2,const double &z,
                                     const double &paipb, const double &alpha,
                                     const double &cos2phi)
{
  const double mw2     = sqr(Flavour(24).Mass());
  const double Qprime2 = Q2-mw2;
  m_c =  4*k2*mw2*pow(paipb,2) + 4*k2*pow(paipb,2)*Qprime2 - 2*alpha*k2*paipb*pow(Qprime2,2) +
         pow(paipb,2)*pow(Qprime2,2) - 4*k2*pow(paipb,2)*Qprime2*z -
         2*pow(paipb,2)*pow(Qprime2,2)*z + pow(paipb,2)*pow(Qprime2,2)*pow(z,2);
  m_d = -4*pow(mw2,2)*pow(paipb,2) - 8*mw2*pow(paipb,2)*Qprime2 +
         4*alpha*mw2*paipb*pow(Qprime2,2) -
         8*alpha*cos2phi*mw2*paipb*pow(Qprime2,2) - 4*pow(paipb,2)*pow(Qprime2,2) +
         4*alpha*paipb*pow(Qprime2,3) - 8*alpha*cos2phi*paipb*pow(Qprime2,3) -
         pow(alpha,2)*pow(Qprime2,4) + 8*mw2*pow(paipb,2)*Qprime2*z +
         8*pow(paipb,2)*pow(Qprime2,2)*z - 4*alpha*paipb*pow(Qprime2,3)*z +
         8*alpha*cos2phi*paipb*pow(Qprime2,3)*z - 4*pow(paipb,2)*pow(Qprime2,2)*pow(z,2);
  m_e = -2*mw2*pow(paipb,2)*Qprime2 - 2*pow(paipb,2)*pow(Qprime2,2) +
         alpha*paipb*pow(Qprime2,3) - 2*alpha*cos2phi*paipb*pow(Qprime2,3) +
         2*mw2*pow(paipb,2)*Qprime2*z + 4*pow(paipb,2)*pow(Qprime2,2)*z -
         alpha*paipb*pow(Qprime2,3)*z + 2*alpha*cos2phi*paipb*pow(Qprime2,3)*z -
         2*pow(paipb,2)*pow(Qprime2,2)*pow(z,2);
  m_f =  4*pow(mw2,2)*pow(paipb,2) + 8*mw2*pow(paipb,2)*Qprime2 -
         4*alpha*mw2*paipb*pow(Qprime2,2) + 8*alpha*cos2phi*mw2*paipb*pow(Qprime2,2) +
         4*pow(paipb,2)*pow(Qprime2,2) - 4*alpha*paipb*pow(Qprime2,3) +
         8*alpha*cos2phi*paipb*pow(Qprime2,3) + pow(alpha,2)*pow(Qprime2,4) -
         8*mw2*pow(paipb,2)*Qprime2*z - 8*pow(paipb,2)*pow(Qprime2,2)*z +
         4*alpha*paipb*pow(Qprime2,3)*z - 8*alpha*cos2phi*paipb*pow(Qprime2,3)*z +
         4*pow(paipb,2)*pow(Qprime2,2)*pow(z,2);
  m_g = -(k2*pow(paipb,2)*Qprime2) + k2*pow(paipb,2)*Qprime2*z;
  m_b = complexsqrt((-4.*m_c)/(3.*m_d) + (4.*pow(m_e,2))/pow(m_f,2) - (4.*m_c)/m_f +
        (pow(2,0.3333333333333333)*(16.*pow(m_c,2) - 192.*m_e*m_g -
        192.*m_d*pow(k2,2)*pow(paipb,2)))/(3.*m_d*pow(-128.*pow(m_c,3) +
        2304.*m_c*m_e*m_g + 6912.*m_d*pow(m_g,2) - 4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) -
        6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2) + complexsqrt(-4.*pow(16.*pow(m_c,2) - 192.*m_e*m_g -
        192.*m_d*pow(k2,2)*pow(paipb,2),3) + pow(-128.*pow(m_c,3) + 2304.*m_c*m_e*m_g +
        6912.*m_d*pow(m_g,2) - 4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) -
        6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2),2)),0.3333333333333333)) +
        pow(-128.*pow(m_c,3) + 2304.*m_c*m_e*m_g + 6912.*m_d*pow(m_g,2) -
        4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) -
        6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2) + complexsqrt(-4.*pow(16.*pow(m_c,2) - 192.*m_e*m_g -
        192.*m_d*pow(k2,2)*pow(paipb,2),3) + pow(-128.*pow(m_c,3) + 2304.*m_c*m_e*m_g +
        6912.*m_d*pow(m_g,2) - 4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) -
        6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2),2)),0.3333333333333333)/(3.*pow(
        2,0.3333333333333333)*m_d));
  m_b1 = (pow(2,0.3333333333333333)*(16.*pow(m_c,2) - 192.*m_e*m_g -
         192.*m_d*pow(k2,2)*pow(paipb,2)))/(3.*m_d*pow(-128.*pow(m_c,3) +
         2304.*m_c*m_e*m_g + 6912.*m_d*pow(m_g,2) - 4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) -
         6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2) + complexsqrt(-4.*pow(16.*pow(m_c,2) - 192.*m_e*m_g -
         192.*m_d*pow(k2,2)*pow(paipb,2),3) + pow(-128.*pow(m_c,3) + 2304.*m_c*m_e*m_g +
         6912.*m_d*pow(m_g,2) - 4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) -
         6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2),2)),0.3333333333333333));
  m_b2 = pow(-128.*pow(m_c,3) + 2304.*m_c*m_e*m_g + 6912.*m_d*pow(m_g,2) -
         4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) - 6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2) +
         complexsqrt(-4.*pow(16.*pow(m_c,2) - 192.*m_e*m_g -
         192.*m_d*pow(k2,2)*pow(paipb,2),3) + pow(-128.*pow(m_c,3) + 2304.*m_c*m_e*m_g +
         6912.*m_d*pow(m_g,2) - 4608.*m_c*m_d*pow(k2,2)*pow(paipb,2) -
         6912.*pow(m_e,2)*pow(k2,2)*pow(paipb,2),2)),
         0.3333333333333333)/(3.*pow(2,0.3333333333333333)*m_d);
}

double V_Calculator::GetV1()
{
  const Complex v1 = - m_b/2. - m_e/m_f - 0.5*complexsqrt(4./3*m_c/m_d + 8.*sqr(m_e)/sqr(m_f) -
                     4.*m_c/m_f -
                     (-64.*pow(m_e,3.)/pow(m_f,3.) + 64.*m_c*m_e/sqr(m_f) - 128.*m_g/m_f)/(4.*m_b) -
                     m_b1 - m_b2);
  DEBUG_VAR(v1);
  if(IsZero(v1.imag())) { return v1.real(); }
  else return -1.0;
}

double V_Calculator::GetV2()
{
  const Complex v2 = - m_b/2. - m_e/m_f + 0.5*complexsqrt(4./3*m_c/m_d + 8.*sqr(m_e)/sqr(m_f) -
                     4.*m_c/m_f -
                     (-64.*pow(m_e,3.)/pow(m_f,3.) + 64.*m_c*m_e/sqr(m_f) - 128.*m_g/m_f)/(4.*m_b) -
                     m_b1 - m_b2);
  DEBUG_VAR(v2);
  if(IsZero(v2.imag())) { return v2.real(); }
  else return -1.0;
}

double V_Calculator::GetV3()
{
  const Complex v3 =  m_b/2. - m_e/m_f - 0.5*complexsqrt(4./3*m_c/m_d + 8.*sqr(m_e)/sqr(m_f) -
                      4.*m_c/m_f +
                      (-64.*pow(m_e,3.)/pow(m_f,3.) + 64.*m_c*m_e/sqr(m_f) - 128.*m_g/m_f)/(4.*m_b) -
                      m_b1 - m_b2);
  if(IsZero(v3.imag())) { return v3.real(); }
  DEBUG_VAR(v3);
  else return -1.0;
}

double V_Calculator::GetV4()
{
  const Complex v4 =  m_b/2. - m_e/m_f + 0.5*complexsqrt(4./3*m_c/m_d + 8.*sqr(m_e)/sqr(m_f) -
                      4.*m_c/m_f +
                      (-64.*pow(m_e,3.)/pow(m_f,3.) + 64.*m_c*m_e/sqr(m_f) - 128.*m_g/m_f)/(4.*m_b) -
                      m_b1 - m_b2);
  DEBUG_VAR(v4);
  if(IsZero(v4.imag())) { return v4.real(); }
  else return -1.0;
}

int Kinematics_FF::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{
  switch(m_dipole_case){
    case EXTAMP::CS:
    {
      Parton * spect = split->GetSpect();
      Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
    
      double mij2 = split->Mass2(), mk2 = spect->Mass2();
      if (spect->KScheme()) mk2=p2.Abs2();
      bool nospec(false);
      if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
        mij2=p1.Abs2();
        p2=split->GetPrev()->FixSpec();
        nospec=true;
      }
    
      double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),mi2,mj2,mk2,
    		split->GetFlavour(),flj,1);
      Kin_Args ff(y,split->ZTest(),split->Phi());
      if (ConstructFFDipole(mi2,mj2,mij2,mk2,p1,p2,ff)<0 ||
          !ValidateDipoleKinematics(mi2, mj2, mk2, ff)) return -1;
    
      split->SetMomentum(ff.m_pi);
      if (mi2) SetFixVec(split,ff.m_pi);
      if (!nospec) {
        if (mk2) SetFixVec(spect,ff.m_pk);
        spect->SetMomentum(ff.m_pk);
      }
      else if (!IsEqual(ff.m_pk,p2,1.0e-3))
        msg_Error()<<METHOD<<"(): Error in EW splitting ( y = "<<y
    	       <<" ).\n  Shifted p_k = "
    	       <<p2<<" -> "<<ff.m_pk<<std::endl;
      if (pc==NULL) {
        pc = new Parton(flj,ff.m_pj,pst::FS);
        pc->SetMass2(p_ms->Mass2(flj));
      }
      else {
        if (mj2) SetFixVec(pc,ff.m_pj);
        pc->SetMomentum(ff.m_pj);
      }
      break;
    }
    case EXTAMP::IDa:
    {
      Parton * kinspect = split->GetKinSpect();
      Vec4D paitilde = split->Momentum(), pwtilde = kinspect->Momentum(),
      pb = split->GetSpect()->Momentum();

      const double mw2      = sqr(Flavour(24).Mass());
      const double y_wia    = 1.-split->ZTest();
      const double ztilde_w = 1.-split->YTest(); // YTest() gives vi in ID-case

      Kin_Args ff(y_wia,ztilde_w,split->Phi());

      /* boost into pai+p_ rest-frame in order to be able to use phi_ib and construct additional
         momentum */
      const double Q2 = (paitilde+pwtilde).Abs2();
            ATOOLS::Vec4D pminus = pwtilde - mw2/(Q2-mw2)*paitilde;
      const ATOOLS::Vec4D pboost = paitilde+pminus;
      Poincare bst(pboost);
      bst.Boost(paitilde);
      bst.Boost(pwtilde);
      bst.Boost(pminus);
      bst.Boost(pb);
      ff.m_res    = true;
      ff.m_pb     = pb;
      ff.m_pminus = pminus;
      if (ConstructFFDipole(mw2,0.,mw2,0.,pwtilde,paitilde,ff)<0 ||
          !ValidateDipoleKinematics(mw2, 0., 0., ff)) return -1;
      Vec4D pi = ff.m_pj;
      Vec4D pw = ff.m_pi;
      Vec4D pa = ff.m_pk;
      bst.BoostBack(pi);
      bst.BoostBack(pw);
      bst.BoostBack(pa);

      split->SetMomentum(pa);
      kinspect->SetMomentum(pw);
      // check: on-shellness and momentum conservation of constructed momenta
      if (pc==NULL) {
        pc = new Parton(flj,pi,pst::FS);
      }
      else {
        pc->SetMomentum(pi);
      }
      break;
    }
  }
  return 1;
}

double Kinematics_FI::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &mi2,const double &mj2,const double &ma2,
			     const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc) const
{
  double pipj=-(Q2-ma2-mi2-mj2)*(1.0-y)/y;
  if (m_evolscheme==0 || m_evolscheme==2) {
    double kt2=pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+mj2;
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    double kt2=pipj*z*(1.0-z);
    if (fla.IsFermion()) kt2=pipj*(flc.IsVector()?(1.0-z):z);
    else if (flc.IsFermion()) kt2=pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+mj2;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;

}

double Kinematics_FI::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &mi2,const double &mj2,const double &ma2,
			   const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=mi2+mj2+ma2)) return -1.0;
  double kt2=_kt2;
  if (m_evolscheme==2 || m_evolscheme==3) kt2=kt2-mi2-mj2;
  if (m_evolscheme==0 || m_evolscheme==2) {
    return 1.0/(1.0-(kt2/(z*(1.0-z))+mi2*(1.0-z)/z+mj2*z/(1.0-z))/(Q2-ma2-mi2-mj2));
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    if (fla.IsFermion()) {
      if (flc.IsFermion()) return 1.0/(1.0-kt2/z/(Q2-ma2-mi2-mj2));
      return 1.0/(1.0-kt2/(1.0-z)/(Q2-ma2-mi2-mj2));
    }
    if (flc.IsFermion()) return 1.0/(1.0-kt2/(Q2-ma2-mi2-mj2));
    return 1.0/(1.0-kt2/(z*(1.0-z))/(Q2-ma2-mi2-mj2));
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

int Kinematics_FI::MakeKinematics
(Parton *const split,const double & mi2,const double & mj2,
 const ATOOLS::Flavour &flj,Parton *&pc,const int mode)
{ 
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;

  double ma2 = spect->Mass2(), mij2 = split->Mass2(); 
  bool nospec(false);
  if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
    mij2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    nospec=true;
  }
  
  double Q2((p1-p2).Abs2());
  double y=GetY(Q2,split->KtTest(),split->ZTest(),mi2,mj2,ma2,
		split->GetFlavour(),flj,1);
  Kin_Args fi(1.0-y,split->ZTest(),split->Phi(),8);
  if (ConstructFIDipole(mi2,mj2,mij2,ma2,p1,p2,fi)<0 ||
      !ValidateDipoleKinematics(mi2, mj2, ma2, fi)) return -1;

  split->SetMomentum(fi.m_pi);
  if (mi2) SetFixVec(split,fi.m_pi);
  if (!nospec) {
    if (ma2) SetFixVec(spect,fi.m_pk);
    spect->SetMomentum(fi.m_pk);
  }
  else if (!IsEqual(fi.m_pk,p2,1.0e-3))
    msg_Error()<<METHOD<<"(): Error in EW splitting ( y = "<<y
	       <<" ).\n  Shifted p_k = "
	       <<p2<<" -> "<<fi.m_pk<<std::endl;
  if (pc==NULL) {
    pc = new Parton(flj,fi.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(flj));
  }
  else {
    if (mj2) SetFixVec(pc,fi.m_pj);
    pc->SetMomentum(fi.m_pj);
  }
  
  return 1;
}

double Kinematics_IF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &ma2,const double &mi2,const double &mk2,
			     const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const
{
  double pipj=(Q2-ma2-mi2-mk2)*y/z;
  if (m_evolscheme==0 || m_evolscheme==2) {
    double kt2=-pipj*(1.0-z)-mi2-sqr(1.0-z)*ma2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+ma2;
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    double kt2=-pipj*(1.0-z);
    if (flc.IsFermion()) kt2=-pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+ma2;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

double Kinematics_IF::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mk2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2>=ma2+mi2+mk2)) return -1.0;
  double kt2=_kt2;
  if (m_evolscheme==2 || m_evolscheme==3) kt2=kt2-mi2-ma2;
  if (m_evolscheme==0 || m_evolscheme==2) {
    return -z/(Q2-ma2-mi2-mk2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    if (flc.IsFermion()) return -z/(Q2-ma2-mi2-mk2)*kt2;
    return -z/(Q2-ma2-mi2-mk2)*kt2/(1.0-z);
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

int Kinematics_IF::MakeKinematics
(Parton *const split,const double & ma2,const double & mi2,
 const ATOOLS::Flavour &fli,Parton *&pc,const int mode)
{
  Parton *b(NULL);
  for (PLiter pit(split->GetSing()->begin());pit!=split->GetSing()->end();++pit)
    if ((*pit)->GetType()==pst::IS && *pit!=split) {
      b=*pit;
      break;
    }
  if (b==NULL) THROW(fatal_error,"Corrupted singlet");
  double mb2(b->Mass2());

  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();

  double mk2 = spect->Mass2(), mai2 = split->Mass2(); 
  if (spect->KScheme()) mk2=p2.Abs2();
  bool nospec(false);
  if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
    mai2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    nospec=true;
  }

  double y=GetY((p2-p1).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mk2,
		split->GetFlavour(),fli,1);
  Kin_Args ifp(y,split->ZTest(),split->Phi(),split->Kin());
  if (dabs(y-split->ZTest())<Kin_Args::s_uxeps) ifp.m_mode=1;
  if (ConstructIFDipole(ma2,mi2,mai2,mk2,mb2,p1,p2,b->Momentum(),ifp)<0 ||
      !ValidateDipoleKinematics(ma2, mi2, ifp.m_mk2>=0.0?ifp.m_mk2:mk2, ifp)) return -1;

  split->SetLT(ifp.m_lam);
  ifp.m_lam.Invert();
  split->SetMomentum(ifp.m_lam*ifp.m_pi);
  if (ma2) SetFixVec(split,split->Momentum());
  if (!nospec) {
    spect->SetMomentum(ifp.m_lam*ifp.m_pk);
    if (mk2) SetFixVec(spect,spect->Momentum());
  }
  if (pc==NULL) {
    pc = new Parton(fli,ifp.m_lam*ifp.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(fli));
  }
  else {
    pc->SetMomentum(ifp.m_lam*ifp.m_pj);
    if (mi2) SetFixVec(pc,pc->Momentum());
  }
  
  return 1;
}

double Kinematics_II::GetKT2(const double &Q2,const double &y,const double &z,
			     const double &ma2,const double &mi2,const double &mb2,
			     const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const
{
  double pipj=(Q2-ma2-mi2-mb2)*y/z;
  if (m_evolscheme==0 || m_evolscheme==2) {
    double kt2=pipj*(1.0-z)-mi2-sqr(1.0-z)*ma2;
    if (m_evolscheme==0) return kt2;
    if (m_evolscheme==2) return kt2+mi2+ma2;
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    double kt2=pipj*(1.0-z);
    if (flc.IsFermion()) kt2=pipj;
    if (m_evolscheme==1) return kt2;
    if (m_evolscheme==3) return kt2+mi2+ma2;
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

double Kinematics_II::GetY(const double &Q2,const double &_kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mb2,
			   const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
			   const bool force) const
{
  if (!force && (z<=0.0 || z>=1.0 || Q2<=ma2+mi2+mb2)) return -1.0;
  double kt2=_kt2;
  if (m_evolscheme==2 || m_evolscheme==3) kt2=kt2-mi2-ma2;
  if (m_evolscheme==0 || m_evolscheme==2) {
    return z/(Q2-ma2-mb2-mi2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
  }
  else if (m_evolscheme==1 || m_evolscheme==3) {
    if (flc.IsFermion()) return z/(Q2-ma2-mb2-mi2)*kt2;
    return z/(Q2-ma2-mb2-mi2)*kt2/(1.0-z);
  }
  else THROW(fatal_error, "Not implemented");
  return 0.0;
}

int Kinematics_II::MakeKinematics
(Parton *const split,const double & ma2,const double & mi2,
 const ATOOLS::Flavour &newfl,Parton *&pc,const int mode)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  
  double mai2 = split->Mass2(), mb2 = spect->Mass2();
  bool nospec(false);
  if (mode && split->GetPrev() && split->GetPrev()->KScheme()) {
    mai2=p1.Abs2();
    p2=split->GetPrev()->FixSpec();
    nospec=true;
  }

  double y=GetY((p1+p2).Abs2(),split->KtTest(),split->ZTest(),ma2,mi2,mb2,
		split->GetFlavour(),newfl,1);
  Kin_Args ii(y,split->ZTest(),split->Phi(),split->Kin());
  if (ConstructIIDipole(ma2,mi2,mai2,mb2,p1,p2,ii)<0 ||
      !ValidateDipoleKinematics(ma2, mi2, mb2, ii)) return -1;

  split->SetLT(ii.m_lam);
  ii.m_lam.Invert();
  split->SetMomentum(ii.m_lam*ii.m_pi);
  if (ma2) SetFixVec(split,split->Momentum());
  if (!nospec) {
    spect->SetMomentum(ii.m_lam*ii.m_pk);
    if (mb2) SetFixVec(spect,spect->Momentum());
  }
  if (pc==NULL) {
    pc = new Parton(newfl,ii.m_lam*ii.m_pj,pst::FS);
    pc->SetMass2(p_ms->Mass2(newfl));
  }
  else {
    pc->SetMomentum(ii.m_lam*ii.m_pj);
    if (mi2) SetFixVec(pc,pc->Momentum());
  }

  return 1;
}
