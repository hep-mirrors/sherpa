#include "CFPSHOWER++/SplittingFunctions/SF_Base.H"
#include "CFPSHOWER++/SplittingFunctions/Functions.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FF_FFFF_Coll : public SF_FF3_Coll {
  protected:
    double m_z1, m_z2, m_xia, m_za_tilde, m_zi_tilde, m_zj_tilde;
    double m_t, m_sai, m_saj, m_sij, m_Q2, m_saij, m_norm, m_cosphi_aj, m_cosphi_ai;
    inline virtual const double B2(const Splitting & split);
    inline virtual const double R() const { return 0.; } 
    inline virtual const double S() const { return 0.; } 
    inline virtual const double I() const { return 0.; } 
    inline virtual const double SymmetryFactor(const Splitting & split) const { return 1.; }
  public:
    FF_FFFF_Coll(const Kernel_Info & info);
    virtual double operator()(const Splitting & split);
    virtual double Integral(const Splitting & split) const;
    virtual double OverEstimate(const Splitting & split) const;
    virtual void   GeneratePoint(Splitting & split) const;
  };

  
  class q_qqbarq: public FF_FFFF_Coll {
  private:
    inline const double R() const { return 0.; } 
    inline const double S() const { return 0.; } 
    inline const double I() const { return 0.; } 
  public:
    q_qqbarq(const Kernel_Info & info) : FF_FFFF_Coll(info) { SetName("q->qqbarq"); }
  };

  class q_qbarqq: public FF_FFFF_Coll {
  private:
    inline const double R() const {
      return Functions::RqqprimeF(m_za_tilde,m_zi_tilde,m_zj_tilde,
				  m_sai,m_saj,m_sij,m_saij);
      return Functions::RqqbarF(m_za_tilde,m_zi_tilde,m_zj_tilde,
				m_sai,m_saj,m_sij,m_saij);
    }
    inline const double S() const {
      return Functions::SqqprimeF(m_za_tilde,m_zi_tilde,m_zj_tilde,
				  m_sai,m_saj,m_sij,m_saij,m_cosphi_aj);
      return Functions::SqqbarF(m_za_tilde,m_zi_tilde,m_zj_tilde,
				m_sai,m_saj,m_sij,m_saij,m_cosphi_aj,m_cosphi_ai);
    }
    inline const double I() const {
      return Functions::DeltaIqqprimeF(m_za_tilde,m_zi_tilde,m_zj_tilde);
      return Functions::DeltaIqqbarF(m_za_tilde,m_zi_tilde,m_zj_tilde);
    }
    inline const double SymmetryFactor(const Splitting & split) const {
      double z = split.z(0), xi = z/split.z(1);
      return (xi-z)/(1.-z);
    }
  public:
    q_qbarqq(const Kernel_Info & info) : FF_FFFF_Coll(info) { SetName("q->qbarqq"); }
  };

  class q_QbarQq: public FF_FFFF_Coll {
  private:
    inline const double R() const {
      return Functions::RqqprimeF(m_za_tilde,m_zi_tilde,m_zj_tilde,
				  m_sai,m_saj,m_sij,m_saij);
    }
    inline const double S() const {
      return Functions::SqqprimeF(m_za_tilde,m_zi_tilde,m_zj_tilde,
				  m_sai,m_saj,m_sij,m_saij,m_cosphi_aj);
    }
    inline const double I() const {
      return Functions::DeltaIqqprimeF(m_za_tilde,m_zi_tilde,m_zj_tilde);
    } 
  public:
    q_QbarQq(const Kernel_Info & info) : FF_FFFF_Coll(info) { SetName("q->QbarQq"); }
  };
}



/*
  The logic in the 1->3 splittings is that they come in two modes:
  either fully differential (m_mode=0), or as subtraction of the end-point 
  only (m_mode=1).  The latter effectively being described by a 1->2 
  kinematics, where one of the two outgoing objects splits into two 
  entirely collinear partons.  Which of the two is picked is defined at 
  random, with equal probability, in "GeneratePoint".

  The splitting function consists of three parts (similar to MC@NLO,
  but without the B and V terms): a collinear subtraction I, a real
  part R, and a subtraction term S. 
 */


using namespace CFPSHOWER;
using namespace ATOOLS;

FF_FFFF_Coll::FF_FFFF_Coll(const Kernel_Info & info) : SF_Base(info) {}

double FF_FFFF_Coll::operator()(const Splitting & split) {
  double value = B2(split), help = value;
  value *= (*split.GetKernel()->GetGauge())(split)/(2.*M_PI); 
  value *= SymmetryFactor(split);
  //msg_Out()<<METHOD<<", B2 * PS = "<<help<<", "
  //	   <<"as = "<<((*split.GetKernel()->GetGauge())(split)/(2.*M_PI))
  //	   <<" --> "<<split.z()<<"*"<<value<<" "
  //	   <<"(sym = "<<SymmetryFactor(split)<<")\n";
  return split.ztilde(0) * value;
}

const double FF_FFFF_Coll::B2(const Splitting & split) {
  if (m_flavs[m_tags[0]].Kfcode()>split.GetKernel()->GetGauge()->NF(split)) return 0.;
  m_z1 = split.z(0); m_z2 = split.z(1); m_xia = m_z1/m_z2;
  m_t  = split.t(0); m_sai = split.sai(); m_Q2 = split.Q2(); m_saij = m_t/m_xia+m_sai;
  m_norm     = m_Q2/(m_Q2-m_saij);
  m_za_tilde = m_norm*m_z1; m_zi_tilde = m_norm*(m_xia-m_z1);
  m_zj_tilde = 1.-m_za_tilde-m_zi_tilde;
  double PSweight = -2.*log(m_za_tilde)/(1.-m_sai/m_saij);
  if (split.IsEndPoint()) return PSweight*I();
  /*
    m_saj = (m_moms[m_tags[0]]+m_moms[m_tags[2]]).Abs2();
    m_sij = (m_moms[m_tags[1]]+m_moms[m_tags[2]]).Abs2();
    m_cosphi_ai = CosPhi(m_moms[m_tags[0]],m_moms[m_tags[1]],m_moms[m_tags[2]],m_specmom);
    m_cosphi_aj = CosPhi(m_moms[m_tags[0]],m_moms[m_tags[2]],m_moms[m_tags[1]],m_specmom);
  */
  // msg_Out()<<METHOD<<"("<<m_name<<", diff): saij = "<<m_saij<<", yaij = "<<(m_saij/m_Q2)<<", "
  // 	   <<"z1 = "<<m_za_tilde<<", z2 = "<<m_zi_tilde<<", z3 = "<<m_zj_tilde<<",\n"
  // 	   <<"     t = "<<split.t()<<", sai = "<<m_sai<<", "
  // 	   <<"saj = "<<m_saj<<", sij = "<<m_sij<<", "
  // 	   <<"cp13 = "<<m_cosphi_aj<<", cp12 = "<<m_cosphi_ai<<", "
  // 	   <<"PS = "<<PSweight<<", \n"
  // 	   <<"     R = "<<R()<<", S = "<<S()<<".\n";
  return PSweight*(R()-S());
}

double FF_FFFF_Coll::Integral(const Splitting & split) const {
  double invkappa = split.Q2()/split.tcut();
  double SFInt    = 20./9. * log(1.+invkappa);
  return SFInt * split.GetKernel()->GetGauge()->OverEstimate(split)/(2.*M_PI); 
}

double FF_FFFF_Coll::OverEstimate(const Splitting & split) const {
  double kappa = split.tcut()/split.Q2();
  double SFEst = 20./9. * 1./(split.z(0)+kappa);
  return SFEst * split.GetKernel()->GetGauge()->OverEstimate(split)/(2.*M_PI); 
}

void FF_FFFF_Coll::GeneratePoint(Splitting & split) const {
  double kappa = split.tcut()/split.Q2();
  double z1    = pow((1.+kappa)/kappa, -ran->Get()) * (1.+kappa) - kappa;
  double phi1  = 2.*M_PI*ran->Get();
  double z2    = pow(z1, ran->Get());
  double help  = ran->Get();
  double sai   = help/(1.-help) * (z2/z1 * split.t(0) + split.m2(0));
  double phi2  = 2.*M_PI*ran->Get();
  bool   endp  = true; //ran->Get()<0.5;
  if (endp) {
    split.SetEndPoint(true);
    //split.SetMode(splitting_mode::coll);
  }
  split.Set_z(0,z1);
  split.Set_sai(sai);
  split.Set_z(1,z2);
  split.Set_phi(0,phi1);
  split.Set_phi(1,phi2);
}

DECLARE_GETTER(FF_FFFF_Coll,"FF_FFFF",SF_Base,Kernel_Info);

SF_Base * ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFFF_Coll>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FF &&
      info.LogType()==log_type::coll &&
      info.GetFlavs().size()==3 &&
      info.GetSplit().IsFermion() && 
      info.GetFlavs()[0].IsFermion() &&
      info.GetFlavs()[1].IsFermion() &&
      info.GetFlavs()[2].IsFermion()) {
    ATOOLS::Flavour split = info.GetSplit(); 
    // no top pair production
    if (split.Kfcode()!=6 &&
	(abs(info.GetFlavs()[1].Kfcode())==6 ||
	 abs(info.GetFlavs()[2].Kfcode())==6) ) return NULL;
    // checking for q -> q qbar q  ==> not yet implemented.
    if (split==info.GetFlavs()[info.TagSequence()[0]]) return NULL;
    // checking for q -> qbar q q, avoid double counting of symmetric q[1]q[2]
    if (split==info.GetFlavs()[info.TagSequence()[0]].Bar() &&
	((split.IsAnti() && !info.GetFlavs()[info.TagSequence()[0]].IsAnti()) ||
	 (!split.IsAnti() && info.GetFlavs()[info.TagSequence()[0]].IsAnti())) &&
	info.TagSequence()[1]<info.TagSequence()[2]){
      return new q_qbarqq(info);
    }
    // checking for q -> Q Qbar q  or   q -> Qbar Q q
    if ((split!=info.GetFlavs()[info.TagSequence()[0]] &&
	 split!=info.GetFlavs()[info.TagSequence()[0]].Bar()) &&
	info.GetFlavs()[info.TagSequence()[0]]==info.GetFlavs()[info.TagSequence()[1]].Bar()) {
      return new q_QbarQq(info);
    }
  }
  return NULL;
}

void ATOOLS::Getter<SF_Base,Kernel_Info,FF_FFFF_Coll>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"FFFF splitting Function (FF)";
}

