#include "ALPACA/EventGeneration/Elastic_Kernel.H"

#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


Elastic_Kernel::Elastic_Kernel(shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler,
                               shared_ptr<Analysis_Handler> ptr_analysis_handler) :
  p_dyn_quant_handler(ptr_dyn_quant_handler), p_analysis_handler(ptr_analysis_handler),
  m_mes(Simple_MEs()), m_only_gluons(HIPars.OnlyGluons()),
  m_fixed_sigma(HIPars.FixedSigma()),
  mu2_dyn{0.1,0.1}
{
}

Elastic_Kernel::~Elastic_Kernel() 
{
}

double Elastic_Kernel::operator()(const ATOOLS::Flavour in1,  const ATOOLS::Flavour in2, 
		    const ATOOLS::Flavour out1, const ATOOLS::Flavour out2, const double shat) {
  if (shat < 0.) return 0.;

  double sigma;
  
  if(!m_fixed_sigma.first){
    sigma = M_PI*sqr(p_dyn_quant_handler->GetAlphaS())*m_mes(in1,in2,out1,out2,shat)/(shat*shat);
  } else{
    sigma = m_fixed_sigma.second;
  }

  if(sigma < 0 || std::isnan(sigma) || std::isinf(sigma)){
    msg_Out()  << METHOD << ": ERROR, sigma = " << sigma << endl;
  }

  return sigma;   
}


double Elastic_Kernel::FindScatteringXsec(const ATOOLS::Flavour in1,  const ATOOLS::Flavour in2, const double shat) {
  double sigma(0.);
  double id_factor = 1./2.; //Only add id_factor when both in and out are indistinguishable, one factor already in collision kernel
  //id_factor = 1.;
  if(m_fixed_sigma.first){
    id_factor = 1.;
  }

  Flavour u(kf_u), d(kf_d), s(kf_s), g(kf_gluon);
  if (in1.IsQuark() && in2.IsQuark()) {
    if (in1==in2) 
      sigma += ((*this)(u,u,u,u,shat))*id_factor; //*2?
      //sigma += ((*this)(u,u,u,u,shat)); //*2?
    else if (in1==in2.Bar()) {
      sigma += (*this)(u,u.Bar(),u,u.Bar(),shat);
      sigma += (*this)(u,u.Bar(),d,d.Bar(),shat);
      sigma += (*this)(u,u.Bar(),s,s.Bar(),shat); //new
      //double temp_sigma = ((*this)(u,u.Bar(),g,g,shat))*id_factor;
      double temp_sigma = ((*this)(u,u.Bar(),g,g,shat));
      if(temp_sigma >= 0){
        sigma += temp_sigma; //*2?
      } else{
        sigma = -pow(10,10);
      }
    }
    else if (in1!=in2)       
      sigma += (*this)(u,d,u,d,shat);
      //sigma += (*this)(u,s,u,s,shat); //remove?
  }
  else if ((in1.IsGluon() && in2.IsQuark()) || (in1.IsQuark() && in2.IsGluon()))
    sigma += (*this)(u,g,u,g,shat);
  else if (in1.IsGluon() && in2.IsGluon()) {
    sigma += ((*this)(g,g,g,g,shat))*id_factor; //*2?
    if(!m_only_gluons){
      sigma += (*this)(g,g,u,u.Bar(),shat);
      sigma += (*this)(g,g,d,d.Bar(),shat); //new
      sigma += (*this)(g,g,s,s.Bar(),shat); //new
    }
  }
  else THROW(fatal_error,"Did not find process");
  msg_Debugging()<<METHOD<<" for "<<in1<<" + "<<in2<<" at shat = "<<shat<<" : cross section is "<<sigma<<" /GeV^2.\n";
  
  return sigma;
}


double Elastic_Kernel::SelectThat(const ATOOLS::Flavour in1, const ATOOLS::Flavour in2, const double shat,
				  ATOOLS::Flavour & out1, ATOOLS::Flavour & out2) {
  //msg_Out()<<METHOD<<" for "<<in1<<" + "<<in2<<" at shat = "<<shat<<".\n";
  double that_fixed_sigma = 1.; //If fixed sigma is true, return 1
  double id_factor = 1./2.; //Only add id_factor when both in and out are indistinguishable, one factor already in collision kernel
  if(m_fixed_sigma.first){
    id_factor = 1.;
  }
  Flavour u(kf_u), d(kf_d), s(kf_s), g(kf_gluon);
  if (in1.IsQuark() && in2.IsQuark()) {
    if (in1==in2) {
      out1=in1;
      out2=in2;
      //msg_Out() << "\nqiqi->qiqi" << endl;
      if(!m_fixed_sigma.first){
        return SelectThat_qiqi_qiqi(shat); //*2?
      } else{
        return that_fixed_sigma;
      }
    }
    else if (in1==in2.Bar()) {
      double sigmatot(0.), sigma1, sigma2, sigma3;
      sigmatot += sigma1 = (*this)(u,u.Bar(),u,u.Bar(),shat);
      sigmatot += sigma2 = (*this)(u,u.Bar(),d,d.Bar(),shat);
      sigmatot += sigma3 = (*this)(u,u.Bar(),s,s.Bar(),shat); //new
      //sigmatot += (*this)(u,u.Bar(),g,g,shat)*id_factor; //*2?
      sigmatot += (*this)(u,u.Bar(),g,g,shat); //*2?

      double R(ran->Get());
      if (R < sigma1/sigmatot) {
        out1=in1;
        out2=in2;
        //msg_Out() << "\nqiqbi->qiqbi" << endl;
        if(!m_fixed_sigma.first){
          return SelectThat_qiqbi_qiqbi(shat);
        } else{
          return that_fixed_sigma;
        }
      }
      else if (R < (sigma1+sigma2+sigma3)/sigmatot) {
        Flavour outflav;
        do {
          double R2(ran->Get());
          if (R2 < 1./3.)      outflav = u;
          else if (R2 < 2./3.) outflav = d;
          else                 outflav = s;
        } while (outflav == in1 || outflav==in2);
        out1=outflav;
        out2=outflav.Bar();
        //msg_Out() << "\nqiqbi->qjqbj" << endl;
        if(!m_fixed_sigma.first){
          return SelectThat_qiqbi_qjqbj(shat);
        } else{
          return that_fixed_sigma;
        }
      }
      else {
        out1=Flavour(kf_gluon);
        out2=Flavour(kf_gluon);
        //msg_Out() << "\nqiqbi->gg" << endl;
        if(!m_fixed_sigma.first){
          return SelectThat_qiqbi_gg(shat);
        } else{
          return that_fixed_sigma;
        }
      }
    }
    else if (in1!=in2) {
      out1=in1;
      out2=in2;
      //msg_Out() << "\nqiqj->qiqj" << endl;
      if(!m_fixed_sigma.first){
        return SelectThat_qiqj_qiqj(shat);
      } else{
        return that_fixed_sigma;
      }
    }
  }
  else if ((in1.IsGluon() && in2.IsQuark()) || (in1.IsQuark() && in2.IsGluon())) {
      out1=in1;
      out2=in2;
      //msg_Out() << "\nqig->qig" << endl;
      if(!m_fixed_sigma.first){
        return SelectThat_qig_qig(shat);
      } else{
        return that_fixed_sigma;
      }
  }
  else if (in1.IsGluon() && in2.IsGluon()) {
    double sigmatot(0.), sigma1;
    sigmatot += sigma1 = (*this)(g,g,g,g,shat)*id_factor; //*2?
    //sigmatot += sigma1 = (*this)(g,g,g,g,shat); //*2?
    if(!m_only_gluons){
      sigmatot += (*this)(g,g,u,u.Bar(),shat);
      sigmatot += (*this)(g,g,d,d.Bar(),shat); //new
      sigmatot += (*this)(g,g,s,s.Bar(),shat); //new
    }
    double R(ran->Get());
    if(sigmatot > 0. && !std::isnan(sigmatot) && !std::isinf(sigmatot)){
      if (R < sigma1/sigmatot) {
        out1=in1;
        out2=in2;
        //msg_Out() << "\ngg->gg" << endl;
        if(!m_fixed_sigma.first){
          return SelectThat_gg_gg(shat);
        } else{
          return that_fixed_sigma;
        }
      } else {
        Flavour outflav;
        double R2(ran->Get());
        if (R2 < 1./3.)      outflav = u;
        else if (R2 < 2./3.) outflav = d;
        else                 outflav = s;
        out1=outflav;
        out2=outflav.Bar();
        //msg_Out() << "\ngg->qiqbi" << endl;
        if(!m_fixed_sigma.first){
          return SelectThat_gg_qiqbi(shat);
        } else{
          return that_fixed_sigma;
        }
      }
    } else{
      msg_Out() << METHOD << ": ERROR: " << sigmatot << endl;
      return 1.;
    }
  }
  else{
    msg_Out() << "Fatal error \n";
    THROW(fatal_error,"Did not find process");
  }

  return 1.;
}

double Elastic_Kernel::SelectThat_qiqi_qiqi(const double shat) {
  int ntrial(0);
  Flavour u(kf_u);
  double that,R,weight,mefull,meover;
  do {
    R=ran->Get();
    /*
    if (ran->Get() < 0.5) {
      that = mu2_dyn[1]*shat*(R-1.)/(mu2_dyn[1]+R*shat);
    }
    else {
      that = (mu2_dyn[1]+shat)*shat*(R-1.)/(mu2_dyn[1]-shat*(R-1.));
    }
    if (that < -shat || that > 0.)
      msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = 8.*shat*shat/9.*(1./sqr(that-mu2_dyn[1]) + 1./sqr(-that-shat-mu2_dyn[1]));
    meover = (1./36.)*(32*shat*shat*(1./sqr(that-mu2_dyn[1]) + 1./sqr(-that-shat-mu2_dyn[1])));
    */
    that = -R*shat;
    meover = (1./36.)*(32*shat*shat*(1./sqr(-mu2_dyn[1]) + 1./sqr(-shat-mu2_dyn[1])));
    mefull = m_mes(u,u,u,u,shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": gg_qiqbi ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, qiqi_qiqi, selecting that, ntrial = " << ntrial << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  return that;
}

double Elastic_Kernel::SelectThat_qiqj_qiqj(const double shat) {
  int ntrial(0);
  Flavour u(kf_u), d(kf_d);
  double that,R,weight,mefull,meover;
  do {
    R=ran->Get();
    /*
    that = mu2_dyn[1]*shat*(R-1.)/(mu2_dyn[1]+R*shat);
    if (that < -shat || that > 0.)
      msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = 8.*shat*shat/(9.*sqr(that-mu2_dyn[1]));
    meover = (1./36.)*32.*shat*shat/sqr(that-mu2_dyn[1]);
    */
    that = -R*shat;
    meover = (1./36.)*32.*shat*shat/sqr(-mu2_dyn[1]);
    mefull = m_mes(u,d,u,d,shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": qiqj_qiqj ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, qiqj_qiqj, selecting that, ntrial = " << ntrial << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  return that;
}

double Elastic_Kernel::SelectThat_qiqbi_qiqbi(const double shat) {
  int ntrial(0);
  Flavour u(kf_u);
  double that,R,weight,mefull,meover,int1,int2;
  int1 = 8.*pow(shat,3)/(9.*mu2_dyn[1]*(shat+mu2_dyn[1]));
  int2 = 8.*pow(shat,3)*(4.*mu2_dyn[1]+shat)/(27.*mu2_dyn[1]*sqr(shat+mu2_dyn[1]));
  do {
    R=ran->Get();
    /*
    if (ran->Get() < int1/(int1+int2)) {
      that = mu2_dyn[1]*shat*(R-1.)/(mu2_dyn[1]+R*shat);
    }
    else {
      that = -R*shat;
    }
    if (that < -shat || that > 0.)
      msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = 8.*shat*shat/(9.*sqr(that-mu2_dyn[1])) + 8.*shat*shat*(4.*mu2_dyn[1]+shat)/(27.*mu2_dyn[1]*sqr(shat+mu2_dyn[1]));
    meover = (1./36.)*(32.*shat*shat*( 1./sqr(that-mu2_dyn[1]) + (4.*mu2_dyn[1]+shat)/(3.*mu2_dyn[1]*sqr(shat+mu2_dyn[1])) ));
    */
    that = -R*shat;
    meover = (1./36.)*(32.*shat*shat*( 1./sqr(-mu2_dyn[1]) + (4.*mu2_dyn[1]+shat)/(3.*mu2_dyn[1]*sqr(shat+mu2_dyn[1])) ));
    mefull = m_mes(u,u.Bar(),u,u.Bar(),shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": qiqbi_qiqbi ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, qiqbi_qiqbi, selecting that, ntrial = " << ntrial << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  return that;
}

double Elastic_Kernel::SelectThat_qiqbi_qjqbj(const double shat) {
  int ntrial(0);
  Flavour u(kf_u), d(kf_d);
  double that,R,weight,mefull,meover;
  do {
    R=ran->Get();
    that = -R*shat;
    if (that < -shat || that > 0.)
      msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = 8.*shat*shat/(9.*sqr(shat+mu2_dyn[1]));
    meover = (1./36.)*32.*shat*shat/sqr(shat+mu2_dyn[1]);
    mefull = m_mes(u,u.Bar(),d,d.Bar(),shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": qiqbi_qjqbj ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, qiqbi_qjqbj, selecting that, ntrial = " << ntrial << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  return that;
}

double Elastic_Kernel::SelectThat_qiqbi_gg(const double shat) {
  int ntrial(0);
  Flavour u(kf_u), g(kf_gluon);
  double that,R,weight,mefull,meover;
  do {
    R=ran->Get();
    /*
    //that = mu2_dyn[1] + (2*mu2_dyn[1]+shat)/(pow((shat+mu2_dyn[1])/mu2_dyn[1],2.*R-1.)-1.);
    that = mu2_dyn[1] - (2*mu2_dyn[1]+shat)/(pow((shat+mu2_dyn[1])/mu2_dyn[1],2.*R-1.)+1.); // m_f^2or m_g^2 first term here?
    if (that < -shat || that > 0.)
      msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = 32.*shat*shat/(27*(-that-shat-mu2_dyn[0])*(that-mu2_dyn[0]));
    meover = (1./36.)*256.*shat*shat/(3.*(-that-shat-mu2_dyn[0])*(that-mu2_dyn[0]));
    */
    that = -R*shat;
    meover = (1./36.)*256.*shat*shat/(3.*(-shat-mu2_dyn[0])*(-mu2_dyn[0]));
    mefull = m_mes(u,u.Bar(),g,g,shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": qiqbi_gg ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, qiqbi_gg, selecting that, ntrial = " << ntrial << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  return that;
}

double Elastic_Kernel::SelectThat_qig_qig(const double shat) {
  int ntrial(0);
  Flavour u(kf_u), g(kf_gluon);
  double that,R,weight,mefull,meover;
  do {
    R=ran->Get();
    /*
    if (ran->Get() < 9./13.) {
      that = mu2_dyn[1]*shat*(R-1.)/(mu2_dyn[1]+R*shat);
    }
    else {
      that = -R*shat;
    }
    if (that < -shat || that > 0.)
      msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = 2.*shat*shat/sqr(that-mu2_dyn[1]) + 8.*shat*shat/(9.*mu2_dyn[0]*(shat+mu2_dyn[0]));
    meover = (1./96.)*(256.*shat*shat/(3.*mu2_dyn[0]*(shat+mu2_dyn[0])) + 192.*shat*shat/sqr(that-mu2_dyn[1]));
    */
    that = -R*shat;
    meover = (1./96.)*(256.*shat*shat/(3.*mu2_dyn[0]*(shat+mu2_dyn[0])) + 192.*shat*shat/sqr(mu2_dyn[1]));
    mefull = m_mes(u,g,u,g,shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": qig_qig ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, qig_qig, selecting that, ntrial = " << ntrial << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  //msg_Out() << "that = " << that << ", mu2_dyn.size() = " << mu2_dyn.size() << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << "\n";
  return that;
}

double Elastic_Kernel::SelectThat_gg_gg(const double shat) {
  int ntrial(0);
  Flavour g(kf_gluon);
  double that,R,weight,mefull,meover;
  do {
    R=ran->Get();
    that = -R*shat;
    if (that < -shat || that > 0.) msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = 9./2.*(3.+sqr(shat/mu2_dyn[1]));
    //meover = (1./256.)*1152.*(3.+sqr(shat/(that-mu2_dyn[1])) - shat*that/sqr(-shat-that-mu2_dyn[1]));
    meover = (1./256.)*1152.*(3. + shat*shat*2./(sqr(mu2_dyn[1])) + shat*shat/(sqr(shat+ mu2_dyn[1])));
    mefull = m_mes(g,g,g,g,shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": gg_gg ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, gg->gg, selecting that, ntrial = " << ntrial << ", shat = " << shat << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  return that;
}

double Elastic_Kernel::SelectThat_gg_qiqbi(const double shat) {
  int ntrial(0);
  Flavour u(kf_u), g(kf_gluon);
  double that,R,weight,mefull,meover;
  do {
    R=ran->Get();
    /*
//     that = mu2_dyn[1] + (2*mu2_dyn[1]+shat)/(pow((shat+mu2_dyn[1])/mu2_dyn[1],2.*R-1.)-1.);
    that = mu2_dyn[1] - (2*mu2_dyn[1]+shat)/(pow((shat+mu2_dyn[1])/mu2_dyn[1],2.*R-1.)+1.);
    if (that < -shat || that > 0.)
      msg_Error()<<METHOD<<": selected that out of bounds: that = "<<that<<" for shat = "<<shat<<".\n";
    //meover = shat*shat/(6.*(-that-shat-mu2_dyn[1])*(that-mu2_dyn[1]));
    meover = (1./256.)*256.*shat*shat/(3.*(-that-shat-mu2_dyn[0])*(that-mu2_dyn[0]));
    */
    that = -R*shat;
    meover = (1./256.)*256.*shat*shat/(3.*(-shat-mu2_dyn[0])*(-mu2_dyn[0]));
    mefull = m_mes(g,g,u,u.Bar(),shat,that,-shat-that);
    if(meover < mefull){
      msg_Out() << METHOD << ": gg_qiqbi ME weight > 1! meover = " << meover << ", mefull = " << mefull << endl;
      msg_Out() << METHOD << ": shat = " << shat << " that = " << that << ", mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
    }
    //weight = mefull/meover*sqr(p_dyn_quant_handler->GetRunningAlphaS(that+mu2_dyn[1])/p_dyn_quant_handler->GetRunningAlphaS(mu2_dyn[1]));
    weight = mefull/meover;
    ntrial++;
    if(ntrial>1e7){
      msg_Out() << METHOD << ": WARNING, gg_qiqbi, selecting that, ntrial = " << ntrial << endl;
      return 1.;
    }
  } while (ran->Get() > weight);
  return that;
}

void Elastic_Kernel::SetMu2Dyn(const std::vector<double> mu2_temp) { 
  mu2_dyn = {(exp(5./3.)/4.)*mu2_temp[0], (exp(5./3.)/4.)*mu2_temp[1]} ; //Multipling by experimental constant
  m_mes.SetMu2Dyn({(exp(5./3.)/4.)*mu2_temp[0], (exp(5./3.)/4.)*mu2_temp[1]});
}
