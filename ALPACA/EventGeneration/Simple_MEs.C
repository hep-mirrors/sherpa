#include "ALPACA/EventGeneration/Simple_MEs.H"

#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace std;
using namespace ALPACA;
using namespace ATOOLS;

Simple_MEs::Simple_MEs():
  mu2_dyn{0.001,0.01}
{ }
Simple_MEs::~Simple_MEs() {}

double Simple_MEs::
operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
	   const ATOOLS::Flavour & out1,const ATOOLS::Flavour & out2,
	   const double & hats,const double & hatt,const double & hatu) {
  double me2(0.);
  if (in1.IsQuark() && in2.IsQuark()){
    if (in1==in2){
      me2 = qiqi_qiqi(hats,hatt,hatu);
    } else if (in1==in2.Bar()){
      if (out1.IsQuark()){
        if (in1==out1 || in1==out2){
          me2 = qiqbi_qiqbi(hats,hatt,hatu);
        } else {
          me2 = qiqbi_qjqbj(hats,hatt,hatu);
        }
      } else{
	      me2 = qqb_gg(hats,hatt,hatu);
      }
    } else if (in1!=in2){
      me2 = qiqj_qiqj(hats,hatt,hatu);
    }
  }
  else if (in1.IsGluon() && in2.IsQuark())
    me2 = gq_gq(hats,hatt,hatu);
  else if (in1.IsQuark() && in2.IsGluon())
    me2 = gq_gq(hats,hatt,hatu);
  else if ((in1.IsGluon() && in2.IsGluon()) && (out1.IsGluon()))
    me2 = gg_gg(hats,hatt,hatu);
  else if ((in1.IsGluon() && in2.IsGluon()) && 
	   (out1.IsQuark() || out2.IsQuark()))
    me2 = gg_qqb(hats,hatt,hatu);
  else THROW(fatal_error,"Did not find process");

//   if (out1==out2) me2 /= 2.;
//   msg_Out()<<METHOD<<" returns "<<me2<<"\n";
  return me2;
}

//Analytical expressions of matrix elements, depending on mu2
double Simple_MEs::
operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
	   const ATOOLS::Flavour & out1,const ATOOLS::Flavour & out2,
	   const double & hats) {
  double me2 = 0.;
  double hatt_min = -hats;
  double hatt_max = 0.;
  if (in1.IsQuark() && in2.IsQuark()){
    if (in1==in2){ 
      me2 = qiqi_qiqi_mu2(hats,hatt_max)-qiqi_qiqi_mu2(hats,hatt_min);
      if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
        me2 = 0.;
        //msg_Out() << METHOD << " ERROR, qiqi_qiqi_mu2, me2 = " << me2 << endl;
      }
    } else if (in1==in2.Bar()){
      if (out1.IsQuark()) {
        if (in1==out1 || in1==out2){
          me2 = qiqbi_qiqbi_mu2(hats,hatt_max)-qiqbi_qiqbi_mu2(hats,hatt_min);
          if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
            me2 = 0.;
            //msg_Out() << METHOD << " ERROR, qiqbi_qiqbi_mu2, me2 = " << me2 << endl;
          }
        } else{ 
          me2 = qiqbi_qjqbj_mu2(hats,hatt_max)-qiqbi_qjqbj_mu2(hats,hatt_min);
          if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
            me2 = 0.;
            //msg_Out() << METHOD << " ERROR, qiqbi_qjqbj_mu2, me2 = " << me2 << endl;
          }
        }
      } else{ 
	      me2 = qqb_gg_mu2(hats,hatt_max)-qqb_gg_mu2(hats,hatt_min);
        //me2 = 0.;
        if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
          //msg_Out() << "\n\n ######## ERROR ######## \n\n" << endl;
          //msg_Out() << METHOD << " ERROR, qqb -> gg, me2 = " << me2 << endl;
          //msg_Out() << METHOD << " hats = " << hats << ", hatt_max = " << hatt_max << ", hatt_min = " << hatt_min << endl;
          //msg_Out() << METHOD << " mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
          me2 = 0.;
        }
      }
    } else if (in1!=in2){
      me2 = qiqj_qiqj_mu2(hats,hatt_max)-qiqj_qiqj_mu2(hats,hatt_min);
      if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
        me2 = 0.;
        //msg_Out() << METHOD << " ERROR, qiqj_qiqj_mu2, me2 = " << me2 << endl;
      }
    }
  } else if (in1.IsGluon() && in2.IsQuark()){
    me2 = gq_gq_mu2(hats,hatt_max)-gq_gq_mu2(hats,hatt_min);
    //me2 = 0.;
    if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
      me2 = 0.;
      //msg_Out() << METHOD << " ERROR, gq_gq_mu2, me2 = " << me2 << endl;
    }
  } else if (in1.IsQuark() && in2.IsGluon()){
    me2 = gq_gq_mu2(hats,hatt_max)-gq_gq_mu2(hats,hatt_min);
    //me2 = 0.;
    if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
      me2 = 0.;
      //msg_Out() << METHOD << " ERROR, gq_gq_mu2, me2 = " << me2 << endl;
    }
  } else if ((in1.IsGluon() && in2.IsGluon()) && (out1.IsGluon())){
    me2 = gg_gg_mu2(hats,hatt_max)-gg_gg_mu2(hats,hatt_min);
    //me2 = 0.;
    if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
      me2 = 0.;
      //msg_Out() << METHOD << " ERROR, gg_gg_mu2, me2 = " << me2 << endl;
    }
  } else if ((in1.IsGluon() && in2.IsGluon()) && 
	   (out1.IsQuark() || out2.IsQuark())){
    me2 = gg_qqb_mu2(hats,hatt_max)-gg_qqb_mu2(hats,hatt_min);
    //me2 = 0.;
    if(me2 < 0 || std::isnan(me2) || std::isinf(me2)){
      //msg_Out() << "\n\n ######## ERROR ######## \n\n" << endl;
      //msg_Out() << METHOD << " ERROR, gg -> qqb, me2 = " << me2 << endl;
      //msg_Out() << METHOD << " hats = " << hats << ", hatt_max = " << hatt_max << ", hatt_min = " << hatt_min << endl;
      //msg_Out() << METHOD << " mu2_dyn[0] = " << mu2_dyn[0] << ", mu2_dyn[1] = " << mu2_dyn[1] << endl;
      me2 = 0.;
    }
  } else {THROW(fatal_error,"Did not find process")};

  return me2;
}

double Simple_MEs::
operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
	   const double & hats,const double & hatt,const double & hatu) {
  double me2(0.);
  if (in1.IsQuark() && in2.IsQuark()) {
    if (in1==in2) 
      me2 = qiqi_qiqi(hats,hatt,hatu);
    else if (in1==in2.Bar()) {
      me2 =  qiqbi_qiqbi(hats,hatt,hatu);
      me2 += qiqbi_qjqbj(hats,hatt,hatu);
      me2 += qqb_gg(hats,hatt,hatu);
    }
    else if (in1!=in2)
      me2 = qiqj_qiqj(hats,hatt,hatu);
  }
  else if (in1.IsGluon() && in2.IsQuark())
    me2 = gq_gq(hats,hatt,hatu);
  else if (in1.IsQuark() && in2.IsGluon())
    me2 = gq_gq(hats,hatt,hatu);
  else if (in1.IsGluon() && in2.IsGluon()) {
    me2 = gg_gg(hats,hatt,hatu);
    me2 += gg_qqb(hats,hatt,hatu);
  }
  else THROW(fatal_error,"Did not find process");

//   if (out1==out2) me2 /= 2.;
  return ATOOLS::Max(0.,me2);
}



double Simple_MEs::
qiqi_qiqi(const double & hats,const double & hatt,const double & hatu) {
    return (1./36.)*(16. * ((hats*hats+hatu*hatu)/((hatt-mu2_dyn[1])*(hatt-mu2_dyn[1]))+
      (hats*hats+hatt*hatt)/((hatu-mu2_dyn[1])*(hatu-mu2_dyn[1]))) -
      32./3. * (hats*hats)/((hatu-mu2_dyn[1])*(hatt-mu2_dyn[1]))  );
}

double Simple_MEs::
qiqbi_qiqbi(const double & hats,const double & hatt,const double & hatu) {
  return (1./36.)*(16. * ((hats*hats+hatu*hatu)/((hatt-mu2_dyn[1])*(hatt-mu2_dyn[1]))+
      (hatt*hatt+hatu*hatu)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1]))) -
      32./3. * (hatu*hatu)/((hats+mu2_dyn[1])*(hatt-mu2_dyn[1])) );
}

double Simple_MEs::
qiqbi_qjqbj(const double & hats,const double & hatt,const double & hatu) {
  return (1./36.)*16.*(hatt*hatt+hatu*hatu)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1]));
}

double Simple_MEs::
qqb_gg(const double & hats,const double & hatt,const double & hatu) {
  return (1./36.)*(128./3. * (hatu*hatu+hatt*hatt)/((hatu-mu2_dyn[0])*(hatt-mu2_dyn[0]))-
      96. * (hatu*hatu+hatt*hatt)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1])) );
}

double Simple_MEs::
gg_qqb(const double & hats,const double & hatt,const double & hatu) {
  return (1./256.)*(128./3. * (hatu*hatu+hatt*hatt)/((hatu-mu2_dyn[0])*(hatt-mu2_dyn[0]))-
      96. * (hatu*hatu+hatt*hatt)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1])) );
}

double Simple_MEs::
qiqj_qiqj(const double & hats,const double & hatt,const double & hatu) {
  return (1./36.)*16.*(hats*hats+hatu*hatu)/((hatt-mu2_dyn[1])*(hatt-mu2_dyn[1]));
}

double Simple_MEs::
gq_gq(const double & hats,const double & hatt,const double & hatu) {
  return (1./96.)*(-(128./3.)*(hats*hats+hatu*hatu)/((hats+mu2_dyn[0])*(hatu-mu2_dyn[0])) +
    96. * (hats*hats+hatu*hatu)/((hatt-mu2_dyn[1])*(hatt-mu2_dyn[1])));
}

double Simple_MEs::
gg_gg(const double & hats,const double & hatt,const double & hatu) {
  return (1./256.)*(1152*(1.-(hatt*hatu)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1])) +
                          1.-(hats*hatu)/((hatt-mu2_dyn[1])*(hatt-mu2_dyn[1])) +
                          1.-(hats*hatt)/((hatu-mu2_dyn[1])*(hatu-mu2_dyn[1])) ));
}

/*
double Simple_MEs::
gg_qqb(const double & hats,const double & hatt,const double & hatu) {
return 
// ATOOLS::sqr(hats/(hats+mu2_dyn[1]))*
1./6.*(hatt*hatt+hatu*hatu)/((hatt-mu2_dyn[1])*(hatu-mu2_dyn[1])) - 3./8.*(hatt*hatt+hatu*hatu)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1]));
}
*/


//Integrated matrix elements, depending on mu2_dyn[0] (viritual fermion) and m2_dyn[1] (viritual boson).
double Simple_MEs::
qiqi_qiqi_mu2(const double & hats, const double & hatt) {

  return (1./36.)*(16./3.)*(
                    6.*(hats+mu2_dyn[1])*log(abs((hatt-mu2_dyn[1])/(hats+hatt+mu2_dyn[1])))
                    +2.*hats*hats*log(abs((hatt-mu2_dyn[1])/(hats+hatt+mu2_dyn[1])))/(hats+2.*mu2_dyn[1])
                    +6.*hatt
                    +3.*(mu2_dyn[1]*mu2_dyn[1]+2.*mu2_dyn[1]*hats+2.*hats*hats)*(1./(mu2_dyn[1]-hatt) - 1./(hats+hatt+mu2_dyn[1]))
                  );
  /* old version
  return (1./36.)*(16./3.)*(
                    6.*(hats+mu2_dyn[1])*log(abs((hatt-mu2_dyn[1])/(hats+hatt+mu2_dyn[1])))
                    +2.*hats*hats*log(abs((hatt-mu2_dyn[1])/(hats+hatt+mu2_dyn[1])))/(hats+mu2_dyn[1])
                    +6.*(hatt-mu2_dyn[1])
                    +3.*(mu2_dyn[1]*mu2_dyn[1]+2.*mu2_dyn[1]*hats+2.*hats*hats)*(1./(mu2_dyn[1]-hatt) - 1./(hats+hatt+mu2_dyn[1]))
                  );
  */
}

double Simple_MEs::
qiqbi_qiqbi_mu2(const double & hats, const double & hatt) {
  return  (1./36.)*16.*(
                -2.*pow(mu2_dyn[1],4)
                +3.*pow(mu2_dyn[1],3)*(hats+2.*hatt)
                +sqr(mu2_dyn[1])*(15.*hats*hats + 9.*hats*hatt -2.*hatt*hatt)
                +mu2_dyn[1]*(18.*hats*hats*hats + 8*hats*hats*hatt +2.*hats*hatt*hatt + 3.*hatt*hatt*hatt)
                +4.*pow(hats+mu2_dyn[1],3)*(mu2_dyn[1]-hatt)*log(abs(mu2_dyn[1]-hatt))
                -2.*(-3*hats*hats*hats*hats + hats*hats*hatt*hatt + hats*hatt*hatt*hatt + hatt*hatt*hatt*hatt)
              )/(3.*sqr(hats+mu2_dyn[1])*(mu2_dyn[1]-hatt));
}

double Simple_MEs::
qiqbi_qjqbj_mu2(const double & hats, const double & hatt) {
  return (1./36.)*16.*hatt*(
                    3.*hats*hats + 3.*hats*hatt + 2.*hatt*hatt
                  )/(3.*(hats+mu2_dyn[1])*(hats+mu2_dyn[1]));
}

double Simple_MEs::
qqb_gg_mu2(const double & hats, const double & hatt) {
  /*
  //Explicitly solved with x = -(t-mf2)
  double mf2 = mu2_dyn[0];
  double mg2 = mu2_dyn[1];
  double x = -(hatt-mf2);
  double s = hats;
  double B = (96./3.)*(pow(s+hatt,3.)+pow(hatt,3.))/pow(mg2+s,2.);
  double A = -(128./3.)*(1./(2.*mf2+s))*( 2.*(2.*mf2+s)*(mf2-x) + (2.*mf2*mf2 + 2.*mf2*s + s*s)*log(x/(2.*mf2+s-x)) );
  //return (1./36.)*(A-B);
  */

  return (1./36.)*(32./3.)*(
                    4.*(pow(hats+mu2_dyn[0],2) + mu2_dyn[0]*mu2_dyn[0])*log(abs((hats+hatt+mu2_dyn[0])/(hatt-mu2_dyn[0])))/(2.*mu2_dyn[0]+hats)
                    //-4.*(pow(hats+mu2_dyn[0],2) + mu2_dyn[0]*mu2_dyn[0])*log(abs((hatt-mu2_dyn[0])/(hats+hatt+mu2_dyn[0])))/(2.*mu2_dyn[0]+hats)
                    -hatt*(8.*pow(hats+mu2_dyn[1],2) + 9.*hats*hats + 9.*hatt*hats + 6.*hatt*hatt)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1]))
                  );
  
}

double Simple_MEs::
gg_qqb_mu2(const double & hats, const double & hatt) {  
  
  /*
  //Explicitly solved with x = -(t-mf2)
  double mf2 = mu2_dyn[0];
  double mg2 = mu2_dyn[1];
  double x = -(hatt-mf2);
  double s = hats;
  double B = (96./3.)*(pow(s+hatt,3.)+pow(hatt,3.))/pow(mg2+s,2.);
  double A = -(128./3.)*(1./(2.*mf2+s))*( 2.*(2.*mf2+s)*(mf2-x) + (2.*mf2*mf2 + 2.*mf2*s + s*s)*log(x/(2.*mf2+s-x)) );
  //return (1./36.)*(A-B);
  */
  return (1./256.)*(32./3.)*(
                    //-4.*(pow(hats+mu2_dyn[0],2) + mu2_dyn[0]*mu2_dyn[0])*log(abs((hatt-mu2_dyn[0])/(hats+hatt+mu2_dyn[0])))/(2.*mu2_dyn[0]+hats)
                    4.*(pow(hats+mu2_dyn[0],2) + mu2_dyn[0]*mu2_dyn[0])*log(abs((hats+hatt+mu2_dyn[0])/(hatt-mu2_dyn[0])))/(2.*mu2_dyn[0]+hats)
                    -hatt*(8.*pow(hats+mu2_dyn[1],2) + 9.*hats*hats + 9.*hatt*hats + 6.*hatt*hatt)/((hats+mu2_dyn[1])*(hats+mu2_dyn[1]))
                  );
}

double Simple_MEs::
qiqj_qiqj_mu2(const double & hats, const double & hatt) {
  return (1./36.)*16.*(
                (mu2_dyn[1]*mu2_dyn[1] + 2.*mu2_dyn[1]*hats + 2.*hats*hats)/(mu2_dyn[1]-hatt)
                +2.*(hats+mu2_dyn[1])*log(abs(mu2_dyn[1]-hatt))
                +hatt
             );
}

double Simple_MEs::
gq_gq_mu2(const double & hats, const double & hatt) {
  double A, B, C, D, retval;
  A = (1./96.)*(64./3.)*(hatt*(hatt+2.*hats-2.*mu2_dyn[0]) +  2.*(mu2_dyn[0]*mu2_dyn[0]+hats*hats)*log(hats+hatt+mu2_dyn[0]))/(mu2_dyn[0]+hats);

  B = (1./96.)*96.*(2.*(mu2_dyn[1]+hats)*log(mu2_dyn[1]-hatt) + (mu2_dyn[1]*mu2_dyn[1] + 2.*mu2_dyn[1]*hats +2.*hats*hats)/(mu2_dyn[1]-hatt) + hatt );

  C = (1./96.)*(32./3.)*(
                    9.*(mu2_dyn[1]*mu2_dyn[1] + 2.*mu2_dyn[1]*hats + 2.*hats*hats)/(mu2_dyn[1]-hatt)
                    +18.*(hats+mu2_dyn[1])*log(abs(hatt-mu2_dyn[1]))
                    +(hatt*(5.*mu2_dyn[0]+13.*hats) +2.*hatt*hatt + 4.*(hats*hats+mu2_dyn[0]*mu2_dyn[0])*log(abs(hats+hatt+mu2_dyn[0])))/(hats+mu2_dyn[0])
                  );

  /*
  D = (1./96.)*(32./3.)*(
                    9.*(mu2_dyn[1]*mu2_dyn[1] + 2.*mu2_dyn[1]*hats + 2.*hats*hats)/(mu2_dyn[1]-hatt)
                    +18.*(hats+mu2_dyn[1])*log(abs(hatt-mu2_dyn[1]))
                    +(hatt*(13.*mu2_dyn[0]+5.*hats) +2.*hatt*hatt + 4.*(hats*hats+mu2_dyn[0]*mu2_dyn[0])*log(abs(hats+hatt+mu2_dyn[0])))/(hats+mu2_dyn[0])
                  ); //old version

  
  //msg_Out() << " A+B = " << A+B << ", C = " << C << ", D = " << D << endl;
  */
  retval = A + B;
  return retval;
  /*
  
  return (1./96.)*(32./3.)*(
                    9.*(mu2_dyn[1]*mu2_dyn[1] + 2.*mu2_dyn[1]*hats + 2.*hats*hats)/(mu2_dyn[1]-hatt)
                    +18.*(hats+mu2_dyn[1])*log(abs(hatt-mu2_dyn[1]))
                    +(hatt*(13.*mu2_dyn[0]+5.*hats) +2.*hatt*hatt + 4.*(hats*hats+mu2_dyn[0]*mu2_dyn[0])*log(abs(hats+hatt+mu2_dyn[0])))/(hats+mu2_dyn[0])
                  );
  */
}

double Simple_MEs::
gg_gg_mu2(const double & hats, const double & hatt) {
  return (1./256.)*1152.*(
                  hats*log(abs((hatt-mu2_dyn[1])/(mu2_dyn[1]+hats+hatt)))
                  +3.*hatt
                  +(2.*hatt*hatt*hatt + 3.*hats*hatt*hatt)/(6.*pow(hats+mu2_dyn[1],2))
                  +hats*(mu2_dyn[1]+hats)*(1./(mu2_dyn[1]-hatt) - 1./(mu2_dyn[1]+hats+hatt))
               );
}

/*
double Simple_MEs::
gg_qqb_mu2(const double & hats, const double & hatt) {
  return ;
}
*/