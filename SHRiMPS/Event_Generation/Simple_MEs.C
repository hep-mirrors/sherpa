#include "SHRiMPS/Event_Generation/Simple_MEs.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;

Simple_MEs::Simple_MEs() {}
Simple_MEs::~Simple_MEs() {}

double Simple_MEs::
operator()(const ATOOLS::Flavour & in1,const ATOOLS::Flavour & in2,
	   const ATOOLS::Flavour & out1,const ATOOLS::Flavour & out2,
	   const double & hats,const double & hatt,const double & hatu,
	   const double & mu2) {
  double me2(0.);
  if (in1.IsQuark() && in2.IsQuark()) {
    if (in1==in2)            me2 = qq_qq(hats,hatt,hatu,mu2);
    else if (in1==in2.Bar()) me2 = qqb_qqb(hats,hatt,hatu,mu2);
    else if (in1!=in2)       me2 = q1q2_q1q2(hats,hatt,hatu,mu2);
  }
  else if (in1.IsGluon() && in2.IsQuark())
    me2 = gq_gq(hats,hatt,hatu,mu2);
  else if (in1.IsQuark() && in2.IsGluon())
    me2 = gq_gq(hats,hatt,hatu,mu2);
  else if ((in1.IsGluon() && in2.IsGluon()) && (out1.IsGluon()))
    me2 = gg_gg(hats,hatt,hatu,mu2);
  else if ((in1.IsGluon() && in2.IsGluon()) && 
	   (out1.IsQuark() || out2.IsQuark()))
    me2 = gg_qqb(hats,hatt,hatu,mu2);

  if (out1==out2) me2 /= 2.;
  // factor of 16 Pi^2 to go from g_S^4 to alpha_S^2
  return ATOOLS::sqr(4.*M_PI)*me2;
}

double Simple_MEs::
qq_qq(const double & hats,const double & hatt,const double & hatu,
      const double & mu2) {
  return 
    ATOOLS::sqr(3./4.) *
    (4./9.  * ((hats*hats+hatu*hatu)/((hatt-mu2)*(hatt-mu2))+
	       (hats*hats+hatt*hatt)/((hatu-mu2)*(hatu-mu2))) -
     8./27. * (hats*hats)/((hatu-mu2)*(hatt-mu2)));
}

double Simple_MEs::
qqb_qqb(const double & hats,const double & hatt,const double & hatu,
	const double & mu2) {
  return 
    //ATOOLS::sqr(3./4.) *
    (4./9.  * ((hats*hats+hatu*hatu)/((hatt-mu2)*(hatt-mu2))+
	       (hatt*hatt+hatu*hatu)/((hats+mu2)*(hats+mu2))) -
     8./27. * (hatu*hatu)/((hats+mu2)*(hatt-mu2)));
}

double Simple_MEs::
q1q2_q1q2(const double & hats,const double & hatt,const double & hatu,
	  const double & mu2) {
  return 
    //ATOOLS::sqr(3./4.) *
    4./9. * (hats*hats+hatu*hatu)/((hatt-mu2)*(hatt-mu2));
}

double Simple_MEs::
gq_gq(const double & hats,const double & hatt,const double & hatu,
      const double & mu2) {
  return
    //3./4.*1./3. *
    (-4./9. * (hats*hats+hatu*hatu)/((hats+mu2)*(hatu-mu2)) +
     (hats*hats+hatu*hatu)/((hatt-mu2)*(hatt-mu2)));
}

double Simple_MEs::
gg_gg(const double & hats,const double & hatt,const double & hatu,
      const double & mu2) {
  return 
    //ATOOLS::sqr(1./3.) * 
    (9./2.*(1.-(hatt*hatu)/((hats+mu2)*(hats+mu2)) +
	    1.-(hats*hatu)/((hatt-mu2)*(hatt-mu2)) +
	    1.-(hats*hatt)/((hatu-mu2)*(hatu-mu2))));
}

double Simple_MEs::
gg_qqb(const double & hats,const double & hatt,const double & hatu,
       const double & mu2) {
  return 
    1./6.*(hatt*hatt+hatu*hatu)/((hatt-mu2)*(hatu-mu2))- 
    1./6.*(hatt*hatt+hatu*hatu)/((hats+mu2)*(hats+mu2));
}
