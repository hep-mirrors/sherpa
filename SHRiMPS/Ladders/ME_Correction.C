#include "SHRiMPS/Ladders/ME_Correction.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

ME_Correction::ME_Correction(const double & tmin) :
  m_tmin(tmin), m_Ecms(rpa->gen.Ecms()), m_S(sqr(m_Ecms))
{
  for (size_t i=0;i<2;i++) {
    p_pdf[i]   = NULL;
    m_Ebeam[i] = rpa->gen.PBeam(i)[0];
  }
}

double ME_Correction::operator()(Ladder * ladder) {
  if (ladder->GetEmissions()->size()==2) return 1.;
  double y1, y2, s_12, t_12, ME_12 = 0.;
  ladder->ExtractHardTwo2Two(y1,y2,s_12,t_12);
  double y12 = (y1+y2)/2., cosh_dy12 = cosh((y1-y2)/2.);
  double kt2 = Max(s_12/sqr(2.*cosh_dy12),m_tmin), kt = sqrt(kt2);
  double x1  = Max(1.001e-5,2.*kt/m_Ecms * exp( y12) * cosh_dy12);
  double x2  = Max(1.001e-5,2.*kt/m_Ecms * exp(-y12) * cosh_dy12);
  if (x1<0.99 && x1>1.e-5 && x2<0.99 && x2>1.e-5) {
    CalcPDFs(x1,x2,kt2);
    ME_12 = ( 1./(2.*s_12) *
	      p_pdf[0]->XPDF(Flavour(kf_gluon)) *
	      p_pdf[1]->XPDF(Flavour(kf_gluon)) *
	      s_12/Max(t_12, m_tmin) );
  }
  double ya, yb, s_ab;
  ladder->ExtractExternalTwo2Two(ya,yb,s_ab);
  double yab   = (ya+yb)/2., cosh_dyab = cosh((ya-yb)/2.);
  double qt2   = s_ab/sqr(2.*cosh_dyab), qt = sqrt(qt2);
  double xa    = 2.*qt/m_Ecms * exp( yab) * cosh_dyab;
  double xb    = 2.*qt/m_Ecms * exp(-yab) * cosh_dyab;
  if (xa<0.99 && xa>1.e-5 && xb<0.99 && xb>1.e-5) {
    CalcPDFs(xa,xb,0.);
    double ME_ab = ( 1./(2.*s_ab) *
		     p_pdf[0]->XPDF(Flavour(kf_gluon)) *
		     p_pdf[1]->XPDF(Flavour(kf_gluon)) *
		     pow(s_ab/m_tmin,0.08) );
    //msg_Out()<<"* |ME^2("<<y1<<", "<<y2<<", "<<s_12<<", "<<t_12<<")| = "<<ME_12<<" "
    //	   <<"for x1 = "<<x1<<", x2 = "<<x2<<", mu = "<<sqrt(kt2)<<"\n"
    //	   <<"  |ME^2("<<ya<<", "<<yb<<", "<<s_ab<<", "<<m_tmin<<")| = "<<ME_ab<<" "
    //	   <<"for xa = "<<xa<<", xb = "<<xb<<", mu = "<<sqrt(qt2)<<"\n"
    //	   <<"*** ratio = "<<ME_12/ME_ab<<"\n";
    return ME_12/ME_ab;
  }
  //msg_Out()<<"* |ME^2("<<y1<<", "<<y2<<", "<<s_12<<", "<<t_12<<")| = "<<ME_12<<" "
  //	   <<"for x1 = "<<x1<<", x2 = "<<x2<<", mu = "<<sqrt(kt2)<<"\n"
  //	   <<"  |ME^2("<<ya<<", "<<yb<<", "<<s_ab<<", "<<m_tmin<<")| = "<<0.<<" "
  //	   <<"for xa = "<<xa<<", xb = "<<xb<<", mu = "<<sqrt(qt2)<<"\n";
  return 0.;
}

void ME_Correction::CalcPDFs(const double & x1,const double & x2,
			    const double & scale) {
  // Calculate both sets of PDFs at the relevant x and Q^2
  p_pdf[0]->Calculate(x1,Max(scale,p_pdf[0]->Q2Min()));
  p_pdf[1]->Calculate(x2,Max(scale,p_pdf[1]->Q2Min()));
}

