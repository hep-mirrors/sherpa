#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "SHRiMPS/Cross_Sections/Sigma_Total.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;
using namespace ATOOLS;


Cross_Sections::Cross_Sections() :
  p_selected(NULL),
  m_xstot(0.), m_slope(0.), m_xsinel(0.), m_xsel(0.),
  m_sigma_inelastic(Sigma_Inelastic()),
  m_sigma_elastic(Sigma_Elastic()),
  m_sigma_D(Sigma_D())
{ }

Cross_Sections::~Cross_Sections()
{ }

void Cross_Sections::CalculateCrossSections()
{
  Sigma_Tot sigma_tot;
  m_xstot  = sigma_tot.Calculate();
  Elastic_Slope slope(m_xstot);
  m_slope  = slope.Calculate();
  m_xsinel = m_sigma_inelastic.Calculate();
  m_xsel   = m_sigma_elastic.Calculate();
  m_sigma_elastic.FillGrids();
  m_sigma_D.FillGrids(&m_sigma_elastic);
  //SaveXS();
  for (size_t i=0;i<2;i++) m_xsSD[i] = m_sigma_D.GetXSec(i);
  m_xsDD = m_sigma_D.GetXSec(2);
  double sigma_qe = m_sigma_D.Calculate();

  msg_Info()<<"===========================================================\n"
	    <<"   sigma_tot                 = "<<m_xstot/1.e9<<" mb, (B = "<<m_slope<<")\n"
	    <<"   sigma_inel                = "<<m_xsinel/1.e9<<" mb\n"   
	    <<"   sigma_el                  = "<<m_xsel/1.e9<<" mb\n"   
	    <<"      test: int dt dsigma/dt = "<<m_sigma_elastic.Summed()/1.e9<<" mb\n"
        <<"      ratio: = "<<m_xsel/m_sigma_elastic.Summed()<<"\n"
        <<"   sigma_SD0                 = "<<m_xsSD[0]/1.e9<<" mb\n"
	    <<"   sigma_SD1                 = "<<m_xsSD[1]/1.e9<<" mb\n"
        <<"   sigma_DD                  = "<<m_xsDD/1.e9<<" mb\n"
        <<"   sigma_QE                  = "<<sigma_qe/1.e9<<" mb\n"
        <<"      test: int dt dsigma/dt = "<<(m_sigma_elastic.Summed()+m_xsSD[0]+m_xsSD[1]+m_xsDD)/1.e9<<" mb\n"
        <<"      ratio: = "<<sigma_qe/(m_sigma_elastic.Summed()+m_xsSD[0]+m_xsSD[1]+m_xsDD)<<"\n."
        <<"===========================================================\n";
  MBpars.SetXSecs(this);
  m_sigma_inelastic.SetSigma(m_xsinel);
  m_sigma_elastic.SetSigma(m_xsel);
  m_sigma_D.SetSigma(m_xsSD[0]+m_xsSD[1]+m_xsDD);
  
}

void Cross_Sections::SaveXS() {
  std::ofstream XSfile, ELASfile, QEfile;
  XSfile.open("./cross_section_273.txt");
  double tmin = m_sigma_elastic.Tmin();
  double tmax = m_sigma_elastic.Tmax();
  size_t steps = m_sigma_elastic.Steps();
  double deltaT = (tmax - tmin) / steps;
  for (int diff = 0; diff < 3; diff++) {
    for (int i = 0; i < steps; i++) {
      double t = tmin + i*deltaT;
      XSfile << diff << "\t" << t << "\t" << m_sigma_D.GetXSvsT(diff,t) << std::endl;
    }
  }
  XSfile.close();
  ELASfile.open("./sig_elas_273.txt");
  for (int i = 0; i < steps; i++) {
    double t = tmin + i*deltaT;
    ELASfile << m_sigma_elastic.GetXSvsT(t) << std::endl;
  }
  ELASfile.close();
  QEfile.open("./sig_QE_273.txt");
  for (int i = 0; i < steps; i++) {
    double t = tmin + i*deltaT;
    double qevalue = 0.;
    for (int diff = 0; diff < 3; diff++) qevalue+=m_sigma_D.GetXSvsT(diff,t);
    qevalue += m_sigma_elastic.GetXSvsT(t);
    QEfile << qevalue << std::endl;
  }
  QEfile.close();
  std::ofstream XSfileB, ELASfileB, QEfileB;
  double bmin = 0.1;
  double bmax = 20.;//m_sigma_elastic.Bmax();
  double delta_b = 0.01;// m_sigma_elastic.deltaB();
  XSfileB.open("./crosssection_273_B.txt");
  ELASfileB.open("./sig_elas_273_B.txt");
  QEfileB.open("./sig_QE_273_B.txt");
  double qe, elastic, diffractives[3];
  for (double b = bmin; b <= bmax; b = b + delta_b) {
    qe = 0.;
    elastic = m_sigma_D.GetCombinedValueEL(b);
    diffractives[0] = m_sigma_D.GetCombinedValueSD0(b);
    diffractives[1] = m_sigma_D.GetCombinedValueSD1(b);
    diffractives[2] = m_sigma_D.GetCombinedValueDD(b);
    qe += elastic;
    for(int diff = 0; diff < 3; diff++) qe += diffractives[diff];
    ELASfileB << elastic << std::endl;
    QEfileB << qe << std::endl;
    for (int diff = 0; diff < 3; diff++){
      XSfileB << diff << "\t" << b << "\t" << diffractives[diff] << std::endl;
    }
  }
  XSfileB.close();
  ELASfileB.close();
  QEfileB.close();
}

void Cross_Sections::Test(const std::string & dirname)
{
  Sigma_Tot sigma_tot;
  Sigma_Inelastic sigma_inel;
  Sigma_Elastic sigma_el;
  double xstot_test(sigma_tot.Test());
  double xsinel_test(sigma_inel.Test());
  double xsel_test(sigma_el.Test());
  msg_Info()<<METHOD<<":\n"
	    <<"   sigma_tot  = "<<m_xstot/1.e9<<" mb "
	    <<"vs. "<<xstot_test/1.e9<<" mb,"
	    <<"   relative deviation = "
	    <<dabs(1.-m_xstot/xstot_test)*100.<<" %\n"   
	    <<"   sigma_el   = "<<m_xsel/1.e9<<" mb "
	    <<"vs. "<<xsel_test/1.e9<<" mb,"
	    <<"   relative deviation = "
	    <<dabs(1.-m_xsel/xsel_test)*100.<<" %\n"   
	    <<"   sigma_inel = "<<m_xsinel/1.e9<<" mb "
	    <<"vs. "<<xsinel_test/1.e9<<" mb,"
	    <<"   relative deviation = "
	    <<dabs(1.-m_xsinel/xsinel_test)*100.<<" %\n";   
  std::ofstream was;
  std::string filename = dirname+std::string("/Sigma.dat");
  was.open(filename.c_str());
  was<<"Test output for checking cross sections:\n"
     <<"   sigma_tot  = "<<m_xstot/1.e9<<" mb "
     <<"vs. "<<xstot_test/1.e9<<" mb,"
     <<"   relative deviation = "
     <<dabs(1.-m_xstot/xstot_test)*100.<<" %\n"   
     <<"   sigma_el   = "<<m_xsel/1.e9<<" mb "
     <<"vs. "<<xsel_test/1.e9<<" mb,"
     <<"   relative deviation = "
     <<dabs(1.-m_xsel/xsel_test)*100.<<" %\n"   
     <<"   sigma_inel = "<<m_xsinel/1.e9<<" mb "
     <<"vs. "<<xsinel_test/1.e9<<" mb,"
     <<"   relative deviation = "
     <<dabs(1.-m_xsinel/xsinel_test)*100.<<" %\n";   
  was.close();
}

//run_mode::code Cross_Sections::SelectCollisionMode() {
//  double random(ran->Get());
//  for (std::map<run_mode::code,double>::iterator miter=m_modemap.begin();
//       miter!=m_modemap.end();miter++) {
//    random -= miter->second;
//    if (random<=0) return miter->first;
//  }
//  return run_mode::unknown;
//}


