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
  m_xstot(0.), m_slope(0.), m_xsinel(0.), m_xsel(0.)
{  }

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
  m_sigma_SD.FillGrids(&m_sigma_elastic);
  for (size_t i=0;i<2;i++) m_xsSD[i] = m_sigma_SD.GetXSec(i);
  m_xsDD = m_sigma_SD.GetXSec(2);

  msg_Info()<<"===========================================================\n"
	    <<"   sigma_tot                 = "<<m_xstot/1.e9<<" mb, (B = "<<m_slope<<")\n"
	    <<"   sigma_inel                = "<<m_xsinel/1.e9<<" mb\n"   
	    <<"   sigma_el                  = "<<m_xsel/1.e9<<" mb\n"   
	    <<"      test: int dt dsigma/dt = "<<m_sigma_elastic.Summed()/1.e9<<" mb\n"
	    <<"   sigma_SD0                 = "<<m_xsSD[0]/1.e9<<" mb\n"
	    <<"   sigma_SD1                 = "<<m_xsSD[1]/1.e9<<" mb\n"
	    <<"   sigma_DD                  = "<<m_xsDD/1.e9<<" mb.\n"
	    <<"===========================================================\n";
  MBpars.SetXSecs(this);
  m_sigma_inelastic.SetSigma(m_xsinel);
  m_sigma_elastic.SetSigma(m_xsel);
  m_sigma_SD.SetSigma(m_xsSD[0]+m_xsSD[1]+m_xsSD[2]);
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


