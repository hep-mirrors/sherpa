#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "SHRiMPS/Cross_Sections/Sigma_Total.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
//#include "SHRiMPS/Cross_Sections/Sigma_SD.H"
//#include "SHRiMPS/Cross_Sections/Sigma_DD.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
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
  Sigma_Inelastic sigma_inel;
  m_xsinel = sigma_inel.Calculate();
  Sigma_Elastic sigma_el;
  m_xsel   = sigma_el.Calculate();
  msg_Info()<<"===========================================================\n"
	    <<"   sigma_tot  = "<<m_xstot/1.e9<<" mb, (B = "<<m_slope<<")\n"
	    <<"   sigma_inel = "<<m_xsinel/1.e9<<" mb\n"   
	    <<"   sigma_el   = "<<m_xsel/1.e9<<" mb\n"   
	    <<"===========================================================\n";
  sigma_el.FillDifferentialGrids();
  MBpars.SetXSecs(this);
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


