#include "HADRON_RESCATTERING/XSecs/BaryonMeson.H"
#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "NPion.H"
#include "NKaon.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

BaryonMeson::BaryonMeson() 
{
    m_test = true;
    if (m_test) 
    {
        BaryonMeson::Tests(); 
        exit(1); 
    }
}
BaryonMeson::~BaryonMeson() {}
double BaryonMeson::XStot(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B, const double & s)
{   
    // msg_Out() << "BaryonMeson::XStot Function is called" << std::endl;
    //Scoped_Settings s{Settings::GetMainSettings()["FormFactors"]};
    //m_bb_xstot = s["BBScattering"].SetDefault(scatmodel::off).Get<scatmodel::code>();
    // if (!(A.IsBaryon() && B.IsBaryon())) return 0.;

    if ( (A.IsMeson() && B.IsBaryon()) || (A.IsBaryon() && B.IsMeson()))
    {

        //p-Kaon^- or p-Kaon^+
        if ((A.Kfcode()==2212 && B.Kfcode()==321) || (A.Kfcode()==321 && B.Kfcode()==2212))
        {

            if ((A.IsAnti() && B.Kfcode() == 2212) || (A.Kfcode()==2212 && B.IsAnti()))
            {
                // msg_Out() << "BaryonMeson::XStot Function is called for proton- Kaon-" << std::endl;
                return m_NK.pKaonMinus(s);
                
                //return m_hpr1r2.xs_tot(hpr1r2::pKMinus,s);

            }
            // msg_Out() << "BaryonMeson::XStot Function is called for proton- Kaon^+" << std::endl;
        }





        if (( A.Kfcode() == 2212 && B.Kfcode()==111 ) || (A.Kfcode() == 111 && B.Kfcode() == 2212))
        {
        //   msg_Out() << "No reference for Proton-Pion0 " << std::endl;
        //   return m_NP.pPiZeroWignerEckart(s);
            return m_NP.pPiZero(s);

        }
        //p-Pi^- or p-Pi^+
        if ((A.Kfcode()==2212 && B.Kfcode()==211) || (A.Kfcode()==211 && B.Kfcode()==2212)) 
        {
            // msg_Out() << "BaryonMeson::XStot Function is called for pPi+" << std::endl;

            if ((A.IsAnti() && B.Kfcode() == 2212) || (A.Kfcode()==2212 && B.IsAnti()))
            {
                // msg_Out() << "BaryonMeson::XStot Function is called for pPi-" << std::endl;
                return m_NP.pPiMinus(s);

            }
            // return m_BM.pPiPlus(s);
            return m_NP.pPiPlus(s);
        }
        //n-Pi^- or n-Pi^+
        if ((A.Kfcode()==2112 && B.Kfcode()==211) || (A.Kfcode()==211 && B.Kfcode()==2112)) 
        {
            msg_Out() << "BaryonMeson::XStot Function is called for nPi+" << std::endl;
            return m_NP.nPiPlus(s);

            if ((A.IsAnti() && B.Kfcode() == 2112) || (A.Kfcode()==2112 && B.IsAnti()))
            {
                msg_Out() << "BaryonMeson::XStot Function is called for nPi-" << std::endl;

                return m_NP.nPiMinus(s);

            }
            // return m_BM.nPiPlus(s);
        }

        if ( A.Kfcode() == 2112 && B.Kfcode()==111 || A.Kfcode() == 111 && B.Kfcode()==2112)
        {
          msg_Out() << "Neutron-Pion0 " << std::endl;
          return m_NP.nPiZero(s);

        }
        // //pK- 
        // if ( (A.Kfcode() == 2212 && B.Kfcode()==321) || (A.Kfcode() == 321 && B.Kfcode()==2212))
        // {
        //     msg_Out() << "BaryonMeson::XStot Function is called for pK-" << std::endl;

        //     return m_BM.pKMinus(s);
        // }
        // //pKBar0
        // if ( (A.Kfcode() == 2212 && B.Kfcode()==311) || (A.Kfcode() == 311 && B.Kfcode()==2212))
        // {
        //     msg_Out() << "BaryonMeson::XStot Function is called for pKBar0" << std::endl;

        //     return m_BM.pKBarZero(s);
        // }

        //nK- 
        if ( (A.Kfcode() == 2112 && B.Kfcode()==321) || (A.Kfcode() == 321 && B.Kfcode()==2112))
        {
            msg_Out() << "BaryonMeson::XStot Function is called for nK-" << std::endl;

            return m_NK.nKaonMinus(s);
        }
        // //nKBar0
        // if ( (A.Kfcode() == 2112 && B.Kfcode()==311) || (A.Kfcode() == 311 && B.Kfcode()==2112))
        // {
        //     msg_Out() << "BaryonMeson::XStot Function is called for nKBar0" << std::endl;

        //     return m_BM.nKBarZero(s);
        // }


    m_s   = s;
    m_mA = A.HadMass(); m_mA2 = sqr(m_mA);
    m_mB = B.HadMass(); m_mB2 = sqr(m_mB);
    msg_Out() << "No test has been called" << std::endl;

    }
    else
    {
        return 0.;
    }
}


// double BaryonMeson::XSel(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B, const double & s)
// {
//     msg_Out() << "BaryonMeson::XSel Function is called" << std::endl;

//     if ( (A.IsMeson() && B.IsBaryon()) || (A.IsBaryon() && B.IsMeson()))
//     {
//         if ( (A.Kfcode() == 2212) && (B.Kfcode() == 211) || (A.Kfcode() == 2112 && (B.Kfcode() ==211)))
//         {
//             return m_BM.pAndN_Pi(s);
//         }
//     }
// }


void BaryonMeson::Tests() {

  msg_Out() << "BaryonMeson::Test  called" << std::endl;
//Must Fix the test function for different interactions. 

  size_t bins = 10000;
  double pmin = 0. , pmax = 30., pinc = (pmax-pmin)/double(bins);
  map<string,Histogram *>  histos;
  histos["protonPionZero_Total"] = new Histogram(0,pmin,pmax,bins);
  Flavour flA(kf_p_plus);
  double  plab, s;
  // Flavour flB(kf_pi_plus);
  Flavour flB(kf_pi);
//   flB = flB.Bar();
 
  // Flavour flA(kf_n);
  // flB = flB.Bar();

  for (int i=0;i<bins;i++) 
  {
    plab   = pmin+i*pinc;
    s      = ( sqr(flA.HadMass()) + sqr(flB.HadMass()) +
	       2.*flA.HadMass()*sqrt(sqr(flB.HadMass())+sqr(plab)) );
    histos["protonPionZero_Total"]->Insert(i,   XStot(flA,flB,s));
  }

  Histogram * histo;
  std::string name;
  for (std::map<std::string,Histogram *>::iterator
	 hit=histos.begin();hit!=histos.end();hit++) {
    histo = hit->second;
    name  = std::string("XSecs_BaryonMesonTest/")+hit->first+std::string(".dat");
    histo->Output(name);
    delete histo;
  }
  histos.clear();
}



















