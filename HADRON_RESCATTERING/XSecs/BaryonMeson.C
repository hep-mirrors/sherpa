#include "HADRON_RESCATTERING/XSecs/BaryonMeson.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "NPion.H"

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

double BaryonMeson::XStot(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B, const double & s)
{   
    // msg_Out() << "BaryonMeson::XStot Function is called" << std::endl;
    //Scoped_Settings s{Settings::GetMainSettings()["FormFactors"]};
    //m_bb_xstot = s["BBScattering"].SetDefault(scatmodel::off).Get<scatmodel::code>();
    // if (!(A.IsBaryon() && B.IsBaryon())) return 0.;

    if ( (A.IsMeson() && B.IsBaryon()) || (A.IsBaryon() && B.IsMeson()))
    {
        //p-Pi^- or p-Pi^+
        if ((A.Kfcode()==2212 && B.Kfcode()==211) || (A.Kfcode()==211 && B.Kfcode()==2212)) 
        {
            // msg_Out() << "BaryonMeson::XStot Function is called for pPi+" << std::endl;

            if ((A.IsAnti() && B.Kfcode() == 2212) || (A.Kfcode()==2212 && B.IsAnti()))
            {
                msg_Out() << "BaryonMeson::XStot Function is called for pPi-" << std::endl;

                // return m_BM.pPiMinus(s); 
            }
            // return m_BM.pPiPlus(s);
            return m_NP.pPiPlus(s);
        }
        // //n-Pi^- or n-Pi^+
        // if ((A.Kfcode()==2112 && B.Kfcode()==211) || (A.Kfcode()==211 && B.Kfcode()==2112)) 
        // {
        //     msg_Out() << "BaryonMeson::XStot Function is called for nPi+" << std::endl;

        //     if ((A.IsAnti() && B.Kfcode() == 2112) || (A.Kfcode()==2112 && B.IsAnti()))
        //     {
        //         msg_Out() << "BaryonMeson::XStot Function is called for nPi-" << std::endl;

        //         return m_BM.nPiMinus(s); 
        //     }
        //     return m_BM.nPiPlus(s);
        // }
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

        // //nK- 
        // if ( (A.Kfcode() == 2112 && B.Kfcode()==321) || (A.Kfcode() == 321 && B.Kfcode()==2112))
        // {
        //     msg_Out() << "BaryonMeson::XStot Function is called for nK-" << std::endl;

        //     return m_BM.nKMinus(s);
        // }
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
  double pmin = 0., pmax = 30., pinc = (pmax-pmin)/double(bins);
  map<string,Histogram *>  histos;
  histos["pPi_total"]      = new Histogram(0,pmin,pmax,bins);
  Flavour flA(kf_p_plus);
  double  plab, s;
  Flavour flB(kf_pi_plus);

  for (int i=0;i<bins;i++) 
  {
    plab   = pmin+i*pinc;
    s      = ( sqr(flA.HadMass()) + sqr(flB.HadMass()) +
	       2.*flA.HadMass()*sqrt(sqr(flB.HadMass())+sqr(plab)) );
    histos["pPi_total"]->Insert(i,   XStot(flA,flB,s));
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



















