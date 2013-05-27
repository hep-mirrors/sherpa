#include "OpenLoops_Virtual.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Poincare.H"


using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace OpenLoops {

OpenLoops_Interface* OpenLoops_Virtual::s_interface=NULL;

OpenLoops_Virtual::OpenLoops_Virtual(const Process_Info& pi,
                                     const Flavour_Vector& flavs,
                                     Amp2Func amp2,
                                     PermutationFunc permutationfunc,
                                     std::vector<int> permutation,
                                     std::string functag) :
  Virtual_ME2_Base(pi, flavs),
  m_amp2(amp2), m_permutationfunc(permutationfunc),
  m_permutation(permutation), m_analyse(true)
{
  Data_Reader reader(" ",";","#","=");
  m_write_points=reader.GetValue<int>("OL_WRITE_POINTS", 0);
  if (m_write_points) {
    m_points_file.open(string("points."+functag+".txt").c_str());
  }
  // There is a bool-bit: m_analyse that triggers the stability analysis.
  // Here you can add more histograms for the stability analysis.  The arguments
  // to initialise them are (0=lin, xmin, xmax, number of bins).  Since the
  // histograms are put into a map, make sure the names differ!
  // To insert values, see below.  The setup is also such that all histograms
  // in the map are automatically written out in the subdirectory
  // "$pwd/Stability_Analysis/", if it exists.
  if (m_analyse) {
    m_histograms[string("log_V_by_B")]     = new Histogram(0,-5.,2.,56);
    m_histograms[string("log_V_by_B_reg")] = new Histogram(0,-5.,2.,56);
    m_histograms[string("log_V_by_B_wt")]  = new Histogram(0,-5.,2.,56);
  }

  m_oew=pi.m_oew;
  m_oqcd=pi.m_oqcd;
}

OpenLoops_Virtual::~OpenLoops_Virtual()
{
  if (m_write_points) {
    m_points_file.close();
  }
  // Here the histograms are written out.  I'll add mean, right variance etc.
  if (m_analyse) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = string("Stability_Analysis/")+hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
}


void OpenLoops_Virtual::Calc(const Vec4D_Vector& momenta) {
  Vec4D_Vector m_moms(momenta);

  double alpha_QED=AlphaQED();
  double alpha_S=AlphaQCD();
  s_interface->OpenLoopsInit(m_mur2, alpha_QED, alpha_S);

  double M2L0;
  vector<double> M2L1(3), M2L2(5), IRL1(3), IRL2(5);

  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
  }
  m_permutationfunc(&m_permutation[0]);
  m_amp2(&m_moms[0][0], &M2L0, &M2L1[0], &IRL1[0], &M2L2[0], &IRL2[0]);
  if (msg_LevelIsDebugging()) {
    timing->Stop();
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
               <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
  }

  double B(M2L0), V_finite(M2L1[0]), V_eps(M2L1[1]), V_eps2(M2L1[2]);

  // factor which by Sherpa convention has to be divided out at this stage
  double factor=B*alpha_S/2.0/M_PI;

  m_born=B;
  m_res.Finite()=(V_finite/factor);
  m_res.IR()=(V_eps/factor);
  m_res.IR2()=(V_eps2/factor);

  if (m_write_points) {
    m_points_file<<"2 -> "<<m_flavs.size()-2<<" 1 "<<alpha_S<<" "<<sqrt(m_mur2)<<" "<<alpha_QED<<" "<<sqrt(m_mur2)<<endl;
    for (size_t i=0; i<momenta.size(); ++i) {
      m_points_file<<m_flavs[i].HepEvt();
      for (size_t j=0; j<4; ++j) {
        m_points_file<<" "<<momenta[i][j];
      }
      m_points_file<<endl;
    }
    m_points_file<<B<<" "<<V_finite<<" "<<V_eps<<" "<<V_eps2<<endl;
  }
  // Here the histograms are filled.  The syntax is with 
  // Insert(xbin, value), where the value, by default is 1.
  if (m_analyse) {
    double logVbyB(log(dabs(V_finite/B))/log(10.));
    m_histograms[string("log_V_by_B")]->Insert(logVbyB);    
    m_histograms[string("log_V_by_B_wt")]->Insert(logVbyB,B);    
    double logVbyBren(log(dabs(V_finite/(B+.00001*dabs(V_finite))))/log(10.));
    m_histograms[string("log_V_by_B_reg")]->Insert(logVbyBren);    
  }
}

}
