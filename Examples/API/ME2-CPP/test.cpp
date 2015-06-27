#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "AddOns/Python/MEProcess.H"

int main(int argc,char* argv[])
{
#ifdef USING__MPI
  MPI::Init(argc,argv);
#endif
  SHERPA::Sherpa *Generator(new SHERPA::Sherpa());
  // initialize the framework
  try {
    Generator->InitializeTheRun(argc,argv);

    // create a MEProcess instance
    MEProcess Process(Generator);
    Process.Initialize();

    for (size_t n(1);n<=Process.NumberOfPoints();++n) {
      // set momenta from file
      Process.SetMomenta(n);

      msg_Out()<<"Calculating matrix element values for phase space point "<<n<<":\n";
      msg_Out()<<*Process.GetAmp()<<std::endl;

      // compute flux factor -- fix
      double flux = Process.GetFlux();

      // get matrix elements
      double me    = Process.MatrixElement();
      double cs_me = Process.CSMatrixElement();

      // info strings
      std::string gen = Process.GeneratorName();

      size_t precision(msg_Out().precision());
      msg_Out().precision(16);
      msg_Out()<<"Matrix element generator:                        "<<gen  <<std::endl;
      msg_Out()<<"Color-summed matrix element:                     "<<cs_me<<std::endl;
      if (gen=="Comix")
        msg_Out()<<"Matrix element for specified color confiuration: "<<me <<std::endl;
      msg_Out()<<"Flux:                                            "<<flux <<std::endl;
      msg_Out().precision(precision);
    }
  }
  catch(ATOOLS::Exception exception) {
    std::terminate();
  }
  delete Generator;
#ifdef USING__MPI
  MPI::Finalize();
#endif
  return 0;
}

