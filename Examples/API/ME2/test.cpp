#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "AddOns/Python/MEProcess.H"

int main(int argc,char* argv[])
{
  try {
    // initialize the framework
    SHERPA::Sherpa *Generator(new SHERPA::Sherpa());
    Generator->InitializeTheRun(argc,argv);

    // create a MEProcess instance
    MEProcess Process(Generator);
    Process.Initialize();

    for (size_t n(1);n<=Process.NumberOfPoints();++n) {
      // generate one specific random color configuration
      // (necessary if comix is used)
      double fac(1.);
      if (Process.HasColorIntegrator()) fac=Process.GenerateColorPoint();

      // set momenta from file
      Process.SetMomenta(n);

      msg_Out()<<"Calculating matrix element values for phase space point "<<n<<":\n";
      msg_Out()<<*Process.GetAmp()<<std::endl;

      // compute flux factor -- fix
      double flux = Process.GetFlux();

      // get matrix elements
      double me    = Process.MatrixElement()/fac;
      double cs_me = Process.CSMatrixElement();

      // info strings
      std::string gen = Process.GeneratorName();

      size_t precision(msg_Out().precision());
      msg_Out().precision(16);
      msg_Out()<<"Matrix element generator:                        "<<gen  <<std::endl;
      msg_Out()<<"Color-summed matrix element:                     "<<cs_me<<std::endl;
      if (gen=="Comix") {
        msg_Out()<<"Matrix element for specified color confiuration: "<<me <<std::endl;
        msg_Out()<<"Number of colour configurations:                 "<<fac<<std::endl;
      }
      msg_Out()<<"Flux:                                            "<<flux <<std::endl;
      msg_Out().precision(precision);
    }

    delete Generator;
  }
  catch(ATOOLS::Exception a) {
    msg_Error()<<a<<std::endl;
    std::terminate();
  }
  catch(char const* a) {
    std::cout << a << std::endl;
    exit(1);
  }
}

