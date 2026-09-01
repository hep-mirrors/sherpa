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
  MPI_Init(&argc, &argv);
#endif
  // initialize the framework
  try {
    SHERPA::Sherpa *Generator(new SHERPA::Sherpa(argc, argv));
    Generator->InitializeTheRun();

    // create a MEProcess instance
    MEProcess Process(Generator);
    Process.Initialize();

    msg_Out()<<"n="<<Process.NumberOfPoints()<<std::endl;
    for (size_t n(1);n<=Process.NumberOfPoints();++n) {
      // set momenta from file
      Process.ReadProcess(n);

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
    delete Generator;
  }
  catch (const ATOOLS::normal_exit& exception) {
    msg_Error() << exception << std::endl;
    ATOOLS::exh->Terminate(0);
  }
  catch (const ATOOLS::Exception& exception) {
    msg_Error() << exception << std::endl;
    ATOOLS::exh->Terminate(1);
  }
  catch (const std::exception& exception) {
    msg_Error() << exception.what() << std::endl;
    ATOOLS::exh->Terminate(1);
  }

#ifdef USING__MPI
  MPI_Finalize();
#endif

  return 0;
}

