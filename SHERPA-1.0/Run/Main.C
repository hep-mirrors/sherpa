#include "Sherpa.H"
#include "prof.hh"

//#define _USE_MPI_
#ifdef _USE_MPI_
#include <mpi++.h>
#endif

extern "C" {
  void apainit_();
  void aparun_();
}

using namespace SHERPA;

int main(int argc,char* argv[]) {  
  /*
  // *AS* simple calls:
  cout<<" apainit "<<endl;
  apainit_();
  for (int i=1;i<=20000;++i) {
    //    if (i%100==0) 
    cout<<" Event "<<i<<endl;
    aparun_();
  }
  */

  /*
  set_prof();
  */
#ifdef _USE_MPI_
  MPI::Init(argc, argv);
#endif

  Sherpa Generator;
  Generator.Init();
  if (Generator.CrossSections())
    Generator.GenerateEvents();

#ifdef _USE_MPI_
  MPI::Finalize();
#endif

  /*
  std::ofstream file("profile.out");
  print_profile( file );
  */
}



