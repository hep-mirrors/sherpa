//****************************************************************************//
//                                                                            //
//    cdemo0_rcl.cpp                                                          //
//    is part of RECOLA2 (REcursive Computation of One Loop Amplitudes)2      //
//                                                                            //
//    Copyright (C) 2016,2017 Ansgar Denner, Jean-Nicolas Lang and            //
//                            Sandro Uccirati                                 //
//                                                                            //
//    RECOLA2 is licenced under the GNU GPL version 3,                        //
//    see COPYING for details.                                                //
//                                                                            //
//****************************************************************************//
// PARTICLES                                                                  //
// Scalars in the SM:    'H', 'G0', 'G+', 'G-'                                //
// Scalars in the HSESM: 'Hl', 'Hh', 'G0', 'G+', 'G-'                         //
// Scalars in the THDM:  'Hl', 'Hh', 'Ha', 'H+', 'H-', 'G0', 'G+', 'G-'       //
// Vector bosons:        'g', 'A', 'Z', 'W+', 'W-'                            //
// leptons:              'nu_e', 'nu_e~', 'e-', 'e+',                         //
//                       'nu_mu', 'nu_mu~', 'mu-', 'mu+',                     //
//                       'nu_tau', 'nu_tau~', 'tau-', 'tau+'                  //
// quarks:               'u', 'u~', 'd', 'd~',                                //
//                       'c', 'c~', 's', 's~',                                //
//                       't', 't~', 'b', 'b~'                                 //
//****************************************************************************//


// Each C++ program, which uses RECOLA must have:
#include <iostream>
#include "recola.hpp"
#include "yaml-cpp/yaml.h"

void readcard(std::vector<std::vector<std::vector<double> > > &moms)
{
  YAML::Node file = YAML::LoadFile("../../Sherpa.yaml");
  for(const auto& t : file){
    std::cout<<"    "<<t.first << ": " << t.second <<std::endl;
    std::string name = t.first.as<std::string>();
    if(name=="MOMENTA")
      // list of list of vectors
      moms = t.second.as<std::vector<std::vector<std::vector<double> > > >();
  }
}


int main(int argc, char *argv[])
{

  std::vector<std::vector<std::vector<double> > > moms;
  readcard(moms);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Step 1                                                                     //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Set input values for the computation.                                      //
// General methods for arbitrary models are defined in "input_rcl.f90",       //
// whereas model specific ones are defined in "recola1_interface_rcl.90".     //
// Methods dedicated to the THDM/HSESM are listed in                          //
// "extended_higgs_interface_rcl.f90".                                        //
// Since all variables have default values, this step is optional.            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

// The standard output is selected
  Recola::set_output_file_rcl("*");

// Let's print the squared amplitude
  Recola::set_print_level_squared_amplitude_rcl(1);

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Step 2                                                                     //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// In step 2 processes to be computed are defined. Optionally, specific powers//
// in fundamental couplings, or polarizations can be enforced using methods   //
// defined in the modules "process_definition_rcl" and "recola1_interface_rcl"//
// The processes are defined by calling "define_process_rcl" which requires a //
// unique process number, a process string and the order of computation. Note //
// that by default all powers in fundamental couplings are selected.          //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

// We define a process at NLO:
  Recola::define_process_rcl(1,"e- e+ -> G+ G-","LO");

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Step 3                                                                     //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// The skeleton of the recursive procedure is built for all defined           //
// processes, by calling the subroutine "generate_processes_rcl".             //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  Recola::generate_processes_rcl();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Step 4                                                                     //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// The fourth step is the actual computation of processes.                    //
// Each process defined at step 2 can be computed at this stage for given     //
// phase-space points provided by the user.                                   //
// The computation of (squared) amplitudes is carried out by calling the      //
// subroutine "compute_process_rcl".                                          //
// In the module "process_computation_rcl" other useful                       //
// methods are defined, which allow to get the value for specific             //
// contributions to amplitudes, squared amplitudes, and                       //
// Born colour- and/or spin-correlated squared amplitudes.                    //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


  for(std::vector<std::vector<double> > mlist : moms){
    // mlist is a list of momenta for a single ps-point
    double p[4][4];
    size_t n_part{4}, n_dim{4};
    for(size_t i{0}; i<n_part;++i){
      for(size_t j{0}; j<n_dim;++j){
	// (the first entry is the flavour code, so you can skip it)
	p[i][j] = mlist[i][j+1];
      }
    }
    Recola::compute_process_rcl(1,p,"LO");
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Step 5                                                                     //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
// Calling reset_recola_rcl (module reset_rcl), deallocates all processes     //
// generated in the previous steps and allows for the next call of Recola.    //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  Recola::reset_recola_rcl();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  return 0;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
