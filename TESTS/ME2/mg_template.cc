#include "rambo.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <fstream>
#include <cstring>
#include <math.h>
#include "CPPProcess.h"

int main(int argc,char* argv[])
{
  try
    {
      double               energy = ${cms};
      double               s      = energy*energy;
      double               weight = 1.0;
      int                  ncalls = ${ncalls};
      std::vector<double*> p;

      ostringstream        o_stream;
      std::ofstream        o_file("mg_results.dat");
      o_stream.flags(std::ios::scientific);
      o_stream.precision(std::numeric_limits<double>::digits10 + 1);
      
      CPPProcess           mg_proc;
      mg_proc.initProc("param_card.dat");
      mg_proc.setInitial( ${in_flav0} , ${in_flav1} );


      for (int i(0);i< ncalls;i++)
	{
	  p = get_momenta(mg_proc.ninitial, energy, mg_proc.getMasses(), weight);
	  mg_proc.setMomenta(p);
	  mg_proc.sigmaKin();
	  double MGME = mg_proc.sigmaHat();
	  o_stream.str("");
	  for(int i=0; i<${n_flavs}; i++)
	    o_stream << p[i][0] << "   " << p[i][1] << "   " << p[i][2] << "   " << p[i][3] << "   ";
	  o_stream << mg_proc.sigmaHat() << std::endl;
	  o_file << o_stream.str();
	}
      o_file.close();
    }
  catch(char const* a)
    {
      std::cout << "Exception caught: " << a << std::endl;
      return 1;
    }
  return 0;
}
