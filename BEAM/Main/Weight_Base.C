#include "BEAM/Main/Weight_Base.H"

using namespace BEAM;
using namespace ATOOLS;

Weight_Base::Weight_Base(Kinematics_Base * kinematics) :
  p_kinematics(kinematics), p_beams(kinematics->GetBeams()),
  m_weight(1.), m_keyid(p_kinematics->KeyId()) {};

Weight_Base::~Weight_Base() {};


/*
// The following only if we ever go back to polarised beams -
// and in my mind this will need to be debugged.

double Beam_Spectra_Handler::Weight(int * pol_types, double *dofs)
{
  double weight = 1.;
  for (int i=0;i<2;++i) {
    if (p_BeamBase[i]->PolarisationOn()) {
      if (pol_types[i]!=99) {
	double hel=(double)pol_types[i];
	double pol=p_BeamBase[i]->Polarisation();
	double dof=dofs[i];
	if (hel*pol>0.)
	  weight*=(1.+dabs(pol)*(dof-1.))/dof;
	else
	  weight*=(1.-dabs(pol))/dof;
	//assuming 2 degrees of freedom
	//	weight*=dabs(hel+pol)/2.;
      }
      else msg_Out()<<"ERROR: unpolarised cross section for polarised beam!! "<<endl;
    } // end of polarisation conditional
  } // end of for loop
  return weight;
}
*/
  
