#include "Coulomb.H"
#include "Run_Parameter.H"
#include "Model.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Coulomb::Coulomb(Single_Amplitude* _first_ampl) : first_ampl(_first_ampl)
{
  short int ngraph = 0;
  Single_Amplitude* m = first_ampl;
  while (m) {ngraph++;  m = m->Next;}
  
  COULM = new double*[ngraph];

  for (short int i=0;i<ngraph;i++) {
    COULM[i] = new double[ngraph];
    for (short int j=0;j<ngraph;j++) COULM[i][j] = 0.;
  }
  Build_Matrix();
}


void Coulomb::Build_Matrix()
{
  icoulomb = 0;
  if (AORGTOOLS::rpa.me.UsingCoulombCorr()) {
    Single_Amplitude* m;
    int wcount;
    //tops later on
    m = first_ampl;
    while (m) {
      list<Pfunc*>* pl = m->GetPlist();
      wcount = 0;
      for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
	Pfunc* p = *pit;
	if ((p->fl==Flavour(kf::W)) ||  (p->fl==Flavour(kf::W).Bar())) wcount++;
      }
      if (wcount==2) {
	icoulomb = 1;
	m->Set_Coulomb();
      }
      m = m->Next;
    }
    Single_Amplitude* m1;
    m = first_ampl;
    int c1,c2;
    c1 = 0;
    while (m) {
      if (m->Get_Coulomb()) {
	m1 = m;
	c2 = c1;
	while (m1) {
	  COULM[c1][c2] = COULM[c2][c1] = double(m1->Get_Coulomb());
	  m1 = m1->Next;
	  c2++;
	}
      }
      m = m->Next;
      c1++;
    }
  }
}

void Coulomb::Calculate(int* b,Vec4D* mom)
{
  Single_Amplitude* m;
  short int i;
  m = first_ampl;
  while (m) {
    if (m->Get_Coulomb()) break;
    m = m->Next;
  }
  double mp,mm;
  mp = mm = 0.;
  list<Pfunc*>* pl = m->GetPlist();
  Vec4D sump;
  for (list<Pfunc*>::iterator pit=pl->begin();pit!=pl->end();++pit) {
    Pfunc* pp = *pit;  
    if (pp->fl==Flavour(kf::W)) {
      sump = Vec4D(0.,0.,0.,0.);
      for (i=1;i<pp->argnum;i++) sump += b[pp->arg[i]]*mom[pp->arg[i]];
      mp = sqrt(sump.Abs2());
    }
    if (pp->fl==Flavour(kf::W).Bar()) {
      sump = Vec4D(0.,0.,0.,0.);
      for (i=1;i<pp->argnum;i++) sump += b[pp->arg[i]]*mom[pp->arg[i]];
      mm = sqrt(sump.Abs2());
    }
  }
  double s      = (mom[0]+mom[1]).Abs2();
  double MW     = AORGTOOLS::rpa.consts.Mass(Flavour(kf::W),s);
  double GW     = AORGTOOLS::rpa.consts.Width(Flavour(kf::W),s);
  double beta   = sqrt((s-sqr(mp-mm))*(s-sqr(mp+mm)))/s;
  Complex betam = sqrt(Complex(1.-4.*MW*MW/s,4*MW*GW/s)); 
  double delta  = dabs(sqr(mp)-sqr(mm))/s;
  double arg    = 2./M_PI*atan((abs((betam+delta)*conj(betam+delta))-sqr(beta))/
			       (2*beta*imag(betam)));
  fcoul = mo->Aqed()*M_PI/(2.*beta)*(1.-arg);
}







