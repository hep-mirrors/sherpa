#include "Standard_Selector.H"
#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Flavour.H"

using namespace APHYTOOLS;
using namespace AMATOOLS;

/*--------------------------------------------------------------------

  Energy Selector

  --------------------------------------------------------------------*/

Energy_Selector::Energy_Selector(int _nin,int _nout, Flavour * _fl) {
  name = std::string("Energy_Selector"); 
  nin  = _nin; nout = _nout; n = nin+nout;
  fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms();
  emin  = new double[n];
  emax  = new double[n];
  value = new double[n];
  for (int i=0;i<n;i++) { emin[i] = 0.; emax[i] = E; }
  sel_log = new Selector_Log(name);
}

Energy_Selector::~Energy_Selector() {
  //  delete [] sel_log;
  delete [] emin;
  delete [] emax;
  delete [] value;
}


bool Energy_Selector::Trigger(const vec4d * mom) 
{
  double ei;
  for (int i=nin;i<n;i++) {
    ei = value[i] = mom[i][0];
    if (sel_log->Hit( ((ei<emin[i]) || (ei>emax[i])) )) return 0;
  }
  return 1;
}

double * Energy_Selector::ActualValue() { return value; }

void Energy_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<n-1;i++) {
    cuts->energymin[i] = Max(emin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(emax[i],cuts->energymax[i]);
    for (int j=i+1;j<n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}

void Energy_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=0;i<n-1;i++) {
    cuts->energymin[i] = Max(emin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(emax[i],cuts->energymax[i]);
    for (int j=i+1;j<n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}
 
void Energy_Selector::SetRange(std::vector<APHYTOOLS::Flavour> crit,double _min, 
			       double _max=0.5*AORGTOOLS::rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in Energy_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=nin;i<n;i++) {
    if ( (crit[0].Includes(fl[i])) || ((crit[0].bar()).Includes(fl[i]) ) ) {
      emin[i] = AMATOOLS::Max(_min,fl[i].mass()); 
      emax[i] = AMATOOLS::Min(_max,0.5*AORGTOOLS::rpa.gen.Ecms());
      AORGTOOLS::msg.Debugging()<<"Set e-Range for "<<fl[i]<<" : "
				<<emin[i]<<" ... "<<emax[i]<<endl;
    }
  }
}



/*--------------------------------------------------------------------

  PT Selector

  --------------------------------------------------------------------*/

PT_Selector::PT_Selector(int _nin,int _nout, Flavour * _fl) {
  name = std::string("PT_Selector"); 
  nin  = _nin; nout = _nout; n = nin+nout;
  fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms();
  ptmin  = new double[n];
  ptmax  = new double[n];
  value = new double[n];
  for (int i=0;i<n;i++) { ptmin[i] = 0.; ptmax[i] = E; }
  sel_log = new Selector_Log(name);
}

PT_Selector::~PT_Selector() {
  //  delete [] sel_log;
  delete [] ptmin;
  delete [] ptmax;
  delete [] value;
}


bool PT_Selector::Trigger(const vec4d * mom) 
{
  double pti;
  for (int i=nin;i<n;i++) {
    pti = value[i] = sqrt(sqr(mom[i][1]) + sqr(mom[i][2]));
    if (sel_log->Hit( ((pti<ptmin[i]) || (pti>ptmax[i])) )) return 0;
  }
  return 1;
}

double * PT_Selector::ActualValue() { return value; }

void PT_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<n-1;i++) {
    cuts->energymin[i] = Max(ptmin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(ptmax[i],cuts->energymax[i]);
    for (int j=i+1;j<n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}

void PT_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=0;i<n-1;i++) {
    cuts->energymin[i] = Max(ptmin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(ptmax[i],cuts->energymax[i]);
    for (int j=i+1;j<n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}
 
void PT_Selector::SetRange(std::vector<APHYTOOLS::Flavour> crit,double _min, 
			       double _max=0.5*AORGTOOLS::rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in PT_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=nin;i<n;i++) {
    if ( (crit[0].Includes(fl[i])) || ((crit[0].bar()).Includes(fl[i]) ) ) {
      ptmin[i] = AMATOOLS::Max(_min,fl[i].mass()); 
      ptmax[i] = AMATOOLS::Min(_max,0.5*AORGTOOLS::rpa.gen.Ecms());
      AORGTOOLS::msg.Debugging()<<"Set PT-Range for "<<fl[i]<<" : "
				<<ptmin[i]<<" ... "<<ptmax[i]<<endl;
    }
  }
}




/*--------------------------------------------------------------------

  Rapidity Selector

  --------------------------------------------------------------------*/


Rapidity_Selector::Rapidity_Selector(int _nin,int _nout, Flavour * _fl) {
  name = std::string("Rapidity_Selector"); 
  nin  = _nin; nout = _nout; n = nin+nout;
  fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms()/2;
  double pl;

  ymin  = new double[n];
  ymax  = new double[n];
  value = new double[n];
  for (int i=0;i<n;i++) {
    pl      = sqrt(E*E-sqr(_fl[i].mass()));
    ymax[i] = log( (E+pl)/(E-pl) );
    ymin[i] = -ymax[i];
  }
  sel_log = new Selector_Log(name);
}

Rapidity_Selector::~Rapidity_Selector() {
  //  delete [] sel_log;
  delete [] ymin;
  delete [] ymax;
  delete [] value;
}


bool Rapidity_Selector::Trigger(const vec4d * mom) 
{
  double yi;
  for (int i=nin;i<n;i++) {
    yi = value[i] = log( (mom[i][0]+mom[i][3])/(mom[i][0]-mom[i][3]) );
    if (sel_log->Hit( ((yi<ymin[i]) || (yi>ymax[i])) )) return 0;
  }
  return 1;
}

double * Rapidity_Selector::ActualValue() { return value; }

void Rapidity_Selector::BuildCuts(Cut_Data * cuts) {}

void Rapidity_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {}
 


void Rapidity_Selector::SetRange(std::vector<APHYTOOLS::Flavour> crit,double _min, 
			       double _max=0.5*AORGTOOLS::rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in Rapidity_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  double E = AORGTOOLS::rpa.gen.Ecms()/2;
  double pl,y;

  for (int i=nin;i<n;i++) {
    pl      = sqrt(E*E-sqr(fl[i].mass())); 
    y       = log((E+pl)/(E-pl));
    ymin[i] = AMATOOLS::Max(_min,-y);
    ymax[i] = AMATOOLS::Min(_max,y);
    AORGTOOLS::msg.Debugging()<<"Set y-Range for "<<fl[i]<<" : "
			      <<ymin[i]<<" ... "<<ymax[i]<<endl;
  }
}





Angle_Selector::Angle_Selector(int _nin,int _nout, Flavour * _fl) {
  name = std::string("Angle_Selector"); 
  nin  = _nin; nout = _nout; n = nin+nout;
  fl   = _fl;

  cosmin = new double*[n];
  cosmax = new double*[n];
  value  = new double[n*n];
  for (int i=0;i<n;i++) { cosmin[i] = new double[n]; cosmax[i] = new double[n]; }
  for (int i=0;i<n;i++) {
    for (int j=i+1;j<n;j++) {
      cosmin[i][j] = cosmin[j][i] = -1.; 
      cosmax[i][j] = cosmax[j][i] =  1.; 
    }
  }
    
  sel_log = new Selector_Log(name);
}

Angle_Selector::~Angle_Selector() {
  //  delete [] sel_log;
  for (int i=0;i<n;i++) {
    delete [] cosmin[i];
    delete [] cosmax[i];
  }
  delete [] cosmin;
  delete [] cosmax;
  delete [] value;
}

bool Angle_Selector::Trigger(const vec4d * mom) 
{
  double cosij;
  for (int i=0;i<n;i++) {
    for (int j=i+1;j<n;j++) {
      cosij = value[n*i+j] = 
	vec3d(mom[i])*vec3d(mom[j])/(vec3d(mom[i]).abs()*vec3d(mom[j]).abs());
      if (sel_log->Hit( ((cosij < cosmin[i][j]) || 
			 (cosij > cosmax[i][j])) )) return 0;
    }
  }
  return 1;
}

double * Angle_Selector::ActualValue() { return value; }

void Angle_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<n-1;i++) {
    for (int j=i+1;j<n;j++) {
      cuts->cosmin[i][j] = cuts->cosmin[j][i] = 
	Max(cosmin[i][j],cuts->cosmin[i][j]);
      cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	Min(cosmax[i][j],cuts->cosmax[i][j]);
      cuts->scut[i][j]   = cuts->scut[j][i]   = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}

void Angle_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=0;i<n-1;i++) {
    for (int j=i+1;j<n;j++) {
      cuts->cosmin[i][j] = cuts->cosmin[j][i] = 
	Max(cosmin[i][j],cuts->cosmin[i][j]);
      cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	Min(cosmax[i][j],cuts->cosmax[i][j]);
      cuts->scut[i][j]   = cuts->scut[j][i]   = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}


void Angle_Selector::SetRange(std::vector<APHYTOOLS::Flavour> crit,
			      double _min, double _max)
{
  if (crit.size() != 2) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in Angle_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  for (int i=nin;i<n;i++) {
    for (int j=i+1;j<n;j++) {
      if ( ((crit[0].Includes(fl[i])) && (crit[1].Includes(fl[j])) ) || 
	   ((crit[0].Includes(fl[j])) && (crit[1].Includes(fl[i])) ) ) {
	cosmin[i][j] = cosmin[j][i] = AMATOOLS::Max(_min,-1.); 
	cosmax[i][j] = cosmax[j][i] = AMATOOLS::Min(_max,1.); 
	AORGTOOLS::msg.Debugging()<<"Set cos-Range for "<<fl[i]<<"/"<<fl[j]<<" : "
				  <<cosmin[i][j]<<" ... "<<cosmax[i][j]<<endl;
      }
    }
  }
}

void Angle_Selector::SetRange(std::vector<APHYTOOLS::Flavour> crit,int beam,
			      double _min, double _max)
{
  if (crit.size() != 1) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in Angle_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  for (int i=nin;i<n;i++) {
    if ( (crit[0].Includes(fl[i])) || ((crit[0].bar()).Includes(fl[i]) ) ) {
      cosmin[i][beam] = cosmin[beam][i] = AMATOOLS::Max(_min,-1.); 
      cosmax[i][beam] = cosmax[beam][i] = AMATOOLS::Min(_max, 1.); 
      AORGTOOLS::msg.Debugging()<<"Set cos-Range for "<<fl[i]<<" : "
				<<cosmin[beam][i]<<" ... "<<cosmax[beam][i]<<endl;
    }
  }
}





Mass_Selector::Mass_Selector(int _nin,int _nout, Flavour * _fl) {
  name = std::string("Mass_Selector"); 
  nin  = _nin; nout = _nout; n = nin+nout;
  fl   = _fl;
  
  massmin = new double*[n];
  massmax = new double*[n];
  value   = new double[n*n];

  for (int i=0;i<n;i++) { massmin[i] = new double[n]; massmax[i] = new double[n];}

  for (int i=nin;i<n;i++) {
    for (int j=i+1;j<n;j++) {
      massmin[i][j] = massmin[j][i] = 0.; 
      massmax[i][j] = massmax[j][i] = sqr(AORGTOOLS::rpa.gen.Ecms()); 
    }
  }
  sel_log = new Selector_Log(name);
}

Mass_Selector::~Mass_Selector() {
  //  delete [] sel_log;
  for (int i=0;i<n;i++) {
    delete [] massmin[i];
    delete [] massmax[i];
  }
  delete [] massmin;
  delete [] massmax;
  delete [] value;
}

bool Mass_Selector::Trigger(const vec4d * mom) 
{
  double massij;
  for (int i=nin;i<n;i++) {
    for (int j=i+1;j<n;j++) {
      massij = value[i*n+j] = (mom[i]+mom[j]).abs2();
      if (sel_log->Hit( ((massij < massmin[i][j]) || 
			 (massij > massmax[i][j])) )) return 0;
    }
  }
  return 1;
}

double * Mass_Selector::ActualValue() { return value; }

void Mass_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<n-1;i++) {
    for (int j=i+1;j<n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],massmin[i][j]);
      if (cuts->scut[i][j] < 
	  cuts->energymax[i]*cuts->energymax[j]*(1.-cuts->cosmax[i][j])) {
	cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	  Min(cuts->cosmax[i][j],1.-cuts->scut[i][j]/(cuts->energymax[i]*cuts->energymax[j]));
      }
    }
  }
}

void Mass_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=0;i<n-1;i++) {
    for (int j=i+1;j<n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],massmin[i][j]);
      if (cuts->scut[i][j] < 
	  cuts->energymax[i]*cuts->energymax[j]*(1.-cuts->cosmax[i][j])) {
	cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	  Min(cuts->cosmax[i][j],1.-cuts->scut[i][j]/(cuts->energymax[i]*cuts->energymax[j]));
      }
    }
  }
}


void Mass_Selector::SetRange(std::vector<APHYTOOLS::Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in Mass_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  for (int i=nin;i<n;i++) {
    for (int j=nin+1;i<n;i++) {
      if ( ((crit[0].Includes(fl[i])) && (crit[1].Includes(fl[j])) ) || 
	   ((crit[0].Includes(fl[j])) && (crit[1].Includes(fl[i])) ) ) {
	massmin[i][j] = massmin[j][i] = AMATOOLS::Max(_min,sqr(fl[i].mass()+fl[j].mass())); 
	massmax[i][j] = massmax[j][i] = AMATOOLS::Min(_max,sqr(AORGTOOLS::rpa.gen.Ecms())); 
	AORGTOOLS::msg.Debugging()<<"Set mass-Range for "<<fl[i]<<"/"<<fl[j]<<" : "
				  <<massmin[i][j]<<" ... "<<massmax[i][j]<<endl;
      }
    }
  }
}



/*--------------------------------------------------------------------

  Summed PT Selector

  --------------------------------------------------------------------*/

Summed_PT_Selector::Summed_PT_Selector(int _nin,int _nout, Flavour * _fl) {
  name = std::string("Summed_PT_Selector"); 
  nin  = _nin; nout = _nout; n = nin+nout;
  fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms();
  ptmin    = 0.; ptmax    = E;
  sel_log = new Selector_Log(name);
}

Summed_PT_Selector::~Summed_PT_Selector() { }


bool Summed_PT_Selector::Trigger(const vec4d * mom) 
{
  double pt;
  for (int i=nin;i<n;i++) {
    if (crit.Includes(fl[i])) {
      pt += sqrt(sqr(mom[i][1]) + sqr(mom[i][2]));
    }
  }
  if (sel_log->Hit( ((pt<ptmin) || (pt>ptmax)) )) return 0;
  return 1;
}

double * Summed_PT_Selector::ActualValue() { return NULL; }

void Summed_PT_Selector::BuildCuts(Cut_Data * cuts) { }

void Summed_PT_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) { }
 
void Summed_PT_Selector::SetRange(std::vector<APHYTOOLS::Flavour> _crit,double _min, 
			       double _max=0.5*AORGTOOLS::rpa.gen.Ecms())
{
  if (_crit.size() != 1) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in Summed_PT_Selector::SetRange : "
			  <<_crit.size()<<endl;
    return;
  }
  crit  = _crit[0];
  ptmin = _min; 
  ptmax = _max;
}




