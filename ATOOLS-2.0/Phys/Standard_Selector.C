#include "Standard_Selector.H"
#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Flavour.H"

using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

/*--------------------------------------------------------------------

  Energy Selector

  --------------------------------------------------------------------*/

Energy_Selector::Energy_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Energy_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms();
  emin  = new double[m_n];
  emax  = new double[m_n];
  value  = new double[m_n];
  for (int i=0;i<m_n;i++) { emin[i] = 0.; emax[i] = E; }
  m_sel_log = new Selector_Log(m_name);
}

Energy_Selector::~Energy_Selector() 
{
  delete [] emin;
  delete [] emax;
  delete [] value;
}


bool Energy_Selector::Trigger(const Vec4D * mom) 
{
  double ei;
  for (int i=m_nin;i<m_n;i++) {
    ei = value[i] = mom[i][0];
    if (m_sel_log->Hit( ((ei<emin[i]) || (ei>emax[i])) )) return 0;
  }
  return 1;
}

double * Energy_Selector::ActualValue() { return value; }

void Energy_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n-1;i++) {
    cuts->energymin[i] = Max(emin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(emax[i],cuts->energymax[i]);
    for (int j=i+1;j<m_n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}

void Energy_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=0;i<m_n-1;i++) {
    cuts->energymin[i] = Max(emin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(emax[i],cuts->energymax[i]);
    for (int j=i+1;j<m_n;j++) {
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

  for (int i=m_nin;i<m_n;i++) {
    if ( (crit[0].Includes(m_fl[i])) || ((crit[0].Bar()).Includes(m_fl[i]) ) ) {
      emin[i] = AMATOOLS::Max(_min,m_fl[i].Mass()); 
      emax[i] = AMATOOLS::Min(_max,0.5*AORGTOOLS::rpa.gen.Ecms());
      AORGTOOLS::msg.Debugging()<<"Set e-Range for "<<m_fl[i]<<" : "
				<<emin[i]<<" ... "<<emax[i]<<endl;
    }
  }
}



/*--------------------------------------------------------------------

  Transverse Energy Selector

  --------------------------------------------------------------------*/

ET_Selector::ET_Selector(int _nin,int _nout, Flavour * _fl) 
{
  m_name = std::string("ET_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms();
  etmin  = new double[m_n];
  etmax  = new double[m_n];
  value  = new double[m_n];
  for (int i=0;i<m_n;i++) { etmin[i] = 0.; etmax[i] = E; }
  m_sel_log = new Selector_Log(m_name);
}

ET_Selector::~ET_Selector() 
{
  delete [] etmin;
  delete [] etmax;
  delete [] value;
}

bool ET_Selector::Trigger(const Vec4D * mom) 
{
  double eti;
  for (int i=m_nin;i<m_n;i++) {
   
    double theta = acos(Vec3D(mom[i])*Vec3D(mom[0])/(Vec3D(mom[i]).Abs()*Vec3D(mom[0]).Abs())); 

    //cout<<"theta "<<theta<<endl;
    
    //if (theta<0.) theta = -theta; 

    eti = value[i] = mom[i][0]*sin(theta);

    //cout<<"et["<<i<<"] : "<<eti<<" range ("<<etmin[i]<<", "<<etmax[i]<<")"<<endl;

    if (m_sel_log->Hit( ((eti<etmin[i]) || (eti>etmax[i])) )) return 0;
  }
  return 1;
}

double * ET_Selector::ActualValue() { return value; }

void ET_Selector::BuildCuts(Cut_Data * cuts) 
{
  
  for (int i=0;i<m_n-1;i++) {
    cuts->energymin[i] = Max(etmin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(etmax[i],cuts->energymax[i]);
  }
}

void ET_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{
  
  for (int i=0;i<m_n-1;i++) {
    cuts->energymin[i] = Max(etmin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(etmax[i],cuts->energymax[i]);
  }
}
 
void ET_Selector::SetRange(std::vector<APHYTOOLS::Flavour> crit,double _min, 
			       double _max=0.5*AORGTOOLS::rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    AORGTOOLS::msg.Error()<<"Wrong number of arguments in ET_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=m_nin;i<m_n;i++) {
    if ( (crit[0].Includes(m_fl[i])) || ((crit[0].Bar()).Includes(m_fl[i]) ) ) {
      etmin[i] = _min; 
      etmax[i] = AMATOOLS::Min(_max,0.5*AORGTOOLS::rpa.gen.Ecms());
      AORGTOOLS::msg.Debugging()<<"Set et-Range for "<<m_fl[i]<<" : "
				<<etmin[i]<<" ... "<<etmax[i]<<endl;
    }
  }
}


/*--------------------------------------------------------------------

  PT Selector

  --------------------------------------------------------------------*/

PT_Selector::PT_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("PT_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms();
  ptmin  = new double[m_n];
  ptmax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) { ptmin[i] = 0.; ptmax[i] = E; }
  m_sel_log = new Selector_Log(m_name);
}

PT_Selector::~PT_Selector() {
  //  delete [] m_sel_log;
  delete [] ptmin;
  delete [] ptmax;
  delete [] value;
}


bool PT_Selector::Trigger(const Vec4D * mom) 
{
  double pti;
  for (int i=m_nin;i<m_n;i++) {
    pti = value[i] = sqrt(sqr(mom[i][1]) + sqr(mom[i][2]));
    if (m_sel_log->Hit( ((pti<ptmin[i]) || (pti>ptmax[i])) )) return 0;
  }
  return 1;
}

double * PT_Selector::ActualValue() { return value; }

void PT_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n-1;i++) {
    cuts->energymin[i] = Max(ptmin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(ptmax[i],cuts->energymax[i]);
    for (int j=i+1;j<m_n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = 
	Max(cuts->scut[i][j],2.*cuts->energymin[i]*cuts->energymin[j]*(1.-cuts->cosmax[i][j]));
    }
  }
}

void PT_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=0;i<m_n-1;i++) {
    cuts->energymin[i] = Max(ptmin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(ptmax[i],cuts->energymax[i]);
    for (int j=i+1;j<m_n;j++) {
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

  for (int i=m_nin;i<m_n;i++) {
    if ( (crit[0].Includes(m_fl[i])) || ((crit[0].Bar()).Includes(m_fl[i]) ) ) {
      ptmin[i] = AMATOOLS::Max(_min,m_fl[i].Mass()); 
      ptmax[i] = AMATOOLS::Min(_max,0.5*AORGTOOLS::rpa.gen.Ecms());
      AORGTOOLS::msg.Debugging()<<"Set PT-Range for "<<m_fl[i]<<" : "
				<<ptmin[i]<<" ... "<<ptmax[i]<<endl;
    }
  }
}




/*--------------------------------------------------------------------

  Rapidity Selector

  --------------------------------------------------------------------*/


Rapidity_Selector::Rapidity_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Rapidity_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms()/2;
  double pl;

  ymin  = new double[m_n];
  ymax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) {
    pl      = sqrt(E*E-sqr(_fl[i].Mass()));
    ymax[i] = log( (E+pl)/(E-pl) );
    ymin[i] = -ymax[i];
  }
  m_sel_log = new Selector_Log(m_name);
}

Rapidity_Selector::~Rapidity_Selector() {
  //  delete [] m_sel_log;
  delete [] ymin;
  delete [] ymax;
  delete [] value;
}


bool Rapidity_Selector::Trigger(const Vec4D * mom) 
{
  double yi;
  for (int i=m_nin;i<m_n;i++) {
    yi = value[i] = 0.5 * log( (mom[i][0]+mom[i][3])/(mom[i][0]-mom[i][3]) );
    if (m_sel_log->Hit( ((yi<ymin[i]) || (yi>ymax[i])) )) return 0;
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

  for (int i=m_nin;i<m_n;i++) {
    pl      = sqrt(E*E-sqr(m_fl[i].Mass())); 
    y       = log((E+pl)/(E-pl));
    ymin[i] = AMATOOLS::Max(_min,-y);
    ymax[i] = AMATOOLS::Min(_max,y);
    AORGTOOLS::msg.Debugging()<<"Set y-Range for "<<m_fl[i]<<" : "
			      <<ymin[i]<<" ... "<<ymax[i]<<endl;
  }
}





Angle_Selector::Angle_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Angle_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;

  cosmin = new double*[m_n];
  cosmax = new double*[m_n];
  value  = new double[m_n*m_n];
  for (int i=0;i<m_n;i++) { cosmin[i] = new double[m_n]; cosmax[i] = new double[m_n]; }
  for (int i=0;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      cosmin[i][j] = cosmin[j][i] = -1.; 
      cosmax[i][j] = cosmax[j][i] =  1.; 
    }
  }
    
  m_sel_log = new Selector_Log(m_name);
}

Angle_Selector::~Angle_Selector() {
  //  delete [] m_sel_log;
  for (int i=0;i<m_n;i++) {
    delete [] cosmin[i];
    delete [] cosmax[i];
  }
  delete [] cosmin;
  delete [] cosmax;
  delete [] value;
}

bool Angle_Selector::Trigger(const Vec4D * mom) 
{
  double cosij;
  for (int i=0;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      cosij = value[m_n*i+j] = 
	Vec3D(mom[i])*Vec3D(mom[j])/(Vec3D(mom[i]).Abs()*Vec3D(mom[j]).Abs());
      if (m_sel_log->Hit( ((cosij < cosmin[i][j]) || 
			 (cosij > cosmax[i][j])) )) return 0;
    }
  }
  return 1;
}

double * Angle_Selector::ActualValue() { return value; }

void Angle_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
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
  for (int i=0;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
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
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	cosmin[i][j] = cosmin[j][i] = AMATOOLS::Max(_min,-1.); 
	cosmax[i][j] = cosmax[j][i] = AMATOOLS::Min(_max,1.); 
	AORGTOOLS::msg.Debugging()<<"Set cos-Range for "<<m_fl[i]<<"/"<<m_fl[j]<<" : "
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
  for (int i=m_nin;i<m_n;i++) {
    if ( (crit[0].Includes(m_fl[i])) || ((crit[0].Bar()).Includes(m_fl[i]) ) ) {
      cosmin[i][beam] = cosmin[beam][i] = AMATOOLS::Max(_min,-1.); 
      cosmax[i][beam] = cosmax[beam][i] = AMATOOLS::Min(_max, 1.); 
      AORGTOOLS::msg.Debugging()<<"Set cos-Range for "<<m_fl[i]<<" : "
				<<cosmin[beam][i]<<" ... "<<cosmax[beam][i]<<endl;
    }
  }
}





Mass_Selector::Mass_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Mass_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  
  massmin = new double*[m_n];
  massmax = new double*[m_n];
  value   = new double[m_n*m_n];

  for (int i=0;i<m_n;i++) { 
    massmin[i] = new double[m_n]; 
    massmax[i] = new double[m_n];
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      massmin[i][j] = massmin[j][i] = 0.; 
      massmax[i][j] = massmax[j][i] = sqr(AORGTOOLS::rpa.gen.Ecms()); 
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

Mass_Selector::~Mass_Selector() {
  //  delete [] m_sel_log;
  for (int i=0;i<m_n;i++) {
    delete [] massmin[i];
    delete [] massmax[i];
  }
  delete [] massmin;
  delete [] massmax;
  delete [] value;
}

bool Mass_Selector::Trigger(const Vec4D * mom) 
{
  double massij;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      massij = value[i*m_n+j] = (mom[i]+mom[j]).Abs2();
      if (m_sel_log->Hit( ((massij < massmin[i][j]) || 
			 (massij > massmax[i][j])) )) return 0;
    }
  }
  return 1;
}

double * Mass_Selector::ActualValue() { return value; }

void Mass_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],massmin[i][j]);
      if (cuts->scut[i][j] < 
	  cuts->energymax[i]*cuts->energymax[j]*(1.-cuts->cosmax[i][j])) {
	cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	  Min(cuts->cosmax[i][j],1.-cuts->scut[i][j]/(cuts->energymax[i]*cuts->energymax[j]));
      }
    }
  }
}

void Mass_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{
  for (int i=0;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
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
  for (int i=m_nin;i<m_n;i++) {
    for (int j=m_nin+1;i<m_n;i++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	massmin[i][j] = massmin[j][i] = AMATOOLS::Max(_min,sqr(m_fl[i].Mass()+m_fl[j].Mass())); 
	massmax[i][j] = massmax[j][i] = AMATOOLS::Min(_max,sqr(AORGTOOLS::rpa.gen.Ecms())); 
	AORGTOOLS::msg.Debugging()<<"Set mass-Range for "<<m_fl[i]<<"/"<<m_fl[j]<<" : "
				  <<massmin[i][j]<<" ... "<<massmax[i][j]<<endl;
      }
    }
  }
}



/*--------------------------------------------------------------------

  Summed PT Selector

  --------------------------------------------------------------------*/

Summed_PT_Selector::Summed_PT_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Summed_PT_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  
  double E = AORGTOOLS::rpa.gen.Ecms();
  ptmin    = 0.; ptmax    = E;
  m_sel_log = new Selector_Log(m_name);
}

Summed_PT_Selector::~Summed_PT_Selector() { }


bool Summed_PT_Selector::Trigger(const Vec4D * mom) 
{
  double pt;
  for (int i=m_nin;i<m_n;i++) {
    if (crit.Includes(m_fl[i])) {
      pt += sqrt(sqr(mom[i][1]) + sqr(mom[i][2]));
    }
  }
  if (m_sel_log->Hit( ((pt<ptmin) || (pt>ptmax)) )) return 0;
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




