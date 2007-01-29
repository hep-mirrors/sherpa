#include "Standard_Selector.H"
#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Flavour.H"

using namespace ATOOLS;
using namespace std;

class Order_Y {
public:
  bool operator()(const Vec4D &a,const Vec4D &b)
  {
    return a.Y()>b.Y();
  }
};

/*--------------------------------------------------------------------

  Energy Selector

  --------------------------------------------------------------------*/

Energy_Selector::Energy_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Energy_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  double E = rpa.gen.Ecms();
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
  delete m_sel_log;
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

void Energy_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n;i++) {
    cuts->energymin[i] = Max(emin[i],cuts->energymin[i]);
    cuts->energymax[i] = Min(emax[i],cuts->energymax[i]);
  }
}

void Energy_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
}
 
void Energy_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max=0.5*rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    msg.Error()<<"Wrong number of arguments in Energy_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  double MaxEmin = 0.;
  for (int i=m_nin;i<m_n;i++) {
    //    if (crit[0].Includes(m_fl[i]) || crit[0].Bar().Includes(m_fl[i])) {
    if (crit[0].Includes(m_fl[i])) {
      emin[i] = Max(_min,m_fl[i].SelMass()); 
      emax[i] = Min(_max,rpa.gen.Ecms());
      if (emin[i]>MaxEmin ) MaxEmin = emin[i];
    }
  }
  m_smin = Max(MaxEmin*MaxEmin,m_smin);
}



/*--------------------------------------------------------------------

  Transverse Energy Selector

  --------------------------------------------------------------------*/

ET_Selector::ET_Selector(int _nin,int _nout, Flavour * _fl) 
{
  m_name = std::string("ET_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  double E = rpa.gen.Ecms();
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
  delete m_sel_log;
}

bool ET_Selector::Trigger(const Vec4D * mom) 
{
  double eti;
  for (int i=m_nin;i<m_n;i++) {   
    eti = value[i] = mom[i][0]*mom[i].PPerp()/mom[i].P();
    if (m_sel_log->Hit( ((eti<etmin[i]) || (eti>etmax[i])) )) return 0;
  }
  return 1;
}

void ET_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n;i++) {
    cuts->energymin[i] = Max(etmin[i],cuts->energymin[i]);
    cuts->cosmax[0][i] = cuts->cosmax[1][i] = cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
      Min(cuts->cosmax[0][i],sqrt(1.-4.*sqr(etmin[i])/m_smax));
    cuts->etmin[i] = Max(etmin[i],cuts->etmin[i]);
  }
}

void ET_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{ }
 
void ET_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max=rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    msg.Error()<<"Wrong number of arguments in ET_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  double MaxEtmin = 0.;
  for (int i=m_nin;i<m_n;i++) {
    //    if (crit[0].Includes(m_fl[i]) || crit[0].Bar().Includes(m_fl[i])) {
    if (crit[0].Includes(m_fl[i])) {
      etmin[i] = _min; 
      etmax[i] = _max;
      if (etmin[i] > MaxEtmin) MaxEtmin = etmin[i];
    }
  }
  m_smin = Max(MaxEtmin*MaxEtmin,m_smin);
}


/*--------------------------------------------------------------------

  PT Selector

  --------------------------------------------------------------------*/

PT_Selector::PT_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("PT_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  double E = rpa.gen.Ecms();
  ptmin  = new double[m_n];
  ptmax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) { ptmin[i] = 0.; ptmax[i] = 10.*E; }
  m_sel_log = new Selector_Log(m_name);
}

PT_Selector::~PT_Selector() {
  delete [] ptmin;
  delete [] ptmax;
  delete [] value;
  delete m_sel_log;
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

void PT_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n;i++) {
    cuts->energymin[i] = Max(sqrt(sqr(ptmin[i])+sqr(m_fl[i].SelMass())),cuts->energymin[i]);
    cuts->cosmax[0][i] = cuts->cosmax[1][i] = cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
      Min(cuts->cosmax[0][i],sqrt(1.-sqr(ptmin[i])/(0.25*m_smax-sqr(m_fl[i].SelMass()))));
    cuts->etmin[i] = Max(sqrt(sqr(ptmin[i])+sqr(m_fl[i].SelMass())*(1.-sqr(cuts->cosmax[0][i]))),cuts->etmin[i]);
  }
}

void PT_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {}
 
void PT_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max=0.5*rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    msg.Error()<<"Wrong number of arguments in PT_Selector::SetRange : "
	       <<crit.size()<<endl;
    return;
  }

  double MaxPTmin = 0.;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      ptmin[i] = _min; 
      ptmax[i] = Min(_max,rpa.gen.Ecms());
      if (ptmin[i]>MaxPTmin) MaxPTmin = ptmin[i];
    }
  }
  m_smin = Max(m_smin,4.*MaxPTmin*MaxPTmin);
}

/*--------------------------------------------------------------------

  PT Selector

  --------------------------------------------------------------------*/

BFKL_PT_Selector::BFKL_PT_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("BFKL_PT_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  double E = rpa.gen.Ecms();
  ptmin  = new double[m_n];
  ptmax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) { ptmin[i] = 0.; ptmax[i] = 10.*E; }
  m_sel_log = new Selector_Log(m_name);
}

BFKL_PT_Selector::~BFKL_PT_Selector() {
  delete [] ptmin;
  delete [] ptmax;
  delete [] value;
  delete m_sel_log;
}

bool BFKL_PT_Selector::GetValue(const std::string &name,double &value)
{
  if (name=="qcut") value=ptmin[0];
  else return false;
  return true;
}

bool BFKL_PT_Selector::Trigger(const Vec4D * mom) 
{
  std::vector<Vec4D> moms(m_nout);
  for (size_t i(0);i<m_nout;++i) moms[i]=mom[m_nin+i];
  std::sort(moms.begin(),moms.end(),Order_Y());
  double pti(0.0);
  Vec4D qt(mom[0]);
  for (int i=m_nin;i<m_n;i++) {
    pti = value[i] = moms[i-m_nin].PPerp();
    if (m_sel_log->Hit( ((pti<ptmin[i]) || (pti>ptmax[i])) )) return 0;
    if (i<m_n-1) {
      qt-=moms[i-m_nin];
      pti = value[i] = qt.PPerp();
      if (m_sel_log->Hit( ((pti<ptmin[i]) || (pti>ptmax[i])) )) return 0;
    }
  }
  return 1;
}

void BFKL_PT_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void BFKL_PT_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {}
 
void BFKL_PT_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max=0.5*rpa.gen.Ecms())
{
  if (crit.size() != 1) {
    msg.Error()<<"Wrong number of arguments in BFKL_PT_Selector::SetRange : "
	       <<crit.size()<<endl;
    return;
  }

  double MaxPTmin = 0.;
  for (int i=m_nin;i<m_n;i++) {
    if (crit[0].Includes(m_fl[i])) {
      ptmin[i] = _min; 
      ptmax[i] = Min(_max,rpa.gen.Ecms());
      if (ptmin[i]>MaxPTmin) MaxPTmin = ptmin[i];
    }
  }
  m_smin = Max(m_smin,4.*MaxPTmin*MaxPTmin);
}

/*--------------------------------------------------------------------

  Z Selector

  --------------------------------------------------------------------*/

X_Selector::X_Selector(int nin,int nout,Flavour *fl) 
{
  m_name="X_Selector"; 
  m_nin=nin; 
  m_nout=nout; 
  m_n=m_nin+m_nout;
  m_fl=fl;
  m_smax=sqr(rpa.gen.Ecms());
  m_smin=0.0;
  zmin = new double[m_nin];
  zmax = new double[m_nin];
  value = new double[m_nin];
  for (int i=0;i<m_nin;i++) { 
    zmin[i]=0.0; 
    zmax[i]=1.0; 
  }
  m_sel_log = new Selector_Log(m_name);
}

X_Selector::~X_Selector() 
{
  delete [] zmin;
  delete [] zmax;
  delete [] value;
  delete m_sel_log;
}


bool X_Selector::Trigger(const Vec4D * mom) 
{
  for (int i=0;i<m_nin;i++) {
    if (i==0) value[i]=mom[i].PPlus()/rpa.gen.PBeam(i).PPlus();
    else value[i]=mom[i].PMinus()/rpa.gen.PBeam(i).PMinus();
    if (m_sel_log->Hit(value[i]<zmin[i] || 
		       value[i]>zmax[i])) return false;
  }
  return true;
}

void X_Selector::BuildCuts(Cut_Data *cuts) 
{
  cuts->energymin[0]=Max(cuts->energymin[0],
			 rpa.gen.PBeam(0).PPlus()*zmin[0]/2.0);
  cuts->energymin[1]=Max(cuts->energymin[1],
			 rpa.gen.PBeam(1).PMinus()*zmin[1]/2.0);
}

void X_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{
}
 
void X_Selector::SetRange(std::vector<Flavour> crit,double min,double max)
{
  for (int i=0;i<m_nin;i++) {
    zmin[i]=Max(min,1.0e-37); 
    zmax[i]=Min(max,1.0);
  }
  m_smin=Max(m_smin,min*rpa.gen.PBeam(0).PPlus()*rpa.gen.PBeam(1).PMinus());
  m_smax=Min(m_smax,max*rpa.gen.PBeam(0).PPlus()*rpa.gen.PBeam(1).PMinus());
}

/*--------------------------------------------------------------------

  Rapidity Selector

  --------------------------------------------------------------------*/


Rapidity_Selector::Rapidity_Selector(int _nin,int _nout, Flavour * _fl) {

  m_name = std::string("Rapidity_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  double E = rpa.gen.Ecms();
  double pl;

  ymin  = new double[m_n];
  ymax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) {
    pl      = sqrt(E*E-sqr(_fl[i].SelMass()));
    ymax[i] = log( (E+pl)/(E-pl) );
    ymin[i] = -ymax[i];
  }
  m_sel_log = new Selector_Log(m_name);
}

Rapidity_Selector::~Rapidity_Selector() {
  delete [] ymin;
  delete [] ymax;
  delete [] value;
  delete m_sel_log;
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

void Rapidity_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n;i++) {
    cuts->cosmax[0][i] = cuts->cosmax[i][0] =  
      Min(cuts->cosmax[0][i],1./sqrt(1.-sqr(m_fl[i].SelMass())/sqr(cuts->energymin[i]))*tanh(ymax[i]));
    cuts->cosmax[1][i] = cuts->cosmax[i][1] = 
      Min(cuts->cosmax[0][i],1./sqrt(1.-sqr(m_fl[i].SelMass())/sqr(cuts->energymin[i]))*tanh(-ymin[i]));
  }
}

void Rapidity_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n-1;i++) {
    if (!ATOOLS::IsZero(m_fl[i].SelMass())) {
      cuts->cosmax[0][i] = cuts->cosmax[i][0] =  
	1./sqrt(1.-sqr(m_fl[i].SelMass())/sqr(cuts->energymin[i]))*tanh(ymax[i]);
      cuts->cosmax[1][i] = cuts->cosmax[i][1] = 
	1./sqrt(1.-sqr(m_fl[i].SelMass())/sqr(cuts->energymin[i]))*tanh(-ymin[i]);
    }  
  }
}
 


void Rapidity_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max)
{
  if (crit.size() != 1) {
    msg.Error()<<"Wrong number of arguments in Rapidity_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  double E = rpa.gen.Ecms()/2;
  double pl,y;

  for (int i=m_nin;i<m_n;i++) {
    //    if (crit[0].Includes(m_fl[i]) || crit[0].Bar().Includes(m_fl[i])) {
    if (crit[0].Includes(m_fl[i])) {
      pl      = sqrt(E*E-sqr(m_fl[i].SelMass())); 
      y       = log((E+pl)/(E-pl));
      ymin[i] = Max(_min,-y);
      ymax[i] = Min(_max,y);
    }
  }
}

/*--------------------------------------------------------------------

  PseudoRapidity Selector

  --------------------------------------------------------------------*/


PseudoRapidity_Selector::PseudoRapidity_Selector(int _nin,int _nout, Flavour * _fl) {

  m_name = std::string("PseudoRapidity_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
 
  etamin  = new double[m_n];
  etamax  = new double[m_n];
  value = new double[m_n];
  for (int i=0;i<m_n;i++) {
      etamax[i] = 100.;
      etamin[i] = -etamax[i];
  }
  m_sel_log = new Selector_Log(m_name);
}

PseudoRapidity_Selector::~PseudoRapidity_Selector() {
  delete [] etamin;
  delete [] etamax;
  delete [] value;
  delete m_sel_log;
}


bool PseudoRapidity_Selector::Trigger(const Vec4D * mom) 
{
  double etai,theta;
  
  for (int i=m_nin;i<m_n;i++) {
    theta = acos(Vec3D(mom[i])*Vec3D(mom[0])/(Vec3D(mom[i]).Abs()*Vec3D(mom[0]).Abs())); 
    etai  = value[i] = -log(tan(theta/2.));
    if (m_sel_log->Hit( ((etai<etamin[i]) || (etai>etamax[i])) )) return 0;
  }
  return 1;
}

void PseudoRapidity_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n;i++) {
    cuts->cosmax[0][i] = cuts->cosmax[i][0] = Min(cuts->cosmax[0][i],tanh(etamax[i]));
    cuts->cosmax[1][i] = cuts->cosmax[i][1] = Min(cuts->cosmax[0][i],tanh(-etamin[i]));
  }
}

void PseudoRapidity_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {}
 


void PseudoRapidity_Selector::SetRange(std::vector<Flavour> crit,double _min, 
			       double _max)
{
  if (crit.size() != 1) {
    msg.Error()<<"Wrong number of arguments in PseudoRapidity_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  
  for (int i=m_nin;i<m_n;i++) {
    //    if ( (crit[0].Includes(m_fl[i])) || ((crit[0].Bar()).Includes(m_fl[i])) ) {
    if (crit[0].Includes(m_fl[i])) {
      etamin[i] = _min;
      etamax[i] = _max;
    }
  }
}


/*--------------------------------------------------------------------

  Angle Selector

  --------------------------------------------------------------------*/

Angle_Selector::Angle_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Angle_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();

  cosmin = new double*[m_n];
  cosmax = new double*[m_n];
  value  = new double[m_n*m_n];
  for (int i=0;i<m_n;i++) { cosmin[i] = new double[m_n]; cosmax[i] = new double[m_n]; }
  for (int i=0;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      //for numerical reason min and max should be slightly larger than -/+1
      cosmin[i][j] = cosmin[j][i] = -1.1; 
      cosmax[i][j] = cosmax[j][i] =  1.1; 
    }
  }
    
  m_sel_log = new Selector_Log(m_name);
}

Angle_Selector::~Angle_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] cosmin[i];
    delete [] cosmax[i];
  }
  delete [] cosmin;
  delete [] cosmax;
  delete [] value;
  delete m_sel_log;
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

void Angle_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=0;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
//       cuts->cosmin[i][j] = cuts->cosmin[j][i] = 
// 	Max(cosmin[i][j],cuts->cosmin[i][j]);
      cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	Min(cosmax[i][j],cuts->cosmax[i][j]);
    }
    if (i==1) {
      for (int j=i+1;j<m_n;j++) {
	cuts->cosmin[i][j] = cuts->cosmin[j][i] = 
	  Max(cuts->cosmin[i][j],-cuts->cosmax[0][j]);
	cuts->cosmax[i][j] = cuts->cosmax[j][i] = 
	  Min(cuts->cosmax[i][j],-cuts->cosmin[0][j]);
      }
    }
  }
}

void Angle_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {}


void Angle_Selector::SetRange(std::vector<Flavour> crit,
			      double _min, double _max)
{
  if (crit.size() != 2) {
    msg.Error()<<"Wrong number of arguments in Angle_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  //for numerical reasons exact +1 or -1 may be prolematic as borders
  if (IsEqual(_min,-1.)) _min = -1.1;
  if (IsEqual(_max,1.))  _max = 1.1;


  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	cosmin[i][j] = cosmin[j][i] = _min; 
	cosmax[i][j] = cosmax[j][i] = _max; 
      }
    }
  }
}

void Angle_Selector::SetRange(std::vector<Flavour> crit,int beam,
			      double _min, double _max)
{
  if (crit.size() != 1) {
    msg.Error()<<"Wrong number of arguments in Angle_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }
  for (int i=m_nin;i<m_n;i++) {
    if ( (crit[0].Includes(m_fl[i])) || ((crit[0].Bar()).Includes(m_fl[i]) ) ) {
      cosmin[i][beam] = cosmin[beam][i] = Max(_min,-1.1); 
      cosmax[i][beam] = cosmax[beam][i] = Min(_max, 1.1); 
    }
  }
}

/*--------------------------------------------------------------------

  Invariant Mass Selector

  --------------------------------------------------------------------*/

Mass_Selector::Mass_Selector(int _nin,int _nout, Flavour * _fl) {
  m_name = std::string("Mass_Selector"); 
  m_nin  = _nin; m_nout = _nout; m_n = m_nin+m_nout;
  m_fl   = _fl;
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
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
      massmax[i][j] = massmax[j][i] = 2.*rpa.gen.Ecms(); 
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

Mass_Selector::~Mass_Selector() {
  for (int i=0;i<m_n;i++) {
    delete [] massmin[i];
    delete [] massmax[i];
  }
  delete [] massmin;
  delete [] massmax;
  delete [] value;
  delete m_sel_log;
}

bool Mass_Selector::Trigger(const Vec4D * mom) 
{
  double massij;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      massij = value[i*m_n+j] = sqrt((mom[i]+mom[j]).Abs2());
      if (m_sel_log->Hit( ((massij < massmin[i][j]) || 
			 (massij > massmax[i][j])) )) return 0;
    }
  }
  return 1;
}

void Mass_Selector::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin;i<m_n-1;i++) {
    for (int j=i+1;j<m_n;j++) {
      cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],sqr(massmin[i][j]));
    }
  }
}

void Mass_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) {}


void Mass_Selector::SetRange(std::vector<Flavour> crit,double _min, double _max)
{
  if (crit.size() != 2) {
    msg.Error()<<"Wrong number of arguments in Mass_Selector::SetRange : "
			  <<crit.size()<<endl;
    return;
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=m_nin+1;j<m_n;j++) {
      if ( ((crit[0].Includes(m_fl[i])) && (crit[1].Includes(m_fl[j])) ) || 
	   ((crit[0].Includes(m_fl[j])) && (crit[1].Includes(m_fl[i])) ) ) {
	massmin[i][j] = massmin[j][i] = Max(_min,m_fl[i].SelMass()+m_fl[j].SelMass()); 
	massmax[i][j] = massmax[j][i] = _max;
	if (sqr(massmin[i][j])>m_smin) m_smin = Max(sqr(massmin[i][j]),m_smin);
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
  m_smin = 0.;
  m_smax = rpa.gen.Ecms()*rpa.gen.Ecms();
  
  double E = rpa.gen.Ecms();
  ptmin = new double*[m_n];
  ptmax = new double*[m_n];

  for (int i=0;i<m_n;i++) { 
    ptmin[i] = new double[m_n]; 
    ptmax[i] = new double[m_n];
  }

  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      ptmin[i][j] = ptmin[j][i] = 0.; 
      ptmax[i][j] = ptmax[j][i] = E; 
    }
  }
  m_sel_log = new Selector_Log(m_name);
}

Summed_PT_Selector::~Summed_PT_Selector() {  
  for (int i=0;i<m_n;i++) {
    delete [] ptmin[i];
    delete [] ptmax[i];
  }
  delete [] ptmin;
  delete [] ptmax;
  delete m_sel_log;
 }


bool Summed_PT_Selector::Trigger(const Vec4D * mom) 
{
  double ptij;
  for (int i=m_nin;i<m_n;i++) {
    for (int j=i+1;j<m_n;j++) {
      ptij = (mom[i]+mom[j]).PPerp();
      if (m_sel_log->Hit( ((ptij < ptmin[i][j]) || 
			   (ptij > ptmax[i][j])) )) return 0;
    }
  }
  return 1;
}

void Summed_PT_Selector::BuildCuts(Cut_Data * cuts) 
{
}

void Summed_PT_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts) { }
 
void Summed_PT_Selector::SetRange(std::vector<Flavour> _crit,double _min, 
			       double _max=0.5*rpa.gen.Ecms())
{
  if (_crit.size() != 2) {
    msg.Error()<<"Wrong number of arguments in Summed_PT_Selector::SetRange : "
			  <<_crit.size()<<endl;
    return;
  }
  for (int i=m_nin;i<m_n-1;i++) {
    for (int j=m_nin+1;j<m_n;j++) {
      if ( ((_crit[0].Includes(m_fl[i])) && (_crit[1].Includes(m_fl[j])) ) || 
	   ((_crit[0].Includes(m_fl[j])) && (_crit[1].Includes(m_fl[i])) ) ) {
	ptmin[i][j] = ptmin[j][i] = _min; 
	ptmax[i][j] = ptmax[j][i] = _max;
	if (sqr(ptmin[i][j])>m_smin) m_smin = Max(sqr(ptmin[i][j]),m_smin);
      }
    }
  }
}




