#include "ALPACA/Tools/Parton.H"
#include "ATOOLS/Org/Message.H"

#include <list>

namespace ALPACA {
  long unsigned int Parton::s_totalnumber=0;
}

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;

Parton::Parton() :
  m_number(-1),
  m_momentum(Vec4D(0.,0.,0.,0.)), m_oldmom(Vec4D(0.,0.,0.,0.)),
  m_initpos(Vec4D(0.,0.,0.,0.)),
  m_flav(Flavour(kf_none)), m_flow(Flow()),
  m_lambda(1.), //m_lambda(HIPars("lambda")),
  m_inittau(0.), p_split_partner(nullptr), p_split_recoil_partner(nullptr),
  p_last_scatter_1(nullptr), p_last_scatter_2(nullptr), p_last_scatter_3(nullptr),
  m_formation_tau(make_pair(-1.,-1.)), m_formation_time(make_pair(-1.,-1.)), m_splitted_merged(false),
  m_N_scatter(0), m_timekeeper(0), m_ahat(Vec4D(0.,0.,0.,0.))
{ 
    a_momentum[0] = 0.;
	a_momentum[1] = 0.;
	a_momentum[2] = 0.;
	a_momentum[3] = 0.;
	a_initpos[0] = 0.;
	a_initpos[1] = 0.;
	a_initpos[2] = 0.;
	a_initpos[3] = 0.;
}

Parton::Parton(Flavour flav, Vec4D mom) :
  m_number(-1),
  m_momentum(mom), m_oldmom(mom), m_initpos(Vec4D(0.,0.,0.,0.)),
  m_flav(flav), m_flow(Flow()),
    m_lambda(1.), //m_lambda(HIPars("lambda")),
  m_inittau(0.), p_split_partner(nullptr), p_split_recoil_partner(nullptr),
  p_last_scatter_1(nullptr), p_last_scatter_2(nullptr), p_last_scatter_3(nullptr),
  m_formation_tau(make_pair(-1.,-1.)), m_formation_time(make_pair(-1.,-1.)), m_splitted_merged(false),
  m_N_scatter(0), m_timekeeper(0), m_ahat(Vec4D(0.,0.,0.,0.))
{ 
  a_momentum[0] = mom[1];
	a_momentum[1] = mom[2];
	a_momentum[2] = mom[3];
	a_momentum[3] = mom[0];
	a_initpos[0] = 0.;
	a_initpos[1] = 0.;
	a_initpos[2] = 0.;
	a_initpos[3] = 0.;
}

Parton::Parton(Flavour flav, Vec4D mom, Vec4D pos, double tau) :
  m_number(-1),
  m_momentum(mom), m_oldmom(mom), m_initpos(pos),
  m_flav(flav), m_flow(Flow()),
    m_lambda(1.), //m_lambda(HIPars("lambda")),
  m_inittau(tau), p_split_partner(nullptr), p_split_recoil_partner(nullptr),
  p_last_scatter_1(nullptr), p_last_scatter_2(nullptr), p_last_scatter_3(nullptr),
  m_formation_tau(make_pair(-1.,-1.)), m_formation_time(make_pair(-1.,-1.)), m_splitted_merged(false),
  m_N_scatter(0), m_timekeeper(0), m_ahat(Vec4D(0.,0.,0.,0.))
{ 
  a_momentum[0] = mom[1];
	a_momentum[1] = mom[2];
	a_momentum[2] = mom[3];
	a_momentum[3] = mom[0];
	a_initpos[0] = pos[1];
	a_initpos[1] = pos[2];
	a_initpos[2] = pos[3];
	a_initpos[3] = pos[0];
}

Parton::Parton(Parton * part) :
  m_number(-1),
  m_momentum(part->m_momentum), m_oldmom(part->m_oldmom), m_initpos(part->m_initpos),
  m_flav(part->m_flav), m_flow(part->m_flow), p_split_recoil_partner(part->GetSplitRecoilPartner()),
    m_lambda(1.), //m_lambda(HIPars("lambda")),
  m_inittau(part->m_inittau), p_split_partner(part->GetSplitPartner()),
  p_last_scatter_1(part->GetLastScatter()), p_last_scatter_2(part->GetSecondLastScatter()), p_last_scatter_3(part->GetThirdLastScatter()),
  m_formation_tau(make_pair(-1.,-1.)), m_formation_time(make_pair(-1.,-1.)), m_splitted_merged(false),
  m_N_scatter(part->GetNScatter()), m_timekeeper(0), m_ahat(Vec4D(0.,0.,0.,0.))
{ 
  	a_momentum[0] = m_momentum[1];
	a_momentum[1] = m_momentum[2];
	a_momentum[2] = m_momentum[3];
	a_momentum[3] = m_momentum[0];
	a_initpos[0] = m_initpos[1];
	a_initpos[1] = m_initpos[2];
	a_initpos[2] = m_initpos[3];
	a_initpos[3] = m_initpos[0];
}

Parton::~Parton() {
}

bool Parton::operator==(Parton part){
	if ((part.m_momentum == m_momentum) && 
		(part.m_initpos == m_initpos) &&
		(part.m_inittau == m_inittau) &&
		(part.m_flav == m_flav)) {
		return true; 
	}
	return false;
}

void Parton::SetNumber()           
{ 
  	m_number=++s_totalnumber;
}

Vec4D Parton::Position(const double tau) const{
	Vec4D pos;
	pos = 2.*m_lambda*m_momentum*(tau-m_inittau) + m_initpos;
	return pos;
}

Vec4D Parton::PositionFromt(const double t) const{
	double tau = (t-m_initpos[0])/(2*m_lambda*m_momentum[0])+m_inittau;
	Vec4D pos;
	pos = 2.*m_lambda*m_momentum*(tau-m_inittau) + m_initpos;
	return pos;
}

double Parton::Tau(const double t) const{
  	return (t-m_initpos[0])/(2*m_lambda*m_momentum[0])+m_inittau; 
}

double Parton::Dist2(const shared_ptr<Parton> part, const double tau) {
	if ((*this)==(*part)) return 0.;
	Vec4D x(Position(tau) - part->Position(tau));
	Vec4D p(m_momentum + part->Momentum());
	return -(x.Abs2() - sqr(x[0]*p[0]-x[1]*p[1]-x[2]*p[2]-x[3]*p[3])/p.Abs2());
}

void Parton::SetLastScatter(std::shared_ptr<Parton> part){
	p_last_scatter_3 = p_last_scatter_2;
	p_last_scatter_2 = p_last_scatter_1;
	p_last_scatter_1 = part;
}

void Parton::SetFormationTime(double t, double formation_time){
	m_formation_time = {t, t + formation_time}; //Decide here how formation t should be set
	m_formation_tau = {Tau(m_formation_time.first), Tau(m_formation_time.second)};
}

void Parton::AddXPHistory(int process_type, ATOOLS::Vec4D x, ATOOLS::Vec4D p, ATOOLS::Flavour flav, double tau){
	m_x_p_history.push_back(std::make_pair(process_type, std::make_pair(std::make_pair(x, p), flav)));
	m_x_p_tau.push_back(tau);
}

double Parton::GetLastXPt(){
	std::pair<int, std::pair<momentum_save, ATOOLS::Flavour>> xp = m_x_p_history.back();
	Vec4D x = xp.second.first.first;
	return x[0];
}


pair<bool, pair<pair<Vec4D,Vec4D>, Flavour>> Parton::GetXP(double t){
	//Function to find the old position and momentum at t by looping through momentum changes
	bool exists = true;
	Flavour flav;
	int n = m_x_p_history.size();
	Vec4D x, p, xp, pp;

	if(t < m_x_p_history[0].second.first.first[0]){
		exists = false;
	} else if(t > m_x_p_history[n-1].second.first.first[0] && m_splitted_merged){
		exists = false;
	}
	
	if(exists){
		exists = false;
		//for(int i = 0; i < n-1; i++){
		for (vector<pair<int,pair<momentum_save, Flavour>>>::reverse_iterator iter = m_x_p_history.rbegin(); iter != m_x_p_history.rend(); ++iter ) { 
			//if(t >= m_x_p_history[n-1-i].first[0]){
			if(t >= iter->second.first.first[0]){
				xp = iter->second.first.first;
				pp = iter->second.first.second;
				flav = iter->second.second;
				exists = true;
				break;
			}
		}

		if(!exists){
			msg_Out() << METHOD << ": ERROR: t = " << t << " not found, should not occur" << endl;
		}
	
		if(exists){
			x = pp*(t-xp[0])/pp[0] + xp;
			p = pp;
		}
	}

	return make_pair(exists, make_pair(make_pair(x,p), flav));
}

std::pair<bool, std::pair<std::pair<ATOOLS::Vec4D, ATOOLS::Vec4D>, ATOOLS::Flavour>> Parton::GetXPTau(double tau){
	//Function to find the old position and momentum at t by looping through momentum changes
	bool exists = true;
	Flavour flav;
	int n = m_x_p_history.size();
	Vec4D p, x;

	if(tau < m_x_p_tau[0]){
		exists = false;
	} else if(tau > m_x_p_tau[n-1] && m_splitted_merged){
		exists = false;
	}
	
	pair<int,pair<momentum_save, Flavour>> temp_xp;
	double temp_tau;
	if(exists){
		exists = false;
		for(int i = m_x_p_history.size()-1; i>=0; i--){
			temp_xp = m_x_p_history[i];
			temp_tau = m_x_p_tau[i];
			if(tau >= temp_tau){
				p = temp_xp.second.first.second;
				exists = true;
				break;
			}
		}

		if(!exists){
			msg_Out() << METHOD << ": ERROR: t = " << tau << " not found, should not occur" << endl;
		}
	}

	return make_pair(exists, make_pair(make_pair(x,p), flav));
}

std::array<double,4> Parton::getPosition(double t) {
	std::array<double,4> pos;
	pos[0] = a_momentum[0]*(t-a_initpos[3])/a_momentum[3]+a_initpos[0];
	pos[1] = a_momentum[1]*(t-a_initpos[3])/a_momentum[3]+a_initpos[1];
	pos[2] = a_momentum[2]*(t-a_initpos[3])/a_momentum[3]+a_initpos[2];
	pos[3] = t;
	for (size_t i = 0; i < 3; ++i) {
        int a = (int)pos[i];
		if (pos[i] > 0.) {
            pos[i] = pos[i];
		}
		else {
            pos[i] = pos[i];
		}
	}
	return pos;
}

std::array<double,4> Parton::getPositionNoBound(double t) {
	std::array<double,4> pos;
	pos[0] = a_momentum[0]*(t-a_initpos[3])/a_momentum[3]+a_initpos[0];
	pos[1] = a_momentum[1]*(t-a_initpos[3])/a_momentum[3]+a_initpos[1];
	pos[2] = a_momentum[2]*(t-a_initpos[3])/a_momentum[3]+a_initpos[2];
	pos[3] = t;
	return pos;
}

std::array<double,4> Parton::getPositionTau(double tau) {
	std::array<double,4> pos;
	//if (m_mass == 0.) {
		pos[0] = 2.*m_lambda*a_momentum[0]*(tau-m_inittau) + a_initpos[0];
		pos[1] = 2.*m_lambda*a_momentum[1]*(tau-m_inittau) + a_initpos[1];
		pos[2] = 2.*m_lambda*a_momentum[2]*(tau-m_inittau) + a_initpos[2];
		pos[3] = 2.*m_lambda*a_momentum[3]*(tau-m_inittau) + a_initpos[3];
	//}
	//else {
	//	pos[0] = a_momentum[0]*(tau-m_inittau)/m_mass + a_initpos[0];
	//	pos[1] = a_momentum[1]*(tau-m_inittau)/m_mass + a_initpos[1];
	//	pos[2] = a_momentum[2]*(tau-m_inittau)/m_mass + a_initpos[2];
	//	pos[3] = a_momentum[3]*(tau-m_inittau)/m_mass + a_initpos[3];
	//}
	for (size_t i = 0; i < 3; ++i) {
        int a = (int)pos[i];
		if (pos[i] > 0.) {
            pos[i] = pos[i] - a*1.;
		}
		else {
            pos[i] = pos[i] - (a-1)*1.;
		}
	}
	return pos;
}

std::array<double,4> Parton::getPositionTauNoBound(double tau) {
	std::array<double,4> pos;
	//if (m_mass == 0.) {
		pos[0] = 2.*m_lambda*a_momentum[0]*(tau-m_inittau) + a_initpos[0];
		pos[1] = 2.*m_lambda*a_momentum[1]*(tau-m_inittau) + a_initpos[1];
		pos[2] = 2.*m_lambda*a_momentum[2]*(tau-m_inittau) + a_initpos[2];
		pos[3] = 2.*m_lambda*a_momentum[3]*(tau-m_inittau) + a_initpos[3];
	//}
	//else {
	//	pos[0] = a_momentum[0]*(tau-m_inittau)/m_mass + a_initpos[0];
	//	pos[1] = a_momentum[1]*(tau-m_inittau)/m_mass + a_initpos[1];
	//	pos[2] = a_momentum[2]*(tau-m_inittau)/m_mass + a_initpos[2];
	//	pos[3] = a_momentum[3]*(tau-m_inittau)/m_mass + a_initpos[3];
	//}
	return pos;
}

double Parton::TauBar(std::shared_ptr<Parton> part) {
  std::array<double,3> shift = {0., 0., 0.};
  return TauBar(part, shift);	
}

double Parton::TauBar(std::shared_ptr<Parton> part, std::array<double,3> shift) {
  if ((*this) == (*part)){
	return m_inittau;
  } 
  double taubar;

  // TEMPORARY, box method will not work with for timekeeper = 1 or 2
  if(m_timekeeper == 0){
	std::array<double,4> shift4 = {shift[0], shift[1], shift[2], 0.};
	std::array<double,4> ipos(subtract(a_initpos, add(part->aPosition(),shift4)));
	std::array<double,4> dmom(subtract(a_momentum, part->aMomentum()));
	std::array<double,4> pi(a_momentum);
	std::array<double,4> pj(part->aMomentum());
	double iposdmom(dotpr(ipos,dmom));
	double pipj(dotpr(pi,pj));
	double taubar;
	//if (m_mass == 0.) {
		taubar =  ((m_inittau + part->Inittau())/2. + (iposdmom)/(4.*m_lambda*pipj));
	//}
	//else {
		//  taubar = (m_inittau + part->Inittau())/2. + m_mass*iposdmom/(2.*(pipj-sqr(m_mass)));
	//}
  } else{
	Vec4D shift_vec(0., shift[0], shift[1], shift[2]);
	double lambda_1 = m_lambda;
	double lambda_2 = part->GetLambda();
	Vec4D p_1 = m_momentum;
	Vec4D p_2 = part->Momentum();
	Vec4D P = p_1 + p_2;
	Vec4D v = 2.*lambda_1*p_1 - 2.*lambda_2*p_2;
	Vec4D v_T = v - (v*P)*P/(P*P);

	double init_tau_1 = m_inittau;
	double init_tau_2 = part->Inittau();
	
	double taubar_1 = init_tau_1 - (m_initpos -  (part->Position(init_tau_1)+shift_vec))*v_T/(v_T*v_T);
	double taubar_2 = init_tau_2 - (Position(init_tau_2) - (part->Position()+shift_vec))*v_T/(v_T*v_T);
	double taubar_3 = (init_tau_1+init_tau_2)/2. - (m_initpos - (part->Position()+shift_vec))*v_T/(v_T*v_T);
	
	if(abs(taubar_1/taubar_2) > 1.001 || abs(taubar_1/taubar_2) < 0.999){
		msg_Out() << METHOD << ": WARNING" "taubar_1 = " << taubar_1 << ", taubar_2 = " << taubar_2 << ", abs(taubar_1-taubar_2) > 0.00000001" << endl;
	}

	taubar = taubar_1;
  }
  
  return taubar;
}

double Parton::getInvDistTau(std::shared_ptr<Parton> part, double tau, std::array<double,3> shift) {
  if ((*this) == (*part)) return 0.;
  std::array<double,4> shift4 = {shift[0], shift[1], shift[2], 0.};
  std::array<double,4> x(subtract(getPositionTauNoBound(tau), add(part->getPositionTauNoBound(tau),shift4)));
  std::array<double,4> p(add(a_momentum, part->aMomentum()));
  //double dij2 = -(abs2(x) - sqr(dotpr(x,p))/abs2(p));
  double dij2 = -(abs2(x) - sqr(dotpr(x,p)/abs_vec(p)));
  //double dij2 = - dotpr(subtract(x,scaltimes(dotpr(x,p)/dotpr(p,p),p)),subtract(x,scaltimes(dotpr(x,p)/dotpr(p,p),p)));
  //if (dij2<=0.) std::cout<<"d_ij^2(tau) = "<<dij2<<"  "<<dotpr(x,p)<<std::endl;
  return sqrt(max(dij2,0.));
}



double Parton::Abs3Momentum(){
    return sqrt(m_momentum[1]*m_momentum[1]+m_momentum[2]*m_momentum[2]+m_momentum[3]*m_momentum[3]);
}


double Parton::TauBound(size_t i, double bound) {
    double pos,mom,taub;
    pos = m_initpos[i];
    mom = m_momentum[i];
    taub = (bound-pos)/(2.*m_lambda*mom)+m_inittau;
    return taub;
}

ostream & ALPACA::
operator<<(ostream & s, const Parton & part) {
  s<<"   "<<part.Flav()<<"  "<<part.Momentum()<<" "
   <<"{"<<part.GetFlow(1)<<" "<<part.GetFlow(2)<<"}"
   <<" at "<<part.Position()<<"(for tau = "<<part.Inittau()<<").\n";
  return s;
}
