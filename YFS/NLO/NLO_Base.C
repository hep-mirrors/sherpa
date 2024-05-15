// #include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/NLO_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"


using namespace YFS;
using namespace MODEL;
using namespace std;

std::ofstream out_recola;
std::ofstream out_sub, out_real, out_finite;

double SqLam(double x,double y,double z)
{
  return abs(x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z);
}

NLO_Base::NLO_Base() {
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  p_nlodipoles = new YFS::Define_Dipoles();
  m_evts = 0;
  m_recola_evts = 0;
  m_realtool = 0;
  m_looptool = 0;
  if(m_isr_debug || m_fsr_debug){
  	m_histograms2d["IFI_EIKONAL"] = new Histogram_2D(0, -1., 1., 20, 0, 5., 20 );
  	m_histograms1d["Real_diff"] = new Histogram(0, -1, 1, 100);
  	m_histograms1d["Real_Flux"] = new Histogram(0, 0, 1, 100);
  	if (!ATOOLS::DirectoryExists(m_debugDIR_NLO)) {
			ATOOLS::MakeDir(m_debugDIR_NLO);
		}
  }
}


NLO_Base::~NLO_Base() {
  if(m_isr_debug || m_fsr_debug || m_check_real_sub){
		Histogram_2D * histo2d;
		string name;
		for (map<string, Histogram_2D *>::iterator hit = m_histograms2d.begin();
		        hit != m_histograms2d.end(); hit++) {
			histo2d = hit->second;
			name  = string(m_debugDIR_NLO) + "/" + hit->first + string(".dat");
			// histo2d->MPISync();
			histo2d->Finalize();
			histo2d->Output(name);
			delete histo2d;
		}
		Histogram * histo1d;
		for (map<string, Histogram *>::iterator hit = m_histograms1d.begin();
		        hit != m_histograms1d.end(); hit++) {
			histo1d = hit->second;
			name  = string(m_debugDIR_NLO) +  "/" + hit->first + string(".dat");
			histo1d->MPISync();
			histo1d->Finalize();
			histo1d->Output(name);
			delete histo1d;
		}
	}
	if(p_yfsFormFact) delete p_yfsFormFact;
	if(p_nlodipoles) delete p_nlodipoles;
	// if(p_real) delete p_real;
	// if(p_virt) delete p_virt;
}


void NLO_Base::InitializeVirtual(const PHASIC::Process_Info& pi) {
	p_virt = new YFS::Virtual(pi);
	m_looptool = true;
}

void NLO_Base::InitializeReal(const PHASIC::Process_Info& pi) {
	p_real = new YFS::Real(pi);
	m_realtool = true;
}

void NLO_Base::Init(Flavour_Vector &flavs, Vec4D_Vector &plab, Vec4D_Vector &born) {
	m_flavs = flavs;
	m_plab = plab;
	m_bornMomenta = born;
}


double NLO_Base::CalculateNLO() {
	double result{0.0};
	result += CalculateVirtual();
	result += CalculateReal();
	return result;
}


double NLO_Base::CalculateVirtual() {
	if (!m_looptool && !m_realvirt) return 0;
	double virt;
	double sub;
	Vec4D_Vector p = m_plab;
	CheckMassReg();
	// for(auto pp: m_plab) PRINT_VAR(pp.Mass());
	if(!HasISR()) virt = p_virt->Calc(m_bornMomenta, m_born);
	else virt = p_virt->Calc(p, m_born);
	if(m_check_virt_born) {
			if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-6)) {
			msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
									<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
									<< "\nSherpa" << ":  "<<m_born<<std::endl
									<<"PhaseSpace Point = ";
			for(auto _p: m_plab) msg_Error()<<_p<<std::endl;
		}
	}	
	sub = p_dipoles->CalculateVirtualSub();
	m_oneloop = (virt - sub * m_born);
	if(IsBad(m_oneloop)){
		msg_Error()<<"YFS Virtual is NaN"<<std::endl
							 <<"Virtual:  "<<m_oneloop<<std::endl
							 <<"Subtraction: "<<sub*m_born<<std::endl;
	}
	return m_oneloop;
}


double NLO_Base::CalculateReal() {
	if (!m_realtool) return 0;
	double real(0);
	if(m_coll_real) return p_dipoles->CalculateEEX()*m_born;
	for (auto k : m_ISRPhotons) {
		if(m_check_real_sub) {
			if(k.E() < 0.2*sqrt(m_s)) continue;
				CheckRealSub(k);
		}
		real+=CalculateReal(k);
	}
	for (auto k : m_FSRPhotons) {
		if(m_check_real_sub) {
			if(k.E() < 0.2*sqrt(m_s)) continue;
				CheckRealSub(k);
		}
		real+=CalculateReal(k,1);
	}
	if(IsBad(real)){
		msg_Error()<<"YFS Real is NaN"<<std::endl;
	}
	return real;
}


double NLO_Base::CalculateReal(Vec4D k, int submode) {
	double norm = 2.*pow(2 * M_PI, 3);
	Vec4D_Vector p(m_plab),pi(m_bornMomenta), pf(m_bornMomenta);
	Vec4D kk = k;
	MapMomenta(p, k);
	m_evts+=1;
	p_nlodipoles->MakeDipoles(m_flavs,p,m_plab);
	p_nlodipoles->MakeDipolesII(m_flavs,p,m_plab);
	p_nlodipoles->MakeDipolesIF(m_flavs,p,m_plab);

	double flux;
	if(m_flux_mode==1) flux = p_nlodipoles->CalculateFlux(k);
	else if(m_flux_mode==2) flux = 0.5*(p_dipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
	else flux = p_dipoles->CalculateFlux(kk);
	double tot,rcoll;
	double subloc = p_nlodipoles->CalculateRealSub(k);
	double subb   = p_dipoles->CalculateRealSubEEX(kk);
	if(IsZero(subb)) return 0;
	if(m_isr_debug || m_fsr_debug) m_histograms2d["IFI_EIKONAL"]->Insert(k.Y(),k.PPerp(), p_nlodipoles->CalculateRealSubIF(k));
	p.push_back(k);
	// if(submode!=1) flux = 1;
	// CheckMasses(p,1);
	CheckMomentumConservation(p);
	double r = p_real->Calc_R(p) / norm * flux;
	if(IsZero(r)) return 0;
	if(IsBad(r) || IsBad(flux)) {
		msg_Error()<<"Bad point for YFS Real"<<std::endl
							 <<"Real ME is : "<<r<<std::endl
							 <<"Flux is : "<<flux<<std::endl;
		return 0;
	}
	m_recola_evts+=1;
	// if(submode) tot = r-subloc*m_born;
	// else tot =  (r-subloc*m_born)/subloc;
	if(m_submode==submode::local) tot =  (r-subloc*m_born)/subloc;
	else if(m_submode==submode::global) tot =  (r-subloc*m_born)/subb;
	else if(m_submode==submode::off) tot =  (r)/subb;
	else msg_Error()<<METHOD<<" Unknown YFS Subtraction Mode "<<m_submode<<std::endl;
  if(m_isr_debug || m_fsr_debug){
		double diff = ((r/subloc - m_born)-( rcoll/subb - m_born))/((r/subloc - m_born)+( rcoll/subb - m_born));
		m_histograms1d["Real_diff"]->Insert(diff);
		m_histograms1d["Real_Flux"]->Insert(flux);
  }
  if(m_no_subtraction) return r/subloc;
  if(IsBad(tot)){
  	msg_Error()<<"NLO real is NaN"<<std::endl
  							<<"R = "<<r<<std::endl
  							<<"Local  S = "<<subloc*m_born<<std::endl
  							<<"GLobal S = "<<subb<<std::endl;
  }
	return tot;// / flux;
}

void NLO_Base::RandomRotate(Vec4D &p){
  Vec4D t1 = p;
  // rotate around x
  p[2] = cos(m_ranTheta)*t1[2] - sin(m_ranTheta)*t1[3];
  p[3] = sin(m_ranTheta)*t1[2] + cos(m_ranTheta)*t1[3];
  t1 = p;
  // rotate around z
  p[1] = cos(m_ranPhi)*t1[1]-sin(m_ranPhi)*t1[2];
  p[2] = sin(m_ranPhi)*t1[1]+cos(m_ranPhi)*t1[2];
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k) {
	Vec4D Q;
	Vec4D QQ, PP;
	Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  double s = (m_plab[0]+m_plab[1]).Abs2();
  double t = (m_plab[0]-m_plab[2]).Abs2();
  m_ranTheta = acos(1.+2.*t/s);
	m_ranPhi = ran->Get()*2.*M_PI;
	// Poincare boostLab(p[0] + p[1]);
	for (int i = 2; i < p.size(); ++i)
	{
		Q += p[i];
	}
	Q += k;
	double sq = Q.Abs2();
	Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
	for (int i = 2; i < p.size(); ++i) {
		boostQ.Boost(p[i]);
		// pRot.Rotate(p[i]);
		// RandomRotate(p[i]);
	}
	boostQ.Boost(k);
	// pRot.Rotate(k);
	// RandomRotate(k);
	double qx(0), qy(0), qz(0);
	for (int i = 2; i < p.size(); ++i)
	{
		qx += p[i][1];
		qy += p[i][2];
		qz += p[i][3];
	}
	if (!IsEqual(k[1], -qx, 1e-5) || !IsEqual(k[2], -qy, 1e-5) || !IsEqual(k[3], -qz, 1e-5) ) {
		if( k[1]> 1e-6 && k[2]> 1e-6 && k[3]> 1e-6 ){
			msg_Error() << "YFS Mapping has failed for ISR\n";
			msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx << std::endl;
			msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy << std::endl;
			msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz << std::endl;
		}
		}
	for (int i = 2; i < p.size(); ++i)
	{
		QQ += p[i];
	}
	QQ+=k;
	double sqq = QQ.Abs2();
	if (!IsEqual(sqq, sq, 1e-8))
	{
		msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
	}
	// if(m_is_isr) QQ = p[0]+p[1];
  // double zz = sqrt(sqq) / 2.;
	// double z = zz * sqrt((sqq - sqr(m_flavs[0].Mass() - m_flavs[1].Mass())) * (sqq - sqr(m_flavs[0].Mass() + m_flavs[1].Mass()))) / sqq;
	double sign_z = (p[0][3] < 0 ? -1 : 1);
	// p[0] = {zz, 0, 0, z};
	// p[1] = {zz, 0, 0, -z};
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamCM = 0.5*sqrt(SqLam(sqq,m1*m1,m2*m2)/sqq);
  double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
  double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
 	p[0] = {E1, 0, 0, sign_z*lamCM};
  p[1] = {E2, 0, 0, -sign_z*lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 	0., 1., 0.));
	for (int i = 0; i < p.size(); ++i)
	{
		// pRot2.Rotate(p[i]);
		boostLab.Boost(p[i]);
	}
	// pRot2.Rotate(k);
	boostLab.Boost(k);
}


void NLO_Base::CheckMasses(Vec4D_Vector &p, int realmode){
	bool allonshell=true;
	std::vector<double> masses;
	Flavour_Vector flavs = m_flavs;
	if(realmode) flavs.push_back(Flavour(kf_photon));
	for (int i = 0; i < p.size(); ++i)
	{
		masses.push_back(flavs[i].Mass());
		if(!IsEqual(p[i].Mass(),flavs[i].Mass(),1e-6)){
			msg_Debugging()<<"Wrong particle masses in YFS Mapping"<<std::endl
								 <<"Flavour = "<<flavs[i]<<", with mass = "<<flavs[i].Mass()<<std::endl
								 <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
			allonshell = false;

		}
	}
	if(!allonshell) m_stretcher.StretchMomenta(p, masses);
	// return true;
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
	for (int i = 0; i < m_plab.size(); ++i)
	{
		if (m_flavs[i].IsChargedLepton()) {
			double sik = (k + m_plab[i]).Abs2();
			if (sik < m_hardmin ) {
				return false;
			}
		}
	}
	return true;
}


bool NLO_Base::CheckMomentumConservation(Vec4D_Vector p){
  Vec4D incoming = p[0]+p[1];
  Vec4D outgoing;
  for (int i = 2; i < p.size(); ++i)
  {
    outgoing+=p[i];
  }
  Vec4D diff = incoming - outgoing;
  if(!IsEqual(incoming,outgoing, 1e-5)){
    msg_Error()<<METHOD<<std::endl<<"Momentum not conserverd in YFS"<<std::endl
               <<"Incoming momentum = "<<incoming<<std::endl
               <<"Outgoing momentum = "<<outgoing<<std::endl
               <<"Difference = "<<diff<<std::endl
               <<"Vetoing Event "<<std::endl;
  }
  return true;
}

void NLO_Base::CheckMassReg(){
	double virt;
	if (m_check_mass_reg==1 && !m_realvirt) {
		out_sub.open("yfs-sub.txt", std::ios_base::app);
		out_recola.open("recola-res.txt", std::ios_base::app); // append instead of overwrite
		out_finite.open("yfs-finite.txt", std::ios_base::app);
		if(!HasISR()) virt = p_virt->Calc(m_bornMomenta, m_born);
		else virt = p_virt->Calc(m_plab, m_born);
		if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-6)) {
			msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
									<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
									<< "Sherpa" << ":  "<<m_born;
		}
		double sub = p_dipoles->CalculateVirtualSub();
		std::cout << setprecision(15);
		out_sub<< setprecision(15) << m_photonMass << "," << -sub*m_born << std::endl;
		out_recola<< setprecision(15) << m_photonMass << "," << virt << std::endl;
		out_finite<< setprecision(15) << m_photonMass << "," << virt - sub*m_born << std::endl;
		out_sub.close();
		out_recola.close();
		exit(0);
	}
}


void NLO_Base::CheckRealSub(Vec4D k){
		// if(k.E() < 20) return;
		// k*=100;
		double real;
		std::string filename="Real_subtracted_";
		for(auto f: m_flavs) {
			filename+=f.IDName();
			filename+="_";
		}
		filename+=".txt";
		if(ATOOLS::FileExists(filename))  ATOOLS::Remove(filename);
		out_sub.open(filename, std::ios_base::app);
		// if(k.E() < 0.8*sqrt(m_s)/2.) return;
		for (double i = 1; i < 20 ; i+=0.25)
		{
			k=k/i;
			real=CalculateReal(k);
			out_sub<<k.E()<<","<<fabs(real)<<std::endl;
			if(k.E() < 1e-10 || real==0) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		exit(0);
}

