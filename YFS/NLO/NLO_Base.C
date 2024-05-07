YFS/NLO/NLO_Base.C// #include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/NLO_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"


using namespace YFS;
using namespace MODEL;

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
  m_realrealtool = 0;
  m_realtool = 0;
  m_looptool = 0;
  m_realvirt = 0;
  // if(m_check_real_sub) {
  // 	m_histograms2d["Real_me_sub"] = new Histogram_2D(110, 1e-6, sqrt(m_s)/2., 50, 1e-6, 100, 50);
  // 	if (!ATOOLS::DirectoryExists(m_debugDIR_NLO)) {
	// 		ATOOLS::MakeDir(m_debugDIR_NLO);
  // 	}
 	// }
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
}


void NLO_Base::InitializeVirtual(const PHASIC::Process_Info& pi) {
	if (m_griff != 0) {
		p_griffin = new Griffin::Griffin_Interface();
		p_griffin->Initialize(pi);
		return;
	}
	p_virt = new YFS::Virtual(pi);
	m_looptool = true;
}

void NLO_Base::InitializeReal(const PHASIC::Process_Info& pi) {
	p_real = new YFS::Real(pi);
	m_realtool = true;
}

void NLO_Base::InitializeRealVirtual(const PHASIC::Process_Info& pi) {
	p_realvirt = new YFS::RealVirtual(pi);
	m_realvirt = true;
	// m_looptool = true;
}

void NLO_Base::InitializeRealReal(const PHASIC::Process_Info& pi) {
	p_realreal = new YFS::RealReal(pi);
	m_realrealtool = true;
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
	result += CalculateRealVirtual();
	result += CalculateRealReal();
	return result;
}


double NLO_Base::CalculateVirtual() {
	if (!m_looptool && !m_realvirt) return 0;
	if (m_griff != 0) {
		return p_griffin->EvaluateLoop(m_plab)*m_born;
	}
	double virt;
	double sub;
	Vec4D_Vector p = m_plab;
	CheckMassReg();
	// for(auto pp: m_plab) PRINT_VAR(pp.Mass());
	if(m_fsrmode==2) virt = p_virt->Calc(m_bornMomenta, m_born);
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
		msg_Error()<<"YFS Virtual is NaN"<<std::endl;
	}
	return m_oneloop;
}


double NLO_Base::CalculateReal() {
	if (!m_realtool) return 0;
	double real(0);
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

double NLO_Base::CalculateRealVirtual() {
	if (!m_realvirt) return 0;
	double real(0), sub(0);
	double norm = 2.*pow(2 * M_PI, 3);
	double flux;
	double mz = Flavour(kf_Z).Mass();
	double gz = Flavour(kf_Z).Width();
	Vec4D Q = m_bornMomenta[0]+m_bornMomenta[1];
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto k : m_FSRPhotons) photons.push_back(k);
	for (auto k : photons) {
		Vec4D_Vector p(m_plab);
		MapMomenta(p, k);
		if(CheckPhotonCollinear(p,k)){
			// if(k.PPerp()<1e-1) return 0;
			flux  = p_nlodipoles->CalculateFlux(k);
			p_nlodipoles->MakeDipoles(m_flavs,p,m_plab);
			p_nlodipoles->MakeDipolesII(m_flavs,p,m_bornMomenta);
			p_nlodipoles->MakeDipolesIF(m_flavs,p,m_plab);
			double subb   = p_dipoles->CalculateRealSub(k);
			double subloc = p_nlodipoles->CalculateRealSub(k);
			p.push_back(k);
			CheckMasses(p,1);
			CheckMassRegRV(p);
			double r = p_realvirt->Calc(p, m_born) / norm*flux;
			if (r == 0 || IsBad(r)) continue;
			double aB = subloc*CalculateVirtual();
			double tot = (r-aB) / subloc;
			real += tot;
		}
	}
	return real;
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

double NLO_Base::CalculateRealReal() {
	if (!m_realrealtool) return 0;
	double real(0), sub(0), subflux1, subflux2;
	Vec4D_Vector photons, p(m_plab);
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto k : m_FSRPhotons) photons.push_back(k);
	// if (NHardPhotons(photons) == 0) return 0;
	double norm = 2.*pow(2 * M_PI, 3);
	double mz = Flavour(kf_Z).Mass();
	double gz = Flavour(kf_Z).Width();
	Vec4D Q = m_bornMomenta[0]+m_bornMomenta[1];
	// Vec4D_Vector photons;
	double flux;
	for (int i = 0; i < photons.size(); ++i) {
		for (int j = 0; j < i; ++j) {
			p = m_plab;
			Vec4D k  = photons[i];
			Vec4D kk = photons[j];
			// if(m_check_rr){
			// 	CheckReaReallSub(k,kk);
			// }
			MapMomenta(p, k, kk);
			CheckMasses(p);
			// flux  = (Q-kk-k).Abs2()/(Q).Abs2();
			double flux1,flux2;
			if(m_flux_mode==1) {
				flux1 = p_dipoles->CalculateFlux(k);
				flux2 = p_dipoles->CalculateFlux(kk);
			}
			// else if(m_flux_mode==2) {
			// 	flux1 = 0.5*(p_dipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
			// 	flux 2= 0.5*(p_dipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
			// }
			else {
				flux1 = p_nlodipoles->CalculateFlux(k);
				flux2 = p_nlodipoles->CalculateFlux(kk);
			}
			// if(k.PPerp()<0.1 || kk.PPerp()<0.1) return 0;
			if(!CheckPhotonForReal(k)||!CheckPhotonForReal(kk)) return 0;

			Vec4D ksum = k + kk;
			p_nlodipoles->MakeDipoles(m_flavs,p,m_plab);
			p_nlodipoles->MakeDipolesII(m_flavs,p,m_bornMomenta);
			p_nlodipoles->MakeDipolesIF(m_flavs,p,m_plab);
			double subloc1 = p_nlodipoles->CalculateRealSub(k);
			double subloc2 = p_nlodipoles->CalculateRealSub(kk);
			double subb1   = p_dipoles->CalculateRealSub(k);
			double subb2   = p_dipoles->CalculateRealSub(kk);	
			p.push_back(k);
			p.push_back(kk);
			CheckMomentumConservation(p);
			if(m_fsrmode==1 || m_noflux!=0){
				Q*=0;
				for (int i = 2; i < m_plab.size(); ++i)
				{
					Q += m_plab[i];
				}
				double sq = (Q).Abs2();
				double sx = (Q+k).Abs2();
				subflux1 = (sq*sq)/(sx*sx)-1;
				sx = (Q+kk).Abs2();
				subflux2 = (sq*sq)/(sx*sx)-1;
				// if(m_ifisub==1 && shifdiff > m_pole_fac*gz) subloc+=p_nlodipoles->CalculateRealSubIF(k);
				if(m_ifisub==1) {
					subloc1+=p_nlodipoles->CalculateRealSubIF(k);
					subloc2+=p_nlodipoles->CalculateRealSubIF(kk);
				}
				if(m_ifisub==2) {
					subloc1+=p_nlodipoles->CalculateRealSubIF(k)*subflux1;
					subloc2+=p_nlodipoles->CalculateRealSubIF(kk)*subflux2;
				}
				if(m_isr_debug || m_fsr_debug) m_histograms2d["IFI_EIKONAL"]->Insert(k.Y(),k.PPerp(), p_nlodipoles->CalculateRealSubIF(k));
			}
			double subb = subloc1*CalculateReal(kk,1);
			subb += subloc2*CalculateReal(k,1);
			subb += subloc1*subloc2*m_born;
			double r = p_realreal->Calc_R(p) / norm / norm *flux1*flux2;
			if(IsNan(r)) r=0;
			if (r == 0) continue;
			real += (r - subb)/subb1/subb2;
		}
	}
	return real;
}


void NLO_Base::MapInitial(Vec4D_Vector &p, Vec4D &k){
	Vec4D Q;
  Q = p[0] + p[1] - k;
  // if(Q.Abs2() > )
  double sp = Q * Q;
  double zz = sqrt(sp) / 2.;
  double z = zz * sqrt((sp - sqr(p[0].Mass() - p[1].Mass())) * (sp - sqr(p[0].Mass() + p[1].Mass()))) / sp;
  p[0] = {zz, 0, 0, z};
  p[1] = {zz, 0, 0, -z};
  ATOOLS::Poincare poin(Q);
  Poincare pRot(p[0], Vec4D(0., 0., 0., 1.));
  for (int i = 0; i < 2; ++i) {
    pRot.Rotate(p[i]);
    poin.BoostBack(p[i]);
  }
	// PHASIC::CE.Isotropic2Momenta(p[0]+p[1]-k, sqr(p[2].Mass()), sqr(p[3].Mass()), p[2], p[3], ran->Get(), ran->Get(), -1, 1);
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


void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k1, Vec4D &k2) {
	Vec4D Q;
	Vec4D QQ, PP;
	Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
	// Poincare boostLab(p[0] + p[1]);
	for (int i = 2; i < p.size(); ++i)
	{
		Q += p[i];
	}
	Q += k1+k2;
	double sq = Q.Abs2();
	Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
	for (int i = 2; i < p.size(); ++i) {
		pRot.Rotate(p[i]);
		boostQ.Boost(p[i]);
	}
	pRot.Rotate(k1);
	boostQ.Boost(k1);
	pRot.Rotate(k2);
	boostQ.Boost(k2);
	double qx(0), qy(0), qz(0);
	for (int i = 2; i < p.size(); ++i)
	{
		qx += p[i][1];
		qy += p[i][2];
		qz += p[i][3];
	}
	Vec4D k = k1+k2;
	if (!IsEqual(k[1], -qx, 1e-6) || !IsEqual(k[2], -qy, 1e-6) || !IsEqual(k[3], -qz, 1e-6) ) {
		if( k[1]> 1e-8 && k[2]> 1e-8 && k[3]> 1e-8 ){
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
	QQ+=k1+k2;
	double sqq = QQ.Abs2();
	if (!IsEqual(sqq, sq, 1e-8))
	{
		msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
	}
	// if(m_is_isr) QQ = p[0]+p[1];
 	double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamCM = 0.5*sqrt(SqLam(sqq,m1*m1,m2*m2)/sqq);
  double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
  double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
 	p[0] = {E1, 0, 0, lamCM};
  p[1] = {E2, 0, 0, -lamCM};
  Poincare pRot2(p[0], Vec4D(0., 	0., 0., 1.));
	for (int i = 0; i < p.size(); ++i)
	{
		boostLab.BoostBack(p[i]);
	}
	boostLab.BoostBack(k1);
	boostLab.BoostBack(k2);
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

bool NLO_Base::CheckPhotonCollinear(const Vec4D_Vector &p, const Vec4D &k){
	for(auto pp: p){
		double m2 = sqr(pp.Mass());
		double sik = (pp+k).Abs2();
		if(sik < 2*m2) return false;
	}
	return true;
}

void NLO_Base::MakeHardMomenta(){
	 // m_plab[0] = 
	m_eikmom = m_plab;
	 for (Dipole_Vector::iterator Dip = p_dipoles->GetDipoleFF()->begin();
         Dip != p_dipoles->GetDipoleFF()->end(); ++Dip) {
	 			for (int i = 0; i < Dip->m_eikmomentum.size(); ++i)
	 			{
	 				m_eikmom[i+2] = Dip->GetEikMomenta(i);
	 			}
	 }
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
		if(m_fsrmode==2) virt = p_virt->Calc(m_bornMomenta, m_born);
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

void NLO_Base::CheckMassRegRV(const ATOOLS::Vec4D_Vector &p){
		double virt;
		if (m_check_mass_reg && m_realvirt) {
			out_sub.open("yfs-RV-sub.txt", std::ios_base::app);
			out_recola.open("recola-RV-res.txt", std::ios_base::app); // append instead of overwrite
			out_finite.open("yfs-RV-finite.txt", std::ios_base::app);
			virt = p_realvirt->Calc(p, m_born);
			double sub = p_dipoles->CalculateRealSub(p.back())*CalculateVirtual();
			std::cout << setprecision(15);
			out_sub<< setprecision(15) << m_photonMass << "," << -sub*m_born << std::endl;
			out_recola<< setprecision(15) << m_photonMass << "," << virt << std::endl;
			out_finite<< setprecision(15) << m_photonMass << "," << virt - sub*m_born << std::endl;
			out_sub.close();
			out_recola.close();
			exit(1);
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

void NLO_Base::CheckRealRealSub(Vec4D k1, Vec4D k2){
		// if(k.E() < 20) return;
		// // k*=100;
		// Vec4D _k1 = k1;
		// Vec4D _k2 = k2;
		// double real;
		// std::string filename1="Real_Real_subtracted_k1_fx";
		// std::string filename2="Real_Real_subtracted_k2_fx";
		// std::string filename3="Real_Real_subtracted_";
		// for(auto f: m_flavs) {
		// 	filename1+=f.IDName();
		// 	filename2+=f.IDName();
		// 	filename3+=f.IDName();
		// 	filename1+="_";
		// 	filename2+="_";
		// 	filename3+="_";
		// }
		// filename1+=".txt";
		// filename2+=".txt";
		// filename3+=".txt";
		// if(ATOOLS::FileExists(filename1))  ATOOLS::Remove(filename1);
		// if(ATOOLS::FileExists(filename2))  ATOOLS::Remove(filename2);
		// if(ATOOLS::FileExists(filename3))  ATOOLS::Remove(filename3);

		// out_sub.open(filename1, std::ios_base::app);
		// if(k1.E() < 0.8*sqrt(m_s)/2.) return;
		// for (double i = 1; i < 20 ; i+=0.25)
		// {
		// 	k1=k1/i;
		// 	real=CalculateRealReal(k1,k2);
		// 	out_sub<<k1.E()<<","<<fabs(real)<<std::endl;
		// 	if(k1.E() < 1e-10 || real==0) break;
		// 	// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		// }
		// out_sub.close();
		// k1 = _k1;
		// out_sub.open(filename2, std::ios_base::app);
		// if(k2.E() < 0.8*sqrt(m_s)/2.) return;
		// for (double i = 1; i < 20 ; i+=0.25)
		// {
		// 	k2=k2/i;
		// 	real=CalculateRealReal(k1,k2);
		// 	out_sub<<k2.E()<<","<<fabs(real)<<std::endl;
		// 	if(k2.E() < 1e-10 || real==0) break;
		// 	// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		// }
		// out_sub.close();
		// k1 = _k1;
		// k2 = _k2;
		// out_sub.open(filename3, std::ios_base::app);
		// if(k2.E() < 0.8*sqrt(m_s)/2.) return;
		// for (double i = 1; i < 20 ; i+=0.25)
		// {
		// 	k2=k2/i;
		// 	k1=k1/i;
		// 	real=CalculateRealReal(k1,k2);
		// 	out_sub<<k2.E()<<","<<fabs(real)<<std::endl;
		// 	if(k2.E() < 1e-10 || real==0) break;
		// 	if(k1.E() < 1e-10 || real==0) break;
		// 	// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		// }
		exit(0);
}

