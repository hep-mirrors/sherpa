#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/NLO_Base.H"


using namespace YFS;

std::ofstream out_recola;
std::ofstream out_sub, out_real, out_finite;

NLO_Base::NLO_Base() {
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  p_nlodipoles = new YFS::Define_Dipoles();
  m_evts = 0;
  m_recola_evts = 0;
  m_realrealtool = 0;
  m_realtool = 0;
  m_looptool = 0;
  m_realvirt = 0;
  if(m_isr_debug || m_fsr_debug){
  	m_histograms2d["Real_me"] = new Histogram_2D(0, -1, 1, 20, 0, sqrt(m_s)/2., 20 );
  	m_histograms1d["Real_diff"] = new Histogram(0, -1, 1, 40);
  }
}


NLO_Base::~NLO_Base() {
	// if (p_real) delete p_real;
	// if (p_virt) delete p_virt;
	// if (p_realvirt) delete p_realvirt;
	// if (p_realreal) delete p_realreal;
	// if (p_griffin)   delete p_griffin;
  if(m_isr_debug || m_fsr_debug){
		Histogram_2D * histo2d;
		string name;
		for (map<string, Histogram_2D *>::iterator hit = m_histograms2d.begin();
		        hit != m_histograms2d.end(); hit++) {
			histo2d = hit->second;
			name  = string("./YFS_NLO_Hist/") + hit->first + string(".dat");
			// histo2d->MPISync();
			histo2d->Finalize();
			histo2d->Output(name);
			delete histo2d;
		}
		Histogram * histo1d;
		for (map<string, Histogram *>::iterator hit = m_histograms1d.begin();
		        hit != m_histograms1d.end(); hit++) {
			histo1d = hit->second;
			name  = string("./YFS_NLO_Hist/") + hit->first + string(".dat");
			histo1d->MPISync();
			histo1d->Finalize();
			histo1d->Output(name);
			delete histo1d;
		}
	}
	PRINT_VAR(m_recola_evts);
	PRINT_VAR(m_evts);
	msg_Out()<<"Percentage of Recola events = "<<m_recola_evts/m_evts*100.<<"% "<<std::endl;
}


void NLO_Base::InitializeVirtual(const PHASIC::Process_Info& pi) {
	if (m_griff != 0) {
		p_griffin = new Griffin::Griffin_Interface();
		p_griffin->Initialize(pi);
		return;
	}
	p_virt = new YFS::Virtual(pi);
	m_looptool = true;
	m_real_only = true;
}

void NLO_Base::InitializeReal(const PHASIC::Process_Info& pi) {
	p_real = new YFS::Real(pi);
	m_realtool = true;
}

void NLO_Base::InitializeRealVirtual(const PHASIC::Process_Info& pi) {
	p_realvirt = new YFS::RealVirtual(pi);
	m_realvirt = true;
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
	if (m_griff != 0) {
		double v =   p_griffin->EvaluateLoop(m_plab);
		// PRINT_VAR(v);
		return v;
	}
	if (!m_looptool) return 0;
	double virt;
	double sub;
	if (m_check_mass_reg) {
		out_sub.open("yfs-sub.txt", std::ios_base::app);
		out_recola.open("recola-res.txt", std::ios_base::app); // append instead of overwrite
		out_finite.open("yfs-finite.txt", std::ios_base::app);
		virt = p_virt->Calc(m_plab, m_born);
		if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-4)) {
			msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
									<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
									<< "Sherpa" << ":  "<<m_born;
		}
		// born = p_virt->p_loop_me->ME_Born();
		sub = p_dipoles->CalculateVirtualSub();
		std::cout << setprecision(16);
		out_sub << m_photonMass << "," << -sub*m_born << std::endl;
		out_recola << m_photonMass << "," << virt << std::endl;
		out_finite << m_photonMass << "," << virt  -sub*m_born << std::endl;
		out_sub.close();
		out_recola.close();
		exit(1);
	}
	// if (m_check_poles) {
	// 	out_sub.open("yfs-sub.txt", std::ios_base::app);
	// 	virt = p_virt->Calc(m_plab, m_born);
	// 	DivArrC  subd = p_yfsFormFact->BVV_full_eps(m_plab[0], m_plab[1], m_photonMass, sqrt(m_s) / 2., 3);
	// 	PRINT_VAR(subd.GetEpsilon());
	// 	out_sub << rpa->gen.Ecms() << "," << -sub*m_born << std::endl;
	// 	out_sub.close();
	// 	exit(1);
	// }
	virt = p_virt->Calc(m_plab, m_born);
	if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-4)) {
		msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
								<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
								<< "\nSherpa" << ":  "<<m_born<<std::endl
								<<"PhaseSpace Point = ";
	for(auto p: m_plab) msg_Error()<<p<<std::endl;
	}
	sub = p_dipoles->CalculateVirtualSub();
	m_oneloop = (virt - sub * m_born);
	return m_oneloop;
}


double NLO_Base::CalculateReal() {
	if (!m_realtool) return 0;
	double real(0), sub(0);
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	Poincare Boost(m_bornMomenta[0]+m_bornMomenta[1]);
	for (auto kk : m_FSRPhotons) { 
		// Boost.BoostBack(kk);
		// photons.push_back(kk);
		// if(kk.E() > 4.56e-05) photons.push_back(kk);
	}
	for (auto k : photons) {
		// Boost.Boost(k);
		double tot = CalculateReal(k);
		real += tot;
	}
	return real;
}


double NLO_Base::CalculateReal(Vec4D &k) {
	double norm = 2.*pow(2 * M_PI, 3);
	m_sp = (m_plab[0]+m_plab[1]).Abs2();
	Vec4D_Vector p(m_plab), pp(m_bornMomenta);
	Vec4D kk = k;
	double flux = (p[2]+p[3]+k).Abs2()/(p[0]+p[1]).Abs2();
	MapMomenta(p, k);
		// PRINT_VAR(kk);
	pp=p;
	// p_virt->Calc(p, m_born);
	// double born = p_virt->p_loop_me->ME_Born();
	double born = m_born;
	double tot,colltot,rcoll;
	p_nlodipoles->MakeDipoles(m_flavs,p,m_bornMomenta);
	p_nlodipoles->MakeDipolesII(m_flavs,p,m_bornMomenta);
	p_nlodipoles->MakeDipolesIF(m_flavs,p,m_bornMomenta);
	// p_dipoles->MakeDipolesIF(m_flavs, m_plab, p);
	// p_dipoles->MakeDipoles(m_flavs, m_plab, p);
	double subb   = p_dipoles->CalculateRealSub(k);
	double subloc = p_nlodipoles->CalculateRealSub(k);
	m_evts+=1;
	if (!CheckPhotonForReal(k)) { 
		rcoll = p_dipoles->CalculateEEXReal(k)*born/M_PI/2.;
		return 0*(rcoll -subloc*born)/subb;
	}
	p.push_back(k);
	double r = p_real->Calc_R(p) / norm;
	if (r == 0) return 0;
	m_recola_evts+=1;
	tot =  ( r - subloc*born)/subb;
  if(m_isr_debug || m_fsr_debug){

		m_histograms2d["Real_me"]->Insert(k.CosTheta(m_bornMomenta[0]), k.E(), r/rcoll);
		double diff = (tot-colltot)/(tot+colltot);
		m_histograms1d["Real_diff"]->Insert(diff);
  }
	return tot;
}

double NLO_Base::CollinearReal(Vec4D k, Vec4D_Vector p){
	// double born = p_virt->p_loop_me->ME_Born();
	Poincare poin(p[0]+p[1]);
	// for(auto _p: p) poin.Boost(_p);
	// poin.Boost(k);
	double born = m_born;
	double p1p2 = p[0]*p[1];
	double a = k*p[1]/p1p2;
	double b = k*p[0]/p1p2;
	double coeff = m_alpha/(4*M_PI*M_PI);
	coeff *= 2.*p[0]*p[1];
	coeff /= k*p[0];
	coeff /= k*p[1];
	double den = sqr(1-a)+sqr(1-b);
	double Af = p[0].Mass()*p[1].Mass()/(2.*p[0]*p[1]);
	coeff *= (1.-Af*((1.-a)*(1.-b))/den*(a/b+b/a));
	double V = 1;//+0.5*p_dipoles->GetDipoleII()->m_gamma;
	// if(m_looptool) V+=0.5*p_dipoles->GetDipoleII()->m_gamma;
	double z = m_sp/m_s;
	double L = log(m_s/p[0].Mass()*p[1].Mass());
	// return m_alpi*((1+z*z)/(1-z)*(L-1.))*born;
	// return coeff*Eikonal(k,p[0],p[1])*(sqr(1.-a)+sqr(1.-b))*born;
	return 0.5*coeff*(sqr(1.-a)+sqr(1.-b))*born;
}

double NLO_Base::CalculateRealVirtual() {
	if (!m_realvirt) return 0;
	double real(0), sub(0);
	double norm = 2.*pow(2 * M_PI, 3);
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	// for (auto k : m_FSRPhotons) photons.push_back(k);
	for (auto k : photons) {
		Vec4D_Vector p(m_plab);
		MapMomenta(p, k);
		// if(k.PPerp()<5) return 0;
		if(k.PPerp()<1) return 0;

		// if(!CheckPhotonForReal(k)) return 0;
		// if(k.E() < 0.01) continue;
		// double subb  = p_dipoles->CalculateRealSub(k);
		// double subloc = 0;
		p_nlodipoles->MakeDipoles(m_flavs,p,m_bornMomenta);
		p_nlodipoles->MakeDipolesII(m_flavs,p,m_bornMomenta);
		p_nlodipoles->MakeDipolesIF(m_flavs,p,m_bornMomenta);
		// p_dipoles->MakeDipolesIF(m_flavs, m_plab, p);
		// p_dipoles->MakeDipoles(m_flavs, m_plab, p);
		double subb   = p_dipoles->CalculateRealSub(k);
		double subloc = p_nlodipoles->CalculateRealSub(k);
		p.push_back(k);
		double r = p_realvirt->Calc(p, m_born);
		if (r == 0 || IsBad(r)) continue;
		double aB = p_dipoles->CalculateRealSub(k)*m_oneloop;
		// PRINT_VAR(aB);
		// PRINT_VAR(r);
		double tot = (r-aB) / p_dipoles->CalculateRealSub(k);
		// tot -= m_born * p_yfsFormFact->BVV_full(p[0], p[1], m_photonMass, sqrt(m_s) / 2., 3);
		real += tot;
	}
	return real/norm;
}


double NLO_Base::CalculateRealReal() {
	if (!m_realrealtool) return 0;
	// if (m_ISRPhotons.size() <= 1 && m_FSRPhotons.size() <= 1) return 0;
	double real(0), sub(0);
	Vec4D_Vector photons, p(m_plab);
	for (auto k : m_ISRPhotons) photons.push_back(k);
	// for (auto k : m_FSRPhotons) photons.push_back(k);
	// if (NHardPhotons(photons) == 0) return 0;
	double norm = 2.*pow(2 * M_PI, 3);
	for (int i = 0; i < photons.size(); ++i) {
		for (int j = 0; j < i; ++j) {
			Vec4D k  = photons[i];
			Vec4D kk = photons[j];
			if(k.PPerp()<0.1 || kk.PPerp()<0.1) return 0;
			// if (!CheckPhotonForReal(k) || !CheckPhotonForReal(kk)) {
			// 	return 0;
			// 	// continue;
			// }
			// PRINT_VAR("HERE");
			Vec4D ksum = k + kk;
			p = m_plab;
			MapMomenta(p, ksum);
			// MapMomenta(p, kk);
			p_nlodipoles->MakeDipoles(m_flavs,p,m_bornMomenta);
			p_nlodipoles->MakeDipolesII(m_flavs,p,p);
			p_nlodipoles->MakeDipolesIF(m_flavs,p,m_bornMomenta);
			double subloc1 = p_nlodipoles->CalculateRealSub(k);
			double subloc2 = p_nlodipoles->CalculateRealSub(kk);
			double subb1   = p_dipoles->CalculateRealSub(k);
			double subb2   = p_dipoles->CalculateRealSub(kk);

			double subb = subloc1*CalculateReal(kk)*subb2;
			subb += subloc2*CalculateReal(k)*subb1;
			subb += subloc1*subloc2*m_born;
			p.push_back(k);
			p.push_back(kk);
			double r = p_realreal->Calc_R(p) / norm / norm;
			if(IsNan(r)) r=0;
			if (r == 0) continue;
			real += (r - subb)/subb1/subb2;
		}
	}
	return real;
}




void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k) {
	Vec4D Q;
	Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
	for (int i = 2; i < p.size(); ++i)
	{
		Q += p[i];
	}
	double sq = Q.Abs2();
	Q += k;
	Poincare boostQ(Q);
	for (int i = 2; i < p.size(); ++i) {
		boostQ.Boost(p[i]);
	}
	// boostQ.Boost(k);
	double qx(0), qy(0), qz(0);
	for (int i = 2; i < p.size(); ++i)
	{
		qx += p[i][1];
		qy += p[i][2];
		qz += p[i][3];
	}
	// if (!IsEqual(k[1], -qx, 1e-3) || !IsEqual(k[2], -qy, 1e-3) || !IsEqual(k[3], -qz, 1e-3) ) {
	// 	if( k[1]> 1e-5 && k[2]> 1e-5 && k[3]> 1e-5 ){
	// 		msg_Error() << "YFS Mapping has failed for ISR\n";
	// 		msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx << std::endl;
	// 		msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy << std::endl;
	// 		msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz << std::endl;
	// 	}
	// }
	Vec4D QQ, PP;
	for (int i = 2; i < p.size(); ++i)
	{
		QQ += p[i];
	}
	double sqq = QQ.Abs2();
	if (!IsEqual(sqq, sq, 1e-4))
	{
		msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
	}
	QQ += k;
	PHASIC::CE.Isotropic2Momenta(QQ, sqr(p[0].Mass()), sqr(p[1].Mass()), p[0], p[1], ran->Get(), ran->Get(), -1, 1);
	for (int i = 0; i < p.size(); ++i)
	{
		boostLab.Boost(p[i]);
	}
	// boostLab.Boost(k);
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
	if (k.E() < m_hardmin) return false;
	if(k.PPerp() < m_hardmin) return false;
	for (int i = 0; i < m_bornMomenta.size(); ++i)
	{
		if (m_flavs[i].IsChargedLepton()) {
			double sik = (k + m_bornMomenta[i]).Abs2();
			if (sqrt(sik) < m_hardmin) {
				return false;
			}
		}
	}
	return true;
}