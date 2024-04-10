#include "YFS/Tools/Debug.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"


using namespace YFS;
using namespace std;

Debug::Debug() {
	InitializeHist();
}

Debug::~Debug() {
	WriteHistograms();
}


void Debug::InitializeHist() {
	if (m_isr_debug && m_fsrmode != 2) {
		if (!ATOOLS::DirectoryExists(m_debugDIR_ISR)) {
			ATOOLS::MakeDir(m_debugDIR_ISR);
		}
		m_histograms_ISR["lep-mass"]  = new Histogram(0, 0., sqrt(m_s), 50);
		m_histograms_ISR["v"]  = new Histogram(0, -0.9, 1.1, 20);
		m_histograms_ISR["K0"]  = new Histogram(0, -10., sqrt(m_s) + 10, 50);
		m_histograms_ISR["jacobian"]  = new Histogram(0, -1., 1, 50);
		m_histograms_ISR["A"]  = new Histogram(0, -1., 3, 50);
		m_histograms_ISR["K2"]  = new Histogram(0, -1., 20., 50);
		m_histograms_ISR["ptk"]  = new Histogram(0, -1., 5., 50);
		m_histograms_ISR["cos(theta)"]  = new Histogram(0, -1., 1., 20);
		m_histograms_ISR["lam"]  = new Histogram(0, 0, 2, 40);
		m_histograms_ISR["NGamma"]  = new Histogram(0, -1, 20, 21);
		m_histograms_ISR["Cutflow"]  = new Histogram(0, -1, 4, 5);
		m_histograms_ISR["nbar"]  = new Histogram(0, -2, 2, 40);
		m_histograms_ISR["weight"]  = new Histogram(0, -50, 50, 200);
		m_histograms_ISR["massweight"]  = new Histogram(0, 0, 2, 100);
		m_histograms_ISR["jacweight"]  = new Histogram(0, 0, 2, 100);
  		m_histograms2d["Form_Factor_FS_Angle"] = new Histogram_2D(0, 0, 30., 60 , 0.99, 1.01, 60 );
	}
	if (m_fsr_debug && m_fsrmode != 0) {

		if (!ATOOLS::DirectoryExists(m_debugDIR_FSR)) {
			ATOOLS::MakeDir(m_debugDIR_FSR);
		}
		m_histograms_FSR["hiddenWeight"] = new Histogram(0, -1, 5, 20);
		if (m_betaorder != 0) {
			m_histograms_FSR["yBet10"] = new Histogram(0, -5, 2, 80);
			m_histograms_FSR["sudY"] = new Histogram(0, -40, 20, 30);
			m_histograms_FSR["sudZ"] = new Histogram(0, -40, 20, 30);
			m_histograms_FSR["sudYFSR"] = new Histogram(0, -40, 20, 30);
			m_histograms_FSR["sudZFSR"] = new Histogram(0, -40, 20, 30);
			m_histograms_FSR["sfac"] = new Histogram(0, -2, 100, 20);
			m_histograms_FSR["hfac"] = new Histogram(0, -2, 100, 20);
			m_histograms_FSR["d10"] = new Histogram(0, -60, 60, 20);
			m_histograms_FSR["cth11"] = new Histogram(0, -1.2, 1.2, 20);
			m_histograms_FSR["cth12"] = new Histogram(0, -1.2, 1.2, 20);
			m_histograms_FSR["cth21"] = new Histogram(0, -1.2, 1.2, 20);
			m_histograms_FSR["cth22"] = new Histogram(0, -1.2, 1.2, 20);
		}
		m_histograms_FSR["beta00"] = new Histogram(0, -2, 6, 20);
		m_histograms_FSR["massWeight"] = new Histogram(0, -1, 5, 20);
		m_histograms_FSR["jacobian"]    = new Histogram(0, 0, 1.2, 20);
		m_histograms_FSR["YFS_IR"]  = new Histogram(0, 0, 5, 20);
		m_histograms_FSR["VolMc"]  = new Histogram(0, 2, 4, 10);
		m_histograms_FSR["wtmass"]  = new Histogram(0, -4, 4, 20);
		m_histograms_FSR["g"]  = new Histogram(0, 5e-2, 7e-2, 20);
		m_histograms_FSR["gp"]  = new Histogram(0, 5e-2, 7e-2, 20);
		m_histograms_FSR["beta"]  = new Histogram(0, 0, 10, 20);
		m_histograms_FSR["massSQ"]  = new Histogram(0, -10, sqrt(m_s) + 10, 50);
		m_histograms_FSR["E1"]  = new Histogram(0, 0.01, sqrt(m_s) / 2. + 10, 20);
		m_histograms_FSR["E2"]  = new Histogram(0, 0.01, sqrt(m_s) / 2. + 10, 20);
		m_histograms_FSR["R_E1"]  = new Histogram(0, -10, sqrt(m_s) / 2. + 10, 20);
		m_histograms_FSR["R_E2"]  = new Histogram(0, -10, sqrt(m_s) / 2. + 10, 20);
		m_histograms_FSR["r1r2"]  = new Histogram(0, 0, sqrt(m_s) + 10, 40);
		m_histograms_FSR["K0"]  = new Histogram(0, 0, sqrt(m_s) + 10, 20);
		m_histograms_FSR["kvec0"]  = new Histogram(0, 0, 15, 20);
		m_histograms_FSR["Eprime"]  = new Histogram(0, 0, sqrt(m_s) + 10, 20);
		m_histograms_FSR["NPhotons"]  = new Histogram(0, -1, 20, 21);
		m_histograms_FSR["NRemoved"]  = new Histogram(0, -1, 20, 21);
		m_histograms_FSR["NPhotons_before_Removal"]  = new Histogram(0, -1, 20, 21);
		m_histograms_FSR["FSRWeight"]  = new Histogram(0, -1, 4, 20);
		m_histograms_FSR["Btil"]  = new Histogram(0, -10, 10, 20);
		m_histograms_FSR["BtilStar"]  = new Histogram(0, -2., 10, 20);
		m_histograms_FSR["BtilQCrude"]  = new Histogram(0, 7, 8, 20);
		m_histograms_FSR["A4"]  = new Histogram(0, -10, 0, 40);
		m_histograms_FSR["A"]  = new Histogram(0, 0, 1, 40);
		m_histograms_FSR["BtilXCrude"]  = new Histogram(0, 7, 8, 20);
		// m_histograms_FSR["Hide-W"] = new Histogram(0,0,2,10);
		m_histograms_FSR["TotalW"]    = new Histogram(0, 0, 2, 10);
		m_histograms_FSR["cos(theta)"]    = new Histogram(0, -1.1, 1.1, 22);
		m_histograms_FSR["sin(theta)"]    = new Histogram(0, -1.1, 1.1, 22);
		m_histograms_FSR["f"]    = new Histogram(0, -10, 10, 20);
		m_histograms_FSR["fbar"]    = new Histogram(0, -10, 50, 20);
		m_histograms_FSR["del1"]    = new Histogram(0, -1, 3, 20);
		m_histograms_FSR["del2"]    = new Histogram(0, -1, 3, 20);
		m_histograms_FSR["DelVol"]    = new Histogram(0, -0.5, 0.5, 20);
		m_histograms_FSR["DelYFS"]    = new Histogram(0, -0.5, 0.5, 20);
		m_histograms_FSR["VoluMC"]    = new Histogram(0, 1, 2, 20);
		m_histograms_FSR["EminQ"]    = new Histogram(0, -10, 10, 40);
		m_histograms_FSR["m_r1"]    = new Histogram(0, 0, 2, 20);
		m_histograms_FSR["m_r2"]    = new Histogram(0, 0, 2, 20);
		m_histograms_FSR["x_r1"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["x_r2"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["x_q1"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["x_q2"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["y_q1"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["y_q2"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["z_q1"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["z_q2"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["x_k"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["y_k"]    = new Histogram(0, -80, 80, 20);
		m_histograms_FSR["delta1"]    = new Histogram(0, 0, 5, 20);
		m_histograms_FSR["BVR_A"]    = new Histogram(0, -10, sqrt(m_s) / 2., 20);
		m_histograms_FSR["pT"]    = new Histogram(0, -10, sqrt(m_s) / 2. + 10, 20);
		m_histograms_FSR["Photon_pT"]    = new Histogram(0, -10, sqrt(m_s) + 10, 20);
		m_histograms_FSR["q1q2"]    = new Histogram(0, -10, sqrt(m_s) + 10, 20);
		m_histograms_FSR["svar"]    = new Histogram(0, -10, sqrt(m_s) + 10, 20);
		m_histograms_FSR["q1_E"]    = new Histogram(0, 0.1, sqrt(m_s) / 2. + 10, 20);
		m_histograms_FSR["q2_E"]    = new Histogram(0, 0.1, sqrt(m_s) / 2. + 10, 20);
		m_histograms_FSR["Q0"]    = new Histogram(0, -10, sqrt(m_s) + 10, 20);
		m_histograms_FSR["phi"]    = new Histogram(0, 0, 2 * 3.14, 10);
		m_histograms_FSR["yy"]  = new Histogram(0, -0.1, 1.1, 20);
		m_histograms_FSR["xfact"]  = new Histogram(0, -0.9, 1.1, 20);
		m_histograms_FSR["Cutflow"]  = new Histogram(0, 0, 5, 5);
		m_histograms_FSR["QQk"]  = new Histogram(0, 0, sqrt(m_s), 20);
		m_histograms_FSR["qboost"]    = new Histogram(0, -10, sqrt(m_s) + 10, 20);
		m_histograms_FSR["pxboost"]    = new Histogram(0, -10, sqrt(m_s) + 10, 20);
		m_histograms_FSR["FSRForm"]    = new Histogram(0, 1., 1.1, 10);

	}
}

void Debug::FillHist(const Vec4D_Vector &plab, YFS::ISR *p_isr, YFS::FSR *p_fsr, double W) {
	if (m_fsr_debug && m_fsrmode != 0) {
		Vec4D_Vector FSRPhotons = p_fsr->p_dipole->GetPhotons();
		Vec4D photonSumFSR = p_fsr->p_dipole->GetPhotonSum();
		double sq = (p_fsr->m_dipole[0] + p_fsr->m_dipole[1]).Abs2();
		double Q0 = (plab[2] + plab[3] + photonSumFSR).E();
		double msq = (plab[2] + plab[3]).Mass();
		double pt = (plab[2] + plab[3]).PPerp();
		double k0(0);
		double ph_pt = photonSumFSR.PPerp();
		// if(FSRPhotons.size()==0) ph_pt = 0;
		int nrem = p_fsr->m_n - FSRPhotons.size();
		for (auto k : FSRPhotons) k0 += k.E();
		double hidecut = 1;
		if (msq > 5) hidecut = 0;
		if (p_fsr->m_fsrWeight != 0 && p_fsr->m_cut == 1) {
			m_histograms_FSR["massSQ"]->Insert(sqrt(p_fsr->m_sQ),W);
			m_histograms_FSR["QQk"]->Insert(sqrt(p_fsr->m_Q * p_fsr->m_photonSum),W);
			m_histograms_FSR["FSRForm"]->Insert(p_fsr->m_fsrform,W);
			m_histograms_FSR["Eprime"]->Insert(sqrt(p_fsr->m_sprim),W);
			m_histograms_FSR["NRemoved"]->Insert(p_fsr->m_NRemoved,W);
			m_histograms_FSR["x_k"]->Insert(photonSumFSR[1],W);
			m_histograms_FSR["y_k"]->Insert(photonSumFSR[2],W);
			m_histograms_FSR["jacobian"]->Insert(p_fsr->m_wt2,W);
			m_histograms_FSR["A4"]->Insert(10 * p_fsr->m_A4,W);
			m_histograms_FSR["A"]->Insert(10 * p_fsr->m_A,W);
			m_histograms_FSR["r1r2"]->Insert(sqrt(p_fsr->m_r1 * p_fsr->m_r2),W);
			m_histograms_FSR["NPhotons"]->Insert(int(FSRPhotons.size()),W);
			m_histograms_FSR["NPhotons_before_Removal"]->Insert(p_fsr->m_n,W);
			m_histograms_FSR["FSRWeight"]->Insert(p_fsr->m_fsrWeight / p_fsr->m_fsrform,W);
			// m_histograms_FSR["E1"]->Insert(p_fsr->m_dipole[0].E(),W);
			// PRINT_VAR(p_dipoles->Old()[0].m_newmomenta[0].E(),W);
			// PRINT_VAR(p_dipoles->Old()[0].m_newmomenta[1].E(),W);
			m_histograms_FSR["E1"]->Insert(p_fsr->m_dipole[0].E(),W);
			m_histograms_FSR["E1"]->Insert(p_fsr->m_dipole[0].E(),W);
			m_histograms_FSR["E2"]->Insert(p_fsr->m_dipole[1].E(),W);
			m_histograms_FSR["R_E1"]->Insert(p_fsr->m_r1.E(),W);
			m_histograms_FSR["R_E2"]->Insert(p_fsr->m_r2.E(),W);
			m_histograms_FSR["x_r1"]->Insert(p_fsr->m_r1[1],W);
			m_histograms_FSR["x_r2"]->Insert(p_fsr->m_r2[1],W);
			m_histograms_FSR["x_q1"]->Insert(plab[2][1],W);
			m_histograms_FSR["x_q2"]->Insert(plab[3][1],W);
			m_histograms_FSR["y_q1"]->Insert(plab[2][2],W);
			m_histograms_FSR["y_q2"]->Insert(plab[3][2],W);
			m_histograms_FSR["z_q1"]->Insert(plab[2][3],W);
			m_histograms_FSR["z_q2"]->Insert(plab[3][3],W);
			m_histograms_FSR["delta1"]->Insert(p_fsr->m_delta1 * 1e7,W);
			m_histograms_FSR["BVR_A"]->Insert(p_fsr->m_bvrA * 1e3,W);
			m_histograms_FSR["pT"]->Insert(pt,W);
			m_histograms_FSR["q1q2"]->Insert(sqrt(p_fsr->m_q1q2),W);
			m_histograms_FSR["svar"]->Insert(sqrt(p_fsr->m_sp),W);
			m_histograms_FSR["q1_E"]->Insert(plab[2].E(),W);
			m_histograms_FSR["q2_E"]->Insert(plab[3].E(),W);
			m_histograms_FSR["Q0"]->Insert(Q0,W);
			std::vector<double> cos   = p_fsr->m_cos;
			std::vector<double> sin   = p_fsr->m_sin;
			std::vector<double> phi   = p_fsr->m_phi_vec;
			std::vector<double> kvec0   = p_fsr->m_k0;
			std::vector<double> dist1 = p_fsr->m_dist1;
			std::vector<double> dist2 = p_fsr->m_dist2;
			std::vector<double> del1 = p_fsr->m_del1;
			std::vector<double> del2 = p_fsr->m_del2;
			for (auto c : cos)   m_histograms_FSR["cos(theta)"]->Insert(c,W);
			for (auto s : sin)   m_histograms_FSR["sin(theta)"]->Insert(s,W);
			for (auto p : phi)   m_histograms_FSR["phi"]->Insert(p,W);
			for (auto p : kvec0)   m_histograms_FSR["kvec0"]->Insert(p * p_fsr->m_xfact,W);
			for (auto d : dist1) m_histograms_FSR["f"]->Insert(d / 1e4,W);
			for (auto fbar : dist2) m_histograms_FSR["fbar"]->Insert(fbar / 1e4,W);
			for (auto d1 : del1) m_histograms_FSR["del1"]->Insert(d1,W);
			for (auto d2 : del2) m_histograms_FSR["del2"]->Insert(d2,W);
			m_histograms_FSR["Btil"]->Insert(p_fsr->m_btil,W);
			m_histograms_FSR["beta00"]->Insert(log(m_born),W);
			m_histograms_FSR["BtilQCrude"]->Insert(p_fsr->m_BtiQcru,W);
			m_histograms_FSR["BtilStar"]->Insert(p_fsr->m_btilStar,W);
			m_histograms_FSR["BtilXCrude"]->Insert(p_fsr->m_BtiXcru,W);
			m_histograms_FSR["VoluMC"]->Insert(p_fsr->m_volumc,W);
			m_histograms_FSR["DelVol"]->Insert(p_fsr->m_delvol,W);
			m_histograms_FSR["DelYFS"]->Insert(p_fsr->m_DelYFS,W);
			m_histograms_FSR["EminQ"]->Insert(log(p_fsr->m_EminQ * 1e6),W);
			m_histograms_FSR["m_r1"]->Insert(p_fsr->m_r1.Mass() * 10,W);
			m_histograms_FSR["m_r2"]->Insert(p_fsr->m_r2.Mass() * 10,W);
			m_histograms_FSR["YFS_IR"]->Insert(p_fsr->m_YFS_IR,W);
			m_histograms_FSR["VolMc"]->Insert(exp(p_fsr->m_volmc),W);
			m_histograms_FSR["wtmass"]->Insert(p_fsr->m_massW,W);
			m_histograms_FSR["hiddenWeight"]->Insert(p_fsr->m_hideW,W);
			m_histograms_FSR["massWeight"]->Insert(p_fsr->m_massW,W);
		}
		m_histograms_FSR["Cutflow"]->Insert(p_fsr->m_cut,W);
		m_histograms_FSR["g"]->Insert(p_fsr->m_g,W);
		m_histograms_FSR["gp"]->Insert(p_fsr->m_gp,W);
		m_histograms_FSR["beta"]->Insert(p_fsr->m_beta1,W);
		// m_histograms_FSR["qboost"]->Insert(m_Q.Mass(),W);
		// m_histograms_FSR["pxboost"]->Insert(m_Q.Mass(),W);
		if (FSRPhotons.size() != 0) {
			m_histograms_FSR["K0"]->Insert(photonSumFSR.E(),W);
			m_histograms_FSR["Photon_pT"]->Insert(p_fsr->m_photonSumPreBoost.PPerp(),W);
			m_histograms_FSR["yy"]->Insert(p_fsr->m_yy,W);
			m_histograms_FSR["xfact"]->Insert(1. / p_fsr->m_xfact,W);
		}
		// }
		// else m_histograms_FSR["Cutflow"]->Insert(0,W);
	}
	if (m_isr_debug && m_fsrmode != 2) {
		Vec4D_Vector ISRPhotons = p_isr->GetPhotons();
		m_histograms_ISR["v"]->Insert(p_isr->m_v,W);
		m_histograms_ISR["K0"]->Insert(p_isr->m_photonSum[0],W);
		m_histograms_ISR["NGamma"]->Insert(int(ISRPhotons.size()),W);
		m_histograms_ISR["nbar"]->Insert(p_isr->m_nbar,W);
		m_histograms_ISR["weight"]->Insert(p_isr->m_weight,W);
		m_histograms_ISR["massweight"]->Insert(p_isr->m_massW,W);
		m_histograms_ISR["jacweight"]->Insert(p_isr->m_jacW,W);
		// m_histograms_ISR["Cutflow"]->Insert(p_isr->m_cut,W);
		// PRINT_VAR(p_isr->m_cut);
		if (p_isr->m_cut == 0) m_histograms_ISR["Cutflow"]->Insert(0.1,W);
		else m_histograms_ISR["Cutflow"]->Insert(1.1,W);
		std::vector<double> cos   = p_isr->m_cos;
		std::vector<double> jac   = p_isr->m_jacvec;
		std::vector<double> scale = p_isr->m_scale;
		std::vector<double> AA = p_isr->m_AA;
		std::vector<double> K2 = p_isr->m_K2;
		std::vector<double> ptk = p_isr->m_PTK;
		for (auto c : cos)   m_histograms_ISR["cos(theta)"]->Insert(c,W);
		for (auto j : jac)   m_histograms_ISR["jacobian"]->Insert(j,W);
		for (auto s : scale) m_histograms_ISR["lam"]->Insert(s,W);
		for (auto A : AA) m_histograms_ISR["A"]->Insert(A,W);
		for (auto k : K2) m_histograms_ISR["K2"]->Insert(k,W);
		for (auto p : ptk) m_histograms_ISR["ptk"]->Insert(p,W);
	}
}

void Debug::FillHist(const std::string &name, const double &x, double weight){
	std::map<string, Histogram *>::iterator itisr  = m_histograms_ISR.find(name);
	if(m_fsr_debug) std::map<string, Histogram *>::iterator itfsr  = m_histograms_FSR.find(name);
	if(m_isr_debug){
		if(itisr != m_histograms_ISR.end()) {
			m_histograms_ISR[name]->Insert(x,weight);
			return;
		}
	}
	if(m_fsr_debug){
		if(itisr != m_histograms_FSR.end()){
			m_histograms_FSR[name]->Insert(x,weight);
			return;
		}
	}
	THROW(fatal_error, "Histogram with key: "+name+" not found in YFS Debug");
}


void Debug::FillHist(const std::string &name, const double &x, const double &y, double weight){
	std::map<string, Histogram_2D *>::iterator it  = m_histograms2d.find(name);
	if(it != m_histograms2d.end()) m_histograms2d[name]->Insert(x,y,weight);
	else THROW(fatal_error, "Histogram2D with key: "+name+" not found in YFS Debug");

}


void Debug::WriteHistograms() {
	Histogram * histo;
	string name;
	if (m_fsr_debug && m_fsrmode != 0) {
		for (map<string, Histogram *>::iterator hit = m_histograms_FSR.begin();
		        hit != m_histograms_FSR.end(); hit++) {
			histo = hit->second;
			name  = string("./" + m_debugDIR_FSR + "/") + hit->first + string(".dat");
			histo->MPISync();
			histo->Finalize();
			histo->Output(name);
			delete histo;
		}
	}
	if (m_isr_debug && m_fsrmode != 2) {
		for (map<string, Histogram *>::iterator hit = m_histograms_ISR.begin();
		        hit != m_histograms_ISR.end(); hit++) {
			histo = hit->second;
			PRINT_VAR(hit->first);
			name  = string("./" + m_debugDIR_ISR + "/") + hit->first + string(".dat");
			histo->MPISync();
			histo->Finalize();
			histo->Output(name);
			delete histo;
		}
		Histogram_2D * histo2d;
		string name;
		for (map<string, Histogram_2D *>::iterator hit = m_histograms2d.begin();
		        hit != m_histograms2d.end(); hit++) {
			histo2d = hit->second;
			name  = string(m_debugDIR_ISR) + "/" + hit->first + string(".dat");
			// histo2d->MPISync();
			histo2d->Finalize();
			histo2d->Output(name);
			delete histo2d;
		}
	}
}
