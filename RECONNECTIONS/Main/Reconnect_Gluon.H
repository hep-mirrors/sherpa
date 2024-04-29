#ifndef RECONNECTIONS_Main_Reconnect_Gluon_H
#define RECONNECTIONS_Main_Reconnect_Gluon_H

#include "RECONNECTIONS/Main/Reconnection_Base.H"
//#include <string>

namespace RECONNECTIONS
{
	class Reconnect_Gluon : public Reconnection_Base
	{
		private:
			std::vector<unsigned int> m_collist;
			std::vector<unsigned int> m_gluon_collist;
			int m_Pmode;
			double m_Q02, m_etaQ, m_R02, m_etaR, m_reshuffle, m_kappa, m_norm;
			//bool m_mom_dist, m_pos_dist, m_col_dist;
			//std::string m_dist_type; // could define in parent class instead so that the handler can include dist type in filename write
			int m_dist_type;
			std::vector<double> m_stringlength_totals;
			std::vector<double> m_hadmasses_before, m_hadmasses_after;
			
			void SetParameters();
			//bool SelectGluonList();
			bool SelectColourPair(const size_t &N, unsigned int &col1, unsigned int &col2);
			bool AttemptSwap(const unsigned int col[2]);
			void UpdateColours();

			double MomDistance(ATOOLS::Particle *trip, ATOOLS::Particle *anti);
			double PosDistance(ATOOLS::Particle *trip, ATOOLS::Particle *anti);
			double ColDistance(ATOOLS::Particle *trip, ATOOLS::Particle *anti);
			double TotalLength();
			void   save_total_length(double total_length);
//			double find_min_totalLength(std::vector<double> total_vec);
//			double find_max_totalLength(std::vector<double> total_vec);
			void   write_stringLengths_to_file(std::vector<double> total_vec);
			void   write_to_file(std::vector<double> data, std::string file_name);
			void   save_masses_before(std::map<unsigned int, ATOOLS::Particle*> cols_before[2]);
			void   save_masses_after(std::map<unsigned int, ATOOLS::Particle*> cols_after[2]);
		
		public:
			Reconnect_Gluon();
			~Reconnect_Gluon();
			
			int    operator()(ATOOLS::Blob_List *const blobs);
			void   Reset();
			double Distance(ATOOLS::Particle *trip, ATOOLS::Particle *anti);
	};
}
#endif
