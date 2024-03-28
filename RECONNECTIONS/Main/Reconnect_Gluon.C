#include "RECONNECTIONS/Main/Reconnect_Gluon.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Histogram.H"
//#include <string>

using namespace RECONNECTIONS;

Reconnect_Gluon::Reconnect_Gluon() : Reconnection_Base()
{
	//need to determine correct parameters to pass to Histogram ...
	//this->m_histomap[std::string("total_stringlength")] = new ATOOLS::Histogram();
	this->m_stringlength_totals.clear();
}

Reconnect_Gluon::~Reconnect_Gluon()
{
	m_collist.clear();
//	msg_Info() << std::endl << "Min. total stringlength: " << this->find_min_totalLength(this->m_stringlength_totals) << std::endl;
//	msg_Info() << std::endl << "Max. total stringlength: " << this->find_max_totalLength(this->m_stringlength_totals) << std::endl;
	write_stringLengths_to_file(m_stringlength_totals);
	m_stringlength_totals.clear();
}

void Reconnect_Gluon::SetParameters()
{
	auto s       = ATOOLS::Settings::GetMainSettings()["COLOUR_RECONNECTIONS"];
	m_Pmode      = s["PMODE"].SetDefault(0).Get<int>();
	m_Q02        = std::sqrt(s["Q_0"].SetDefault(1.00).Get<double>());
	m_etaQ       = std::sqrt(s["eta0"].SetDefault(0.1).Get<double>());
	m_R02        = std::sqrt(s["R_0"].SetDefault(100.).Get<double>());
	m_etaR       = std::sqrt(s["etaR"].SetDefault(0.16).Get<double>());
	m_reshuffle  = s["Reshuffle"].SetDefault(1./9.).Get<double>();
	m_kappa      = s["kappa"].SetDefault(2.).Get<double>();
	//m_mom_dist   = s["mom_dist"].SetDefault(true).Get<bool>();
	//m_pos_dist   = s["pos_dist"].SetDefault(false).Get<bool>();
	//m_col_dist   = s["col_dist"].SetDefault(false).Get<bool>();
	//m_dist_type  = s["dist_type"].SetDefault("mom").Get<std::string>();
//	m_dist_type  = s["dist_type"].SetDefault("mom").Get<string>();
	m_dist_type  = s["dist_type"].SetDefault(1).Get<int>();
	// add new field in runcard settings to specify type of CR to be performed
}

void Reconnect_Gluon::Reset()
{
	m_collist.clear();
	m_gluon_collist.clear();
	Reconnection_Base::Reset();
}

int Reconnect_Gluon::operator()(ATOOLS::Blob_List *const blobs)
{
	if(!HarvestParticles(blobs)) 		   return -1;
	if(m_cols[0].empty() && m_cols[1].empty()) return 0;
	
	// inefficient to call twice, come back and store the result somewhere...
	m_norm = this->TotalLength();
	this->PlotTotalLength(this->TotalLength());

	for(std::map<unsigned int, ATOOLS::Particle *>::iterator cit = m_cols[0].begin(); cit != m_cols[0].end(); cit++)
	{
		m_collist.push_back(cit->first);
		// if gluon, populate gluon list. Use this list downstream 	
		if(cit->second->Flav().IsGluon()) m_gluon_collist.push_back(cit->first);
	}
	std::size_t N = m_gluon_collist.size();
	unsigned int col[2];
	for(std::size_t i=0; i<std::sqrt(N); i++)
	{
		// store copy of original col pair here for before mass plot ...
		if(!this->SelectColourPair(N, col[0], col[1])) break;
		auto copy_cols = m_cols;
		if(!this->AttemptSwap(col)) return false;
		else
		{
			FillMassesInHistogram(copy_cols, "Mass_beforeGluonMove");
			FillMassesInHistogram(m_cols, "Mass_afterGluonMove");
		}
	}
	this->UpdateColours();
	m_collist.clear();
	
	return true;
}


bool Reconnect_Gluon::SelectColourPair(const size_t &N, unsigned int &g1, unsigned int &g2)
{
	msg_Info() << "\n Reconnect_Gluon::SelectColourPair() called...\n";
	unsigned int trials = 0;
	do
	{
		g1 = m_gluon_collist[int(ATOOLS::ran->Get()*N)];
		g2 = m_gluon_collist[int(ATOOLS::ran->Get()*N)];
		if((trials++)==100)
		{
			g1 = g2 = 0;
			return false;
		}
	// does the logic below still hold? Should I have a list of m_cols only gluons instead?
	// index in gluon_collist does not match index in collist, g1/2 is essentially the key/index, so need to do via 'find' rather than via index.
// but this is m_cols, so why do we care?
// this looks fine as it is...
	} while(g1 == g2 || m_cols[0][g1] == m_cols[1][g2] || m_cols[0][g2] == m_cols[1][g1]);
	return true;
}

// attempt swap method ...
bool Reconnect_Gluon::AttemptSwap(const unsigned int gluons[2])
{
	msg_Info() << "\n Reconnect_Gluon::AttemptSwap() called...\n";
	if (m_cols[0].find(gluons[0]) == m_cols[0].end() || 
	    m_cols[0].find(gluons[1]) == m_cols[0].end() ||
	    m_cols[1].find(gluons[0]) == m_cols[1].end() ||
	    m_cols[1].find(gluons[1]) == m_cols[1].end()  )
	{
		msg_Error() << "Error in " << METHOD << ": ill-defined colours. \n";
		return false;
	}
	
	ATOOLS::Particle *part[4];
	for (std::size_t i=0; i<2; i++)
	{
		for(std::size_t j=0; j<2; j++) part[2*i+j] = m_cols[i][gluons[j]];
	}

	double dist0  = Distance(part[0], part[2]), dist1  = Distance(part[1], part[3]);
	double ndist0 = Distance(part[0], part[3]), ndist1 = Distance(part[1], part[3]);
	// take probabilistic approach for moment
	// later, try deterministic approach by computing total stringlength before and after swap, and keep swap if 
	// it reduces.
	double prob = m_reshuffle * std::exp(-m_etaQ * ((ndist0+ndist1)-(dist0+dist1)));
	if (prob > ATOOLS::ran->Get())
	{
		m_cols[1][gluons[0]] = part[3];
		m_cols[1][gluons[1]] = part[2];
	}
	
	return true;
}

void Reconnect_Gluon::UpdateColours()
{
	msg_Info() << "Reconnect_Gluon::UpdateColours() called...";
	for(std::size_t i=0; i<2; i++)
	{
		for(std::map<unsigned int, ATOOLS::Particle*>::iterator cit = m_cols[i].begin();
		   cit != m_cols[i].end(); cit++)
		{
			cit->second->SetFlow(i+1, cit->first);
		}
	}

}

double Reconnect_Gluon::Distance(ATOOLS::Particle *trip, ATOOLS::Particle *anti)
{
	// can't switch on a string
	/*switch(m_dist_type)
	{
		case "mom":
			return MomDistance(trip, anti);
			break; // don't really neee presumably, because return statement
		case "pos":
			return PosDistance(trip, anti);
			break;
		case "col":
			return ColDistance(trip, anti);
			break;

		default:
			return MomDistance(trip, anti);
	}*/
	
	/*if     (m_dist_type == "mom") return MomDistance(trip, anti);
	else if(m_dist_type == "pos") return PosDistance(trip, anti);
	else if(m_dist_type == "col") return ColDistance(trip, anti);
	else			      return MomDistance(trip, anti);*/

	switch(m_dist_type)
        {
                case 1:
                        return MomDistance(trip, anti);
                        break; // don't really neee presumably, because return statement
                case 2:
                        return PosDistance(trip, anti);
                        break;
                case 3:
                        return ColDistance(trip, anti);
                        break;

                default:
                        return MomDistance(trip, anti);
        }
}

double Reconnect_Gluon::MomDistance(ATOOLS::Particle *trip, ATOOLS::Particle *anti)
{
	double p1p2 = trip->Momentum() * anti->Momentum();
	if(trip->Flav().IsGluon()) p1p2 /= 2.;
	if(anti->Flav().IsGluon()) p1p2 /= 2.;
	double m1m2 = trip->Flav().HadMass() * anti->Flav().HadMass();
	return p1p2-m1m2;
}

double Reconnect_Gluon::PosDistance(ATOOLS::Particle *trip, ATOOLS::Particle *anti)
{
	double xdist2 = std::abs((trip->XProd()-anti->XProd()).Abs2());
	return xdist2<m_R02 ? 1. : std::pow(xdist2/m_R02, m_etaR);
}

double Reconnect_Gluon::ColDistance(ATOOLS::Particle *trip, ATOOLS::Particle *anti)
{
	return trip->GetFlow(1)==anti->GetFlow(2) ? 1. : m_reshuffle;
}

// add different distance metrics ...

double Reconnect_Gluon::TotalLength()
{
	double total = 1.;
	ATOOLS::Particle *part1, *part2;
	for(std::map<unsigned int, ATOOLS::Particle*>::iterator cit=m_cols[0].begin();
	   cit!=m_cols[0].end(); cit++)
	{
		part1 = cit->second;
		part2 = m_cols[1].find(cit->first)->second;
		total *= Distance(part1, part2);
	}

	return total / std::pow(m_parts[0].size(), m_kappa);
}

void Reconnect_Gluon::PlotTotalLength(double total_length)
{
	m_stringlength_totals.push_back(total_length);
	// plot this ...
}

void Reconnect_Gluon::write_stringLengths_to_file(std::vector<double> total_vec)
{
	std::ofstream output_file("./Reconnection_Analysis/gluonOnly_stringlength_vec.txt");
	std::ostream_iterator<double> output_iterator(output_file, "\n");
	std::copy(total_vec.begin(), total_vec.end(), output_iterator);
}

// -------

/*
 * idea: take one coloured particle, if it != gluon: skip
 * 	 iterate over each possible alternative partner particle:
 *		check it's a gluon. else: skip
 *		if both gluons: compute new stringlength
 * 		if > previous string-length: discard
 *		... keep swaps that reduce string length
 *		... do we want to keep iterating until we find the minimum stringlength?
 * break down into new functions as appropriate ...
 */

// draw out each data structure and label ...
