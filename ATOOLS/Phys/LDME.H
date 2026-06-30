#ifndef LDME_H
#define LDME_H

#include "ATOOLS/Org/Run_Parameter.H"
#include <unordered_map>
#include <list> 

struct decay_channels{
  std::list<kf_code> decays;
  std::string name;
  kf_code kfc;
};

void msgLDME(kf_code kfc, double LDME);
void LoadLDME();
void LoadLDMEmap();
const std::map<kf_code, double>& GetLDMEmap(kf_code kfc);
void WarnZeroLDME(kf_code kfc);
double GetLDME(kf_code kfc);
kf_code GetQuarkoniumDecayChannel(kf_code octet);


extern std::unordered_map<kf_code, double> ldme__sin_map;
extern std::map<kf_code, std::map<kf_code, double>> ldme_oct_map;

extern decay_channels d_1S0_c, d_3S1_c, d_3P0_c, d_3P1_c, d_3P2_c;
extern std::vector<decay_channels> d_alldecays_charm;


#endif