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
std::map<kf_code, double> GetLDMEmap(kf_code kfc);
double GetTotalLDME(kf_code kfc);
void WarnZeroLDME(kf_code kfc);
double GetLDME(kf_code kfc);
kf_code GetChannel(kf_code octet);


extern std::map<kf_code, double> singlets;
extern std::unordered_map<kf_code, double> ldme_map;
extern std::map<kf_code, std::map<kf_code, double>> LDMEmap;

extern decay_channels d_1S0_c, d_3S1_c, d_3P0_c, d_3P1_c, d_3P2_c;
extern std::list<decay_channels> d_ch_c;


#endif