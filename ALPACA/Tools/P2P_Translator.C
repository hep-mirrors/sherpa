#include "ALPACA/Tools/P2P_Translator.H"
#include "ATOOLS/Org/Message.H"

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


P2P_Translator::P2P_Translator() {}

Parton P2P_Translator::operator()(Particle * particle, Vec4D pos){
    double inittau = 0.;
    Parton parton(particle->Flav(), particle->Momentum(), pos, inittau);
    return parton;
}

Particle P2P_Translator::operator()(Parton * parton){
    Particle particle(-1,parton->Flav(), parton->Momentum());
    return particle;
}
