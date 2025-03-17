#include "ALPACA/Tools/P2P_Translator.H"
#include "ATOOLS/Org/Message.H"

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


P2P_Translator::P2P_Translator() {}

Parton * P2P_Translator::operator()(Particle * particle){
    double inittau = 0.;
    Parton * parton = new Parton(particle->Flav(),
				 particle->Momentum(),
                 particle->XProd(), inittau);
    parton->SetFlow(1,particle->GetFlow(1));
    parton->SetFlow(2,particle->GetFlow(2));
    return parton;
}

Particle * P2P_Translator::operator()(Parton * parton){
    Particle * particle = new Particle(-1,parton->Flav(),
				       parton->Momentum());
    particle->SetNumber();
    particle->SetInfo('F');
    particle->SetFlow(1,parton->GetFlow(1));
    particle->SetFlow(2,parton->GetFlow(2));
    particle->SetPosition(parton->Position());
    return particle;
}
