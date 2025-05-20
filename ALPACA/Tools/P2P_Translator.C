#include "ALPACA/Tools/P2P_Translator.H"
#include "ATOOLS/Org/Message.H"

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;


P2P_Translator::P2P_Translator() {}

shared_ptr<Parton> P2P_Translator::operator()(Particle * particle){
    double inittau = 0.;
    Vec4D p = particle->Momentum();
    //p[0] = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]); //Put on-shell for now
    shared_ptr<Parton> parton = make_shared<Parton>(
        particle->Flav(),
        p,
        particle->XProd(),
        inittau
    );
    parton->SetNumber();
    parton->SetFlow(1,particle->GetFlow(1));
    parton->SetFlow(2,particle->GetFlow(2));
    parton->AddXPHistory(0., particle->XProd(), p, particle->Flav(), 0.);
    return parton;
}

Particle * P2P_Translator::operator()(shared_ptr<Parton> parton){
    Particle * particle = new Particle(-1,parton->Flav(),
				       parton->Momentum());
    particle->SetNumber();
    particle->SetInfo('F');
    particle->SetFlow(1,parton->GetFlow(1));
    particle->SetFlow(2,parton->GetFlow(2));
    particle->SetPosition(parton->Position());
    return particle;
}