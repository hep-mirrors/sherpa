#include "ALPACA/EventGeneration/Colour_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

#include <list>

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;

Colour_Handler::Colour_Handler(shared_ptr<list<shared_ptr<Parton>>> ptr_partons, 
                               shared_ptr<list<pair<int,int>>>      ptr_flow_backtrack):
p_partons(ptr_partons), p_flow_backtrack(ptr_flow_backtrack)
{
}

Colour_Handler::~Colour_Handler()
{
}

void Colour_Handler::UpdateColoursScatter(shared_ptr<Parton> part_in_1, shared_ptr<Parton> part_in_2, 
                                          Flavour flav_out_1, Flavour flav_out_2, 
                                          double t, double mg2, double mf2){
    Flavour flav_in_1 = part_in_1->Flav();
    Flavour flav_in_2 = part_in_2->Flav();
    bool found_process = false;

    if(flav_in_1.IsQuark() && flav_in_2.IsQuark()) {
        if (flav_in_1 == flav_in_2) {
            found_process = UpdateColours_qiqi_qiqi(part_in_1, part_in_2, flav_out_1, flav_out_2, t, mg2);
        } else if(flav_in_1 == flav_in_2.Bar()) {
            if(flav_out_1.IsGluon() && flav_out_2.IsGluon()){
                found_process = UpdateColours_qiqbi_gg(part_in_1, part_in_2, flav_out_1, flav_out_2, t, mg2, mf2);
            } else if(flav_in_1 == flav_out_1 || flav_in_1 == flav_out_2){
                found_process = UpdateColours_qiqbi_qiqbi(part_in_1, part_in_2, flav_out_1, flav_out_2, t, mg2);
            } else{
                found_process = UpdateColours_qiqbi_qjqbj(part_in_1, part_in_2, flav_out_1, flav_out_2);
            }
        } else if(flav_in_1 != flav_in_2) {
            found_process = UpdateColours_qiqj_qiqj(part_in_1, part_in_2, flav_out_1, flav_out_2);
        }
    } else if((flav_in_1.IsGluon() && flav_in_2.IsQuark()) || (flav_in_1.IsQuark() && flav_in_2.IsGluon())) {
        found_process = UpdateColours_qig_qig(part_in_1, part_in_2, flav_out_1, flav_out_2, t, mg2, mf2);
    } else if(flav_in_1.IsGluon() && flav_in_2.IsGluon()) {
        if(flav_out_1.IsGluon() && flav_out_2.IsGluon()){
            found_process = UpdateColours_gg_gg(part_in_1, part_in_2, flav_out_1, flav_out_2, t, mg2);
        } else{
            found_process = UpdateColours_gg_qiqbi(part_in_1, part_in_2, flav_out_1, flav_out_2, t, mg2, mf2);
        }  
    }

    if(!found_process){
        msg_Out() << METHOD << ": ERROR: process not found when updating colors for scattering, will exit." << endl;
        msg_Out() << METHOD << ": " << flav_in_1 << " " << flav_in_2 << " -> " << flav_out_1 << " " << flav_out_2 << endl;
        exit(1.); 
    }
}

void Colour_Handler::UpdateColoursSplit(shared_ptr<Parton> part_in,
                                        shared_ptr<Parton> part_out_1, shared_ptr<Parton> part_out_2){
    Flavour flav_in = part_in->Flav();
    Flavour flav_out_1 = part_out_1->Flav();
    Flavour flav_out_2 = part_out_2->Flav();
    bool found_process = false;

    if(flav_in.IsQuark()){
        // ####  q->gq ####
        if(flav_in.IsAnti()){ //qb -> 
            if(flav_out_1.IsGluon()){ // qb -> g qb
                part_out_1->SetFlow(1, -1);
                part_out_1->SetFlow(2, part_in->GetFlow(2));
                part_out_2->SetFlow(2, part_out_1->GetFlow(1));
                found_process = true;
            } else{ // qb -> qb g
                part_out_2->SetFlow(1, -1);
                part_out_2->SetFlow(2, part_in->GetFlow(2));
                part_out_1->SetFlow(2, part_out_2->GetFlow(1));
                found_process = true;
            }
        } else{ // q ->
            if(flav_out_1.IsGluon()){ // q -> g q
                part_out_1->SetFlow(1, part_in->GetFlow(1));
                part_out_1->SetFlow(2, -1);
                part_out_2->SetFlow(1, part_out_1->GetFlow(2));
                found_process = true;
            } else{ // q -> q g
                part_out_2->SetFlow(1, part_in->GetFlow(1));
                part_out_2->SetFlow(2, -1);
                part_out_1->SetFlow(1, part_out_2->GetFlow(2));
                found_process = true;
            }
        }
    } else if(flav_out_1.IsGluon() && flav_out_2.IsGluon()){
        // ####  g->gg ####
        part_out_1->SetFlow(1, part_in->GetFlow(1));
        part_out_1->SetFlow(2, -1);
        part_out_2->SetFlow(1, part_out_1->GetFlow(2));
        part_out_2->SetFlow(2, part_in->GetFlow(2));
        found_process = true;
    } else if(flav_out_1.IsQuark() && flav_out_2.IsQuark()){
        // ####  g->qqb ####
        if(flav_out_1.IsAnti()){ // g -> qb q
            part_out_1->SetFlow(2, part_in->GetFlow(2));
            part_out_2->SetFlow(1, part_in->GetFlow(1));
            found_process = true;
        } else{ // g -> q qb
            part_out_1->SetFlow(1, part_in->GetFlow(1));
            part_out_2->SetFlow(2, part_in->GetFlow(2));
            found_process = true;
        }
    }

    if(!found_process){
        msg_Out() << METHOD << ": ERROR: process not found when updating colors for splitting, will exit." << endl;
        msg_Out() << METHOD << ": " << flav_in << " -> " << flav_out_1 << " " << flav_out_2 << endl;
        exit(1.); 
    }
}

void Colour_Handler::UpdateColoursMerge(shared_ptr<Parton> part_in_1, shared_ptr<Parton> part_in_2, 
                                        shared_ptr<Parton> part_out){
    Flavour flav_in_1 = part_in_1->Flav();
    Flavour flav_in_2 = part_in_2->Flav();
    Flavour flav_out = part_out->Flav();
    bool found_process = false;
    int i_1, i_2; //Indices to backtrack

    if(flav_out.IsQuark()){
        // ####  gq -> q ####
        if(flav_out.IsAnti()){ // -> qb
            if(flav_in_1.IsGluon()){ // g qb -> qb
                part_out->SetFlow(2, part_in_1->GetFlow(2));
                i_1 = part_in_1->GetFlow(1);
                i_2 = part_in_2->GetFlow(2);
                found_process = true;
            } else{ //qb g -> qb
                part_out->SetFlow(2, part_in_2->GetFlow(2));
                i_1 = part_in_2->GetFlow(1);
                i_2 = part_in_1->GetFlow(2);
                found_process = true;
            }
        } else{ // -> q
            if(flav_in_1.IsGluon()){ // g q -> q
                part_out->SetFlow(1, part_in_1->GetFlow(1));
                i_1 = part_in_1->GetFlow(2);
                i_2 = part_in_2->GetFlow(1);
                found_process = true;
            } else{ //q g -> q
                part_out->SetFlow(1, part_in_2->GetFlow(1));
                i_1 = part_in_2->GetFlow(2);
                i_2 = part_in_1->GetFlow(1);
                found_process = true;
            }
        }
        //Backtrack colors
        if(found_process){
            BackTrackColour(i_1, i_2, part_in_1, part_in_2);
        }
    } else if(flav_in_1.IsQuark() && flav_in_2.IsQuark()){
        // ####  qqb -> g ####
        if(flav_in_1.IsAnti()){ //qb q -> g
            if(part_in_1->GetFlow(2) == part_in_2->GetFlow(1)){
                InfSoftQQbtoG(part_in_2,part_in_1);
            }
            part_out->SetFlow(1, part_in_2->GetFlow(1));
            part_out->SetFlow(2, part_in_1->GetFlow(2));
            found_process = true;
        } else{ // q qb -> g
            if(part_in_1->GetFlow(1) == part_in_2->GetFlow(2)){
                InfSoftQQbtoG(part_in_1,part_in_2);
            }
            part_out->SetFlow(1, part_in_1->GetFlow(1));
            part_out->SetFlow(2, part_in_2->GetFlow(2));
            found_process = true;
        }
    } else if(flav_in_1.IsGluon() && flav_in_2.IsGluon()){
        // ####  gg -> g ####
        if(part_in_1->GetFlow(1) == part_in_2->GetFlow(2) && part_in_1->GetFlow(2) == part_in_2->GetFlow(1)){
            //If all indices match up to create a colour neutral gluon:
            //Exchange an infinitely soft gluon before for part_in_1 to change to colour, then merge
            InfSoftGGtoG(part_in_1,part_in_2);
            part_out->SetFlow(1, part_in_2->GetFlow(1));
            part_out->SetFlow(2, part_in_1->GetFlow(2));
            BackTrackColour(part_in_1->GetFlow(1), part_in_2->GetFlow(2), part_in_1, part_in_2);
        } else if(part_in_1->GetFlow(1) == part_in_2->GetFlow(2)){
            part_out->SetFlow(1, part_in_2->GetFlow(1));
            part_out->SetFlow(2, part_in_1->GetFlow(2));
            BackTrackColour(part_in_1->GetFlow(1), part_in_2->GetFlow(2), part_in_1, part_in_2);
        } else{
            part_out->SetFlow(1, part_in_1->GetFlow(1));
            part_out->SetFlow(2, part_in_2->GetFlow(2));
            BackTrackColour(part_in_1->GetFlow(2), part_in_2->GetFlow(1), part_in_1, part_in_2);
        }
        
         found_process = true;
    }

    if(!found_process){
        msg_Out() << METHOD << ": ERROR: process not found when updating colors for merging, will exit." << endl;
        msg_Out() << METHOD << ": " << flav_in_1 << " " << flav_in_2 << " -> " << flav_out << endl;
        exit(1.); 
    }
}

void Colour_Handler::BackTrackColour(int flow_id_1, int flow_id_2, 
                                     shared_ptr<Parton> part_exc_1, shared_ptr<Parton> part_exc_2){
    // Backtrack the color flow of two indicies which are removed in a scattering, splitting or merging.
    // Sets flow_id_1 = flow_id_2 for the first instance of flow_id_1 or flow_id_2 found in the
    // current list of active partons. Does not affect part_exc_1 or part_exc_2.
                                    
    bool found_index = false;

    if(flow_id_1 == flow_id_2){
        //Already correct colour-pair, no need to backtrack
        found_index = true;
    }

    if(!found_index){
        //If not already colour pair, look through all partons in ALPACA to backtrack
        for (list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
            if ((*iter) == part_exc_1 || (*iter) == part_exc_2) continue;

            if((*iter)->GetFlow(1) == flow_id_1){
                (*iter)->SetFlow(1, flow_id_2);
                found_index = true;
            } else if((*iter)->GetFlow(2) == flow_id_1){
                (*iter)->SetFlow(2, flow_id_2);
                found_index = true;
            } else if((*iter)->GetFlow(1) == flow_id_2){
                (*iter)->SetFlow(1, flow_id_1);
                found_index = true;
            } else if((*iter)->GetFlow(2) == flow_id_2){
                (*iter)->SetFlow(2, flow_id_1);
                found_index = true;
            }

            //if(found_index) msg_Out() << flow_id_1 << " set to " << flow_id_2 << endl;
            if(found_index) break;
        }
    }

    if(!found_index){
        //If not correct colour pair and colour partner not found in ALPACA, will exist in external blob
        msg_Out() << METHOD << ": (Temporary): color indices " << flow_id_1 << " and " << flow_id_1 << "not found when backtracking" << endl;
        for (list<pair<int,int>>::iterator iter=p_flow_backtrack->begin(); iter!=p_flow_backtrack->end(); iter++) {
            if(iter->first == flow_id_1 || iter->second == flow_id_1){
                iter->second = flow_id_2;
                found_index = true;
            } else if(iter->first == flow_id_2 || iter->second == flow_id_2){
                iter->second = flow_id_1;
                found_index = true;
            }

            if(found_index) break;
        }
    }

    if(!found_index){
        msg_Out() << METHOD << ": WARNING! color indices " << flow_id_1 << " and " << flow_id_1 << "not found when backtracking either in current partons or reported non-singlets. Will exit." << endl;
        exit(1);
    }
}

void Colour_Handler::InfSoftGGtoG(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2){
    //Infinitely soft gluon exchange in gg->gg (B_1) or gq->gq (A) to avoid creating a colour neutral gluon in gg->g
    //Find first gluon possible to exchange colour with part_in_1
    bool found_parton = false;
    int out_c_1, out_c_2;
    for (list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
        if ((*iter) == part_in_1 || (*iter) == part_in_2) continue;

        if((*iter)->Flav().IsGluon()){
            //Check if gluon for gg->gg
            BackTrackColour(part_in_1->GetFlow(1), (*iter)->GetFlow(2), part_in_1, *iter);
            out_c_1 = (*iter)->GetFlow(1);
            out_c_2 = part_in_1->GetFlow(2);
            part_in_1->SetFlow(1, out_c_1);
            part_in_1->SetFlow(2, -1);
            (*iter)->SetFlow(1, part_in_1->GetFlow(2));
            (*iter)->SetFlow(2, out_c_2);
            found_parton = true;
        } else{
            //Check if q for gq->gq (not allowing gqb->gqb)
            if(!((*iter)->Flav().IsAnti())){
                BackTrackColour((*iter)->GetFlow(1), part_in_1->GetFlow(2), part_in_1, *iter);
                part_in_1->SetFlow(1, part_in_1->GetFlow(1));
                (*iter)->SetFlow(1, -1);
                (*iter)->SetFlow(2, 0);
                part_in_1->SetFlow(2, (*iter)->GetFlow(1));
                found_parton = true;
            }
        }


        if(found_parton) break;
    }

    if(!found_parton){
        msg_Out() << METHOD << ": ERROR: external gluon or quark not found to make an infintely soft exchange to switch colour indices, exiting." << endl;
        exit(1.); 
    }
}

void Colour_Handler::InfSoftQQbtoG(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2){
    //Infinitely soft gluon exchange in qq->qq (A) or gq->gq (A) to avoid creating a colour neutral gluon in qqb->g
    //Find first quark possible to exchange colour with part_in_1 (this parton is chosen to be q and not qb as input)
    bool found_parton = false;
    int out_c_1, out_c_2;
    for (list<shared_ptr<Parton>>::iterator iter=p_partons->begin(); iter!=p_partons->end(); iter++) {
        if ((*iter) == part_in_1 || (*iter) == part_in_2) continue;

        if((*iter)->Flav().IsQuark() && !(*iter)->Flav().IsAnti()){
            out_c_1 = (*iter)->GetFlow(1);
            out_c_2 = part_in_1->GetFlow(1);
            part_in_1->SetFlow(1, out_c_1);
            (*iter)->SetFlow(1, out_c_2);
            found_parton = true;
        } else if((*iter)->Flav().IsGluon()){
            BackTrackColour(part_in_1->GetFlow(1), (*iter)->GetFlow(2), part_in_1, (*iter));
            (*iter)->SetFlow(1, (*iter)->GetFlow(1));
            part_in_1->SetFlow(1, -1);
            part_in_1->SetFlow(2, 0);
            (*iter)->SetFlow(2, part_in_1->GetFlow(1));
            found_parton = true;
        }

        if(found_parton) break;
    }

    if(!found_parton){
        msg_Out() << METHOD << ": ERROR: external gluon or quark not found to make an infintely soft exchange to switch colour indices, exiting." << endl;
        exit(1.); 
    }
}


bool Colour_Handler::UpdateColours_qiqi_qiqi(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                             Flavour flav_out_1, Flavour flav_out_2, double t, double mg2){
    bool updated_colours = false;
    double s = (part_in_1->Momentum()+part_in_2->Momentum()).Abs2();
    double u = -s - t;
    int c_index = 1;
    if(part_in_1->Flav().IsAnti()){
        c_index = 2;
    }
    
    //Colour differential cross section without the common factor 4*pi*alpha_s^2/(9*s^2)
    double sigma_c_A = (s*s + u*u)/((t-mg2)*(t-mg2)) - s*s/(3.*(t-mg2)*(u-mg2));
    double sigma_c_B = (s*s + t*t)/((u-mg2)*(u-mg2)) - s*s/(3.*(t-mg2)*(u-mg2));
    double sigma_tot = sigma_c_A + sigma_c_B;

    int out_c_1, out_c_2;
    if(ran->Get() < sigma_c_A/sigma_tot){
        out_c_1 = part_in_2->GetFlow(c_index);
        out_c_2 = part_in_1->GetFlow(c_index);
        part_in_1->SetFlow(c_index, out_c_1);
        part_in_2->SetFlow(c_index, out_c_2);
        updated_colours = true;
    } else{
        //No change in colour
        updated_colours = true; 
    }

    return updated_colours;
}


bool Colour_Handler::UpdateColours_qiqj_qiqj(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                             Flavour flav_out_1, Flavour flav_out_2){
    bool updated_colours = false;
    
    int c_index_1 = 1;
    if(part_in_1->Flav().IsAnti()){
        c_index_1 = 2;
    }
    int c_index_2 = 1;
    if(part_in_2->Flav().IsAnti()){
        c_index_2 = 2;
    }

    //No interference terms
    if(c_index_1 == c_index_2){ //Both are qiqj->qiqj or qbiqbj->qbiqbj
        int out_c_1, out_c_2;
        out_c_1 = part_in_2->GetFlow(c_index_1);
        out_c_2 = part_in_1->GetFlow(c_index_1);
        part_in_1->SetFlow(c_index_1, out_c_1);
        part_in_2->SetFlow(c_index_1, out_c_2);
        updated_colours = true;
    } else{ //qiqbj->qiqbj, different colour flow from above
        BackTrackColour(part_in_1->GetFlow(c_index_1), part_in_2->GetFlow(c_index_2), part_in_1, part_in_2);
        if(flav_out_1.IsAnti()){ //Find which of the outgoing is q and qb to assign corresponding colour flow
            part_in_1->SetFlow(2, -1);
            part_in_1->SetFlow(1, 0);
            part_in_2->SetFlow(1, part_in_1->GetFlow(2));
            part_in_2->SetFlow(2, 0);
            updated_colours = true;
        } else{
            part_in_2->SetFlow(2, -1);
            part_in_2->SetFlow(1, 0);
            part_in_1->SetFlow(1, part_in_2->GetFlow(2));
            part_in_1->SetFlow(2, 0);
            updated_colours = true;
        }
        
    }
    

    return updated_colours;
}


bool Colour_Handler::UpdateColours_qiqbi_qiqbi(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                               Flavour flav_out_1, Flavour flav_out_2, double t, double mg2){
    bool updated_colours = false;
    double s = (part_in_1->Momentum()+part_in_2->Momentum()).Abs2();
    double u = -s - t;
    
    int c_index_1 = 1;
    if(part_in_1->Flav().IsAnti()){
        c_index_1 = 2;
    }
    int c_index_2 = 1;
    if(part_in_2->Flav().IsAnti()){
        c_index_2 = 2;
    }
    
    //Colour differential cross section without the common factor 4*pi*alpha_s^2/(9*s^2)
    double sigma_c_A = (s*s + u*u)/((t-mg2)*(t-mg2)) - u*u/(3.*(s+mg2)*(t-mg2));
    double sigma_c_B = (t*t + u*u)/((s+mg2)*(s+mg2)) - u*u/(3.*(s+mg2)*(t-mg2));
    double sigma_tot = sigma_c_A + sigma_c_B;

    if(ran->Get() < sigma_c_A/sigma_tot){
        if(flav_out_1.IsAnti()){
            BackTrackColour(part_in_1->GetFlow(c_index_1), part_in_2->GetFlow(c_index_2), part_in_1, part_in_2);
            part_in_1->SetFlow(2, -1);
            part_in_1->SetFlow(1, 0);
            part_in_2->SetFlow(1, part_in_1->GetFlow(2));
            part_in_2->SetFlow(2, 0);
            updated_colours = true;
        } else{
            BackTrackColour(part_in_1->GetFlow(c_index_1), part_in_2->GetFlow(c_index_2), part_in_1, part_in_2);
            part_in_2->SetFlow(2, -1);
            part_in_2->SetFlow(1, 0);
            part_in_1->SetFlow(1, part_in_2->GetFlow(2));
            part_in_1->SetFlow(2, 0);
            updated_colours = true;
        }
    } else{ //Only update if indices switch
        if(flav_out_1.IsAnti() && part_in_2->Flav().IsAnti()){
            part_in_1->SetFlow(2, part_in_2->GetFlow(2));
            part_in_2->SetFlow(1, part_in_1->GetFlow(1));
            part_in_1->SetFlow(1, 0);
            part_in_2->SetFlow(2, 0);
        } else if(flav_out_2.IsAnti() && part_in_1->Flav().IsAnti()){
            part_in_1->SetFlow(1, part_in_2->GetFlow(1));
            part_in_2->SetFlow(2, part_in_1->GetFlow(2));
            part_in_1->SetFlow(2, 0);
            part_in_2->SetFlow(1, 0);
        }
        updated_colours = true; 
    }

    return updated_colours;
}


bool Colour_Handler::UpdateColours_qiqbi_qjqbj(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                               Flavour flav_out_1, Flavour flav_out_2){
    bool updated_colours = false;

    //No interference terms
    if(flav_out_1.IsAnti() && part_in_2->Flav().IsAnti()){
        part_in_1->SetFlow(2, part_in_2->GetFlow(2));
        part_in_2->SetFlow(1, part_in_1->GetFlow(1));
        part_in_1->SetFlow(1, 0);
        part_in_2->SetFlow(2, 0);
    } else if(flav_out_2.IsAnti() && part_in_1->Flav().IsAnti()){
        part_in_1->SetFlow(1, part_in_2->GetFlow(1));
        part_in_2->SetFlow(2, part_in_1->GetFlow(2));
        part_in_1->SetFlow(2, 0);
        part_in_2->SetFlow(1, 0);
    }

    updated_colours = true;

    return updated_colours;
}


bool Colour_Handler::UpdateColours_qiqbi_gg(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                            Flavour flav_out_1, Flavour flav_out_2, double t, double mg2, double mf2){
    bool updated_colours = false;
    double s = (part_in_1->Momentum()+part_in_2->Momentum()).Abs2();
    double u = -s - t;
    
    std::shared_ptr<Parton> in_q, in_qb;
    if(part_in_1->Flav().IsAnti()){
        in_q = part_in_2;
        in_qb = part_in_1;
    } else{
        in_q = part_in_1;
        in_qb = part_in_2;
    }
    
    //Colour differential cross section without the common factor 32*pi*alpha_s^2/(27*s^2)
    double sigma_c_A = u/(t-mf2) - 9.*u*u/(4.*(s+mg2)*(s+mg2));
    double sigma_c_B = t/(u-mf2) - 9.*t*t/(4.*(s+mg2)*(s+mg2));
    double sigma_tot = sigma_c_A + sigma_c_B;

    int out_c_1, out_c_2;
    if(ran->Get() < sigma_c_A/sigma_tot){
        out_c_1 = in_q->GetFlow(1);
        out_c_2 = in_qb->GetFlow(2);
        part_in_1->SetFlow(1, out_c_1);
        part_in_1->SetFlow(2, -1);
        part_in_2->SetFlow(1, part_in_1->GetFlow(2));
        part_in_2->SetFlow(2, out_c_2);
        updated_colours = true;
    } else{
        out_c_1 = in_qb->GetFlow(2);
        out_c_2 = in_q->GetFlow(1);
        part_in_1->SetFlow(1, -1);
        part_in_1->SetFlow(2, out_c_1);
        part_in_2->SetFlow(1, out_c_2);
        part_in_2->SetFlow(2, part_in_1->GetFlow(1));
        updated_colours = true; 
    }

    return updated_colours;
}


bool Colour_Handler::UpdateColours_qig_qig(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                           Flavour flav_out_1, Flavour flav_out_2, double t, double mg2, double mf2){
    bool updated_colours = false;
    double s = (part_in_1->Momentum()+part_in_2->Momentum()).Abs2();
    double u = -s - t;

    std::shared_ptr<Parton> in_g, in_q, out_g, out_q;
    if(part_in_1->Flav().IsGluon()){
        in_g = part_in_1;
        in_q = part_in_2;
    } else{
        in_g = part_in_2;
        in_q = part_in_1;
    }
    if(flav_out_1.IsGluon()){
        out_g = part_in_1;
        out_q = part_in_2;
    } else{
        out_g = part_in_2;
        out_q = part_in_1;
    }

    
    //Colour differential cross section without the common factor pi*alpha_s^2/s^2
    double sigma_c_A = u*u/((t-mg2)*(t-mg2)) - 4*u/(9*(s+mf2));
    double sigma_c_B = s*s/((t-mg2)*(t-mg2)) - 4*s/(9*(u-mf2));
    //Veto A if incoming colour combination would produce a colour neutral gluon in B
    if(in_q->Flav().IsAnti()){
        if(in_q->GetFlow(2) == in_g->GetFlow(1)) sigma_c_B = 0.;
    } else{
        if(in_q->GetFlow(1) == in_g->GetFlow(2)){
            sigma_c_B = 0.;
            /*
            msg_Out() << "sigma_c_B = 0." << endl;
            msg_Out() << "sigma_c_A = " << sigma_c_A << endl;
            msg_Out() << "s = " << s << ", t = " << t << ", u = " << u << endl;
            */
        } 
    }

    double sigma_tot = sigma_c_A + sigma_c_B;

    if(ran->Get() < sigma_c_A/sigma_tot){
        if(in_q->Flav().IsAnti()){
            BackTrackColour(in_q->GetFlow(2), in_g->GetFlow(1), part_in_1, part_in_2);
            out_g->SetFlow(2, in_g->GetFlow(2));
            out_q->SetFlow(2, -1);
            out_q->SetFlow(1, 0);
            out_g->SetFlow(1, out_q->GetFlow(2));
        } else{
            BackTrackColour(in_q->GetFlow(1), in_g->GetFlow(2), part_in_1, part_in_2);
            out_g->SetFlow(1, in_g->GetFlow(1));
            out_q->SetFlow(1, -1);
            out_q->SetFlow(2, 0);
            out_g->SetFlow(2, out_q->GetFlow(1));
        }
        updated_colours = true;
    } else{
        int out_q_c, out_g_c_1, out_g_c_2;
        if(in_q->Flav().IsAnti()){
            out_q_c = in_g->GetFlow(2);
            out_g_c_1 = in_g->GetFlow(1);
            out_g_c_2 = in_q->GetFlow(2);
            out_q->SetFlow(2, out_q_c);
            out_q->SetFlow(1, 0);
        } else{
            out_q_c = in_g->GetFlow(1);
            out_g_c_1 = in_q->GetFlow(1);
            out_g_c_2 = in_g->GetFlow(2);
            out_q->SetFlow(1, out_q_c);
            out_q->SetFlow(2, 0);
        }
        out_g->SetFlow(1, out_g_c_1);
        out_g->SetFlow(2, out_g_c_2);
        updated_colours = true; 
    }

    return updated_colours;
}


bool Colour_Handler::UpdateColours_gg_gg(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                         Flavour flav_out_1, Flavour flav_out_2, double t, double mg2){
    bool updated_colours = false;
    double s = (part_in_1->Momentum()+part_in_2->Momentum()).Abs2();
    double u = -s - t;

    //Colour differential cross section without the common factor 9*pi*alpha_s^2/(4*s^2)
    double sigma_c_A = s*s/((t-mg2)*(t-mg2)) + 2.*s/(t-mg2) + 3. + 2.*t/(s+mg2) + t*t/((s+mg2)*(s+mg2));
    double sigma_c_B = u*u/((s+mg2)*(s+mg2)) + 2.*u/(s+mg2) + 3. + 2.*s/(u-mg2) + s*s/((u-mg2)*(u-mg2));
    double sigma_c_C = t*t/((u-mg2)*(u-mg2)) + 2.*t/(u-mg2) + 3. + 2.*u/(t-mg2) + u*u/((t-mg2)*(t-mg2));

    //Veto A + B if incoming colour combination would produce a colour neutral gluon in C
    if(part_in_1->GetFlow(1) == part_in_2->GetFlow(2) || part_in_1->GetFlow(2) == part_in_2->GetFlow(1)){
        sigma_c_C = 0.;
    }

    double sigma_tot = sigma_c_A + sigma_c_B + sigma_c_C;
    double R = ran->Get();
    
    int out_c_1, out_c_2, out_c_3, out_c_4;
    if(R < sigma_c_A/sigma_tot){
        if(ran->Get() <= 0.5){ //50/50 chance of A_1 or A_2
            BackTrackColour(part_in_1->GetFlow(1), part_in_2->GetFlow(2), part_in_1, part_in_2);
            out_c_1 = part_in_1->GetFlow(2);
            out_c_2 = part_in_2->GetFlow(1);
            part_in_1->SetFlow(1, -1);
            part_in_1->SetFlow(2, out_c_1);
            part_in_2->SetFlow(1, out_c_2);
            part_in_2->SetFlow(2, part_in_1->GetFlow(1));
            
        } else{
            BackTrackColour(part_in_1->GetFlow(2), part_in_2->GetFlow(1), part_in_1, part_in_2);
            out_c_1 = part_in_1->GetFlow(1);
            out_c_2 = part_in_2->GetFlow(2);
            part_in_1->SetFlow(1, out_c_1);
            part_in_1->SetFlow(2, -1);
            part_in_2->SetFlow(1, part_in_1->GetFlow(2));
            part_in_2->SetFlow(2, out_c_2);
        }
        updated_colours = true;
    } else if(R < (sigma_c_A + sigma_c_B)/sigma_tot){
        if(ran->Get() <= 0.5){ //50/50 chance of B_1 or B_2
            BackTrackColour(part_in_1->GetFlow(1), part_in_2->GetFlow(2), part_in_1, part_in_2);
            out_c_1 = part_in_2->GetFlow(1);
            out_c_2 = part_in_1->GetFlow(2);
            part_in_1->SetFlow(1, out_c_1);
            part_in_1->SetFlow(2, -1);
            part_in_2->SetFlow(1, part_in_1->GetFlow(2));
            part_in_2->SetFlow(2, out_c_2);
        } else{
            BackTrackColour(part_in_1->GetFlow(2), part_in_2->GetFlow(1), part_in_1, part_in_2);
            out_c_1 = part_in_2->GetFlow(2);
            out_c_2 = part_in_1->GetFlow(1);
            part_in_1->SetFlow(1, -1);
            part_in_1->SetFlow(2, out_c_1);
            part_in_2->SetFlow(1, out_c_2);
            part_in_2->SetFlow(2, part_in_1->GetFlow(1)); 
        }
        updated_colours = true;
    } else{
        if(ran->Get() <= 0.5){ //50/50 chance of C_1 or C_2
            out_c_1 = part_in_2->GetFlow(1);
            out_c_2 = part_in_1->GetFlow(2);
            out_c_3 = part_in_1->GetFlow(1);
            out_c_4 = part_in_2->GetFlow(2);
        } else{
            out_c_1 = part_in_1->GetFlow(1);
            out_c_2 = part_in_2->GetFlow(2);
            out_c_3 = part_in_2->GetFlow(1);
            out_c_4 = part_in_1->GetFlow(2);
        }
        part_in_1->SetFlow(1, out_c_1);
        part_in_1->SetFlow(2, out_c_2);
        part_in_2->SetFlow(1, out_c_3);
        part_in_2->SetFlow(2, out_c_4);
        updated_colours = true; 
    }

    return updated_colours;
}


bool Colour_Handler::UpdateColours_gg_qiqbi(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                            Flavour flav_out_1, Flavour flav_out_2, double t, double mg2, double mf2){
    bool updated_colours = false;
    double s = (part_in_1->Momentum()+part_in_2->Momentum()).Abs2();
    double u = -s - t;
    
    std::shared_ptr<Parton> out_q, out_qb;
    if(flav_out_1.IsAnti()){
        out_q = part_in_2;
        out_qb = part_in_1;
    } else{
        out_q = part_in_1;
        out_qb = part_in_2;
    }
    
    //Colour differential cross section without the common factor pi*alpha_s^2/(6*s^2)
    double sigma_c_A = u/(t-mf2) - 9.*u*u/(4.*(s+mg2)*(s+mg2));
    double sigma_c_B = t/(u-mf2) - 9.*t*t/(4.*(s+mg2)*(s+mg2));
    double sigma_tot = sigma_c_A + sigma_c_B;

    int out_q_c, out_qb_c;
    if(ran->Get() < sigma_c_A/sigma_tot){
        BackTrackColour(part_in_1->GetFlow(2), part_in_2->GetFlow(1), part_in_1, part_in_2);
        out_q_c = part_in_1->GetFlow(1);
        out_qb_c = part_in_2->GetFlow(2);
        
    } else{
        BackTrackColour(part_in_1->GetFlow(1), part_in_2->GetFlow(2), part_in_1, part_in_2);
        out_q_c = part_in_2->GetFlow(1);
        out_qb_c = part_in_1->GetFlow(2);
    }

    out_q->SetFlow(1, out_q_c);
    out_q->SetFlow(2, 0);
    out_qb->SetFlow(2, out_qb_c);
    out_qb->SetFlow(1, 0);
    updated_colours = true;

    return updated_colours;
}