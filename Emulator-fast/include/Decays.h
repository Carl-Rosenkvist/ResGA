#pragma once
#include "smash/setup_particles_decaymodes.h"
#include "smash/particles.h"
#include "smash/decaymodes.h"
#include "smash/integrate.h"
#include "smash/particledata.h"
#include "smash/action.h"
#include "smash/collidermodus.h"
#include "smash/config.h"
#include "smash/experiment.h"
#include "smash/forwarddeclarations.h"
#include "smash/library.h"
#include "smash/decaymodes.h"
#include "smash/crosssections.h"
#include "smash/kinematics.h"
#include "smash/clebschgordan.h"
#include "smash/pow.h"
#include "smash/isoparticletype.h"
#include <chrono>

using namespace smash;
struct Parameters {
    std::string name;
    double sqrt_s;
    double spec_mass;
    double spec_charge;
    //Decaymodes pdgs
    PdgCode pdgCode1;
    PdgCode pdgCode2;
    bool two_to_one = false;
    Parameters* next;

    void changeEnergy(double new_sqrt_s){
        sqrt_s = new_sqrt_s;
        if(next){
            next->changeEnergy(new_sqrt_s);
        }
    }

    Parameters(
        const std::string& name, 
        const double& sqrt_s, 
        const double& spec_mass, 
        const int& spec_charge, 
        const std::string& pdg1, 
        const std::string& pdg2) : 
        name(name),
        sqrt_s(sqrt_s),
        spec_mass(spec_mass),
        spec_charge(spec_charge)
        
    {
        pdgCode1 = PdgCode(pdg1); 
        pdgCode2 = PdgCode(pdg2);
        next = nullptr;
    }
    
};


double two_to_one(const ParticleList& incoming_particles_,const Parameters& params);

double formation(

    const ParticleList& incoming_particles_,
    const ParticleType& type_resonance,
    const Parameters& params);


double probability_transit_high(const double sqrt_s_,
    const double region_lower, const double region_upper);



double expBR(
    const double sqrt_s, 
    const double spec_mass, 
    const smash::ParticleTypePtr& ptype, 
    const smash::DecayBranch* mode,
    bool nosamp = false);

 
 

double getBR(const ParticleTypePtr& ptype,const Parameters& params);

double nn_to_resonance_matrix_element(double sqrts,
                                                     const ParticleType& type_a,
                                                     const ParticleType& type_b,
                                                     const int twoI);




double find_nn_xsection_from_type(
    const ParticleList& incoming_particles_,
    const ParticleTypePtrList& list_res_1,
    const ParticleTypePtrList& list_res_2, const Parameters& params);  
 



double getXsection(const ParticleList& incoming_particles_, const Parameters& params);




void initialize(std::string particles_path, std::string decaymodes_path); 



 
