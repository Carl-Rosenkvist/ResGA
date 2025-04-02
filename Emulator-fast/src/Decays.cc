#include "Decays.h"

using namespace smash;


double probability_transit_high(const double sqrt_s_,
    const double region_lower, const double region_upper)  {
  if (sqrt_s_ < region_lower) {
    return 0.;
  }

  if (sqrt_s_ > region_upper) {
    return 1.;
  }

  double x = (sqrt_s_ - 0.5 * (region_lower + region_upper)) /
             (region_upper - region_lower);
  assert(x >= -0.5 && x <= 0.5);
  double prob = 0.5 * (std::sin(M_PI * x) + 1.0);
  assert(prob >= 0. && prob <= 1.);

  return prob;
}

double expBR(
    const double sqrt_s, 
    const double spec_mass, 
    const smash::ParticleTypePtr& ptype, 
    const smash::DecayBranch* mode,
bool noSamp)
{
    using namespace smash;


    double m_min = 0;
    double m_max = sqrt_s-spec_mass;
    if(noSamp) return ptype->partial_width(sqrt_s,mode)/ptype->total_width(sqrt_s); 

    static Integrator integrate;
    double range_prob =  integrate(m_min,m_max, [&](double m) {return ptype->spectral_function(m);});
    double norm = 1.0/range_prob;
   
       double exp = norm*integrate(m_min, m_max, [&](double m){
        return //(1.0/mass_integral)*pCM(sqrt_s,0.938,m)*
        (ptype->partial_width(m,mode)/ptype->total_width(m))*ptype->spectral_function(m);
    });
  
//       std::cout<<ptype->name() <<" , "<<exp<<std::endl;

    return exp>really_small ? exp : 0.0;


}

double getBR(const ParticleTypePtr& ptype,const Parameters& params){
double sqrts = params.sqrt_s;
double spec_mass = params.spec_mass;
const PdgCode pdg1 = params.pdgCode1;
const PdgCode pdg2 = params.pdgCode2;
const auto &modes = ptype->decay_modes().decay_mode_list();
    bool has1 = false;
    bool has2 = false;
    for (const auto &decay_branch : modes) {
        for (const ParticleTypePtr decay_into : decay_branch->particle_types()) {
            if (decay_into->pdgcode() == pdg1){ 
                has1=true;
            }
            if (decay_into->pdgcode() == pdg2){ 
                has2=true;
            }
            
            //std::cout<<has1<< " , " << has2<< " , "<<decay_into->name()<<std::endl;
            if(has1 && has2){ 
                double BR = expBR(sqrts,spec_mass,ptype,decay_branch.get(),params.two_to_one);
                if(params.next != nullptr){ 
                    BR *= getBR(decay_into,*params.next);
                
                }
                return BR;

                
            }

        }
    }
    return 0.0;
}



double two_to_one(const ParticleList& incoming_particles_,const Parameters& params){
  CollisionBranchList resonance_process_list;
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  const double m1 = incoming_particles_[0].type().mass();
  const double m2 = incoming_particles_[1].type().mass();

  ParticleTypePtrList possible_resonances =
      list_possible_resonances(&type_particle_a, &type_particle_b);

  // Find all the possible resonances
    double result = 0;
for (const ParticleTypePtr type_resonance : possible_resonances) {
    result += formation(incoming_particles_,*type_resonance,params);

    }
  return result;}


double formation(

    const ParticleList& incoming_particles_,
    const ParticleType& type_resonance,
    const Parameters& params){

    double sqrt_s_ = params.sqrt_s;
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  const double m1 = incoming_particles_[0].type().mass();
  const double m2 = incoming_particles_[1].type().mass();
  const double cm_momentum_sqr = pCM_sqr(sqrt_s_, m1, m2);

  // Calculate partial in-width.
  const double partial_width = type_resonance.get_partial_in_width(
      sqrt_s_, incoming_particles_[0], incoming_particles_[1]);
  if (partial_width <= 0.) {
    return 0.;
  }

  assert(type_resonance.charge() ==
         type_particle_a.charge() + type_particle_b.charge());
  assert(type_resonance.baryon_number() ==
         type_particle_a.baryon_number() + type_particle_b.baryon_number());

  // Calculate spin factor
  const double spinfactor =
      static_cast<double>(type_resonance.spin() + 1) /
      ((type_particle_a.spin() + 1) * (type_particle_b.spin() + 1));
  const int sym_factor =
      (type_particle_a.pdgcode() == type_particle_b.pdgcode()) ? 2 : 1;
  /** Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude.
   * See Eq. (176) in \iref{Buss:2011mx}. */
    double string_prob = probability_transit_high(sqrt_s_,1.9,2.2);
    double result = (1.0-string_prob)*getBR(&type_resonance,params)*spinfactor * sym_factor * 2. * M_PI * M_PI / cm_momentum_sqr *
         type_resonance.spectral_function(sqrt_s_) * partial_width * hbarc *
         hbarc / fm2_mb;
  return result;
}





double nn_to_resonance_matrix_element(double sqrts,
                                                     const ParticleType& type_a,
                                                     const ParticleType& type_b,
                                                     const int twoI) {
  const double m_a = type_a.mass();
  const double m_b = type_b.mass();
  const double msqr = 2. * (m_a * m_a + m_b * m_b);
  


    /* If the c.m. energy is larger than the sum of the pole masses of the
   * outgoing particles plus three times of the sum of the widths plus 3 GeV,
   * the collision will be neglected.
   *
   * This can be problematic for some final-state cross sections, but at
   * energies that high strings are used anyway.
   */
  const double w_a = type_a.width_at_pole();
  const double w_b = type_b.width_at_pole();
  const double uplmt = m_a + m_b + 3.0 * (w_a + w_b) + 3.0;
  if (sqrts > uplmt) {
    return 0.;
  }
  /// NN → NΔ: fit sqrt(s)-dependence to OBE model [\iref{Dmitriev:1986st}]
  if (((type_a.is_Delta() && type_b.is_nucleon()) ||
       (type_b.is_Delta() && type_a.is_nucleon())) &&
      (type_a.antiparticle_sign() == type_b.antiparticle_sign())) {
    return 68. / std::pow(sqrts - 1.104, 1.951);
    /** All other processes use a constant matrix element,
     *  similar to \iref{Bass:1998ca}, equ. (3.35). */
  } else if (((type_a.is_Nstar() && type_b.is_nucleon()) ||
              (type_b.is_Nstar() && type_a.is_nucleon())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → NN*
    if (twoI == 2) {
      return 4.5 / msqr;
    } else if (twoI == 0) {
      const double parametrization = 14. / msqr;
      /** pn → pnη cross section is known to be larger than the corresponding
       * pp → ppη cross section by a factor of 6.5 [\iref{Calen:1998vh}].
       * Since the eta is mainly produced by an intermediate N*(1535) we
       * introduce an explicit isospin asymmetry for the production of N*(1535)
       * produced in pn vs. pp similar to [\iref{Teis:1996kx}], eq. 29. */
      if (type_a.is_Nstar1535() || type_b.is_Nstar1535()) {
        return 6.5 * parametrization;
      } else {
        return parametrization;
      }
    }
  } else if (((type_a.is_Deltastar() && type_b.is_nucleon()) ||
              (type_b.is_Deltastar() && type_a.is_nucleon())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → NΔ*
    return 15. / msqr;
  } else if ((type_a.is_Delta() && type_b.is_Delta()) &&
             (type_a.antiparticle_sign() == type_b.antiparticle_sign())) {
    // NN → ΔΔ
    if (twoI == 2) {
      return 45. / msqr;
    } else if (twoI == 0) {
      return 120. / msqr;
    }
  } else if (((type_a.is_Nstar() && type_b.is_Delta()) ||
              (type_b.is_Nstar() && type_a.is_Delta())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → ΔN*
    return 7. / msqr;
  } else if (((type_a.is_Deltastar() && type_b.is_Delta()) ||
              (type_b.is_Deltastar() && type_a.is_Delta())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → ΔΔ*
    if (twoI == 2) {
      return 15. / msqr;
    } else if (twoI == 0) {
      return 25. / msqr;
    }
  } else if ((type_a.is_deuteron() && type_b.pdgcode().is_pion()) ||
             (type_b.is_deuteron() && type_a.pdgcode().is_pion())) {
    /* This parametrization is the result of fitting d+pi->NN cross-section.
     * Already Breit-Wigner-like part provides a good fit, exponential fixes
     * behaviour around the treshold. The d+pi experimental cross-section
     * was taken from Fig. 2 of [\iref{Tanabe:1987vg}]. */
    return 0.055 / (pow_int(sqrts - 2.145, 2) + pow_int(0.065, 2)) *
           (1.0 - std::exp(-(sqrts - 2.0) * 20.0));
  }

  // all cases not listed: zero!
  return 0.;
}


double find_nn_xsection_from_type(
    const ParticleList& incoming_particles_,
    const ParticleTypePtrList& list_res_1,
    const ParticleTypePtrList& list_res_2, const Parameters& params)  {
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  const double sqrts  = params.sqrt_s; 
  const double s = sqrts*sqrts;
  const double spec_charge = params.spec_charge;
    double mass_sum = incoming_particles_[0].type().mass()+ incoming_particles_[1].type().mass();
        
      double region_lower = 3.5;
      double region_upper = 4.5;
    double string_prob = probability_transit_high(sqrts,region_lower,region_upper);


    double result = 0.0;
static Integrator integrate;
//std::cout<<type_particle_a.name() <<" , "<< type_particle_b.name()<<std::endl;
  // Loop over specified first resonance list
  for (ParticleTypePtr type_res_1 : list_res_1) {
    // Loop over specified second resonance list
    for (ParticleTypePtr type_res_2 : list_res_2) {
      // Check for charge conservation.
        if (type_res_2->charge() != spec_charge) continue;
      if (type_res_1->charge() + type_res_2->charge() !=
          type_particle_a.charge() + type_particle_b.charge()) {
        continue;
        }
    double BR = getBR(type_res_1,params);
      if (params.spec_mass == 0.0) continue;

      // loop over total isospin
      for (const int twoI : I_tot_range(type_particle_a, type_particle_b)) {
        const double isospin_factor = isospin_clebsch_gordan_sqr_2to2(
            type_particle_a, type_particle_b, *type_res_1, *type_res_2, twoI);
        // If Clebsch-Gordan coefficient is zero, don't bother with the rest.
        if (std::abs(isospin_factor) < really_small) {
          continue;
        }

        // Integration limits.
        const double lower_limit = type_res_1->min_mass_kinematic();
        const double upper_limit = sqrts - type_res_2->mass();
        /* Check the available energy (requiring it to be a little above the
         * threshold, because the integration will not work if it's too close).
         */
        if (upper_limit - lower_limit < 1E-3) {
          continue;
        }

        // Calculate matrix element.
        const double matrix_element = nn_to_resonance_matrix_element(
            sqrts, *type_res_1, *type_res_2, twoI);
        if (matrix_element <= 0.) {
          continue;
        }

       

        double resonance_integral =  integrate(lower_limit,upper_limit, [&](double m) {return type_res_1->spectral_function(m) *pCM_from_s(sqrts*sqrts,m,type_res_2->mass());});

        
            const double spin_factor =(type_res_1->spin() + 1.0) * (type_res_2->spin() + 1.0);
        
            
            const double xsection = (1.0-string_prob)*BR*isospin_factor * spin_factor * matrix_element *
                                resonance_integral / (s * pCM_from_s(s,0.938,0.938)) ;

        result += xsection; 
//        std::cout<<type_res_1->name()<< " , "<< type_res_2->name()<<" : "<<xsection<<std::endl;
        }
    }
  }
return incoming_particles_[0].xsec_scaling_factor()*incoming_particles_[1].xsec_scaling_factor()*result;
}


double getXsection(const ParticleList& incoming_particles_, const Parameters& params){


  /* Find whether colliding particles are nucleons or anti-nucleons;
   * adjust lists of produced particles. */
  bool both_antinucleons =
      (incoming_particles_[0].type().antiparticle_sign() == -1) &&
      (incoming_particles_[1].type().antiparticle_sign() == -1);
  const ParticleTypePtrList& nuc_or_anti_nuc =
      both_antinucleons ? ParticleType::list_anti_nucleons()
                        : ParticleType::list_nucleons();
  const ParticleTypePtrList& delta_or_anti_delta =
      both_antinucleons ? ParticleType::list_anti_Deltas()
                        : ParticleType::list_Deltas();
     

    if(!params.two_to_one){
        double result =  find_nn_xsection_from_type(
        incoming_particles_,
        ParticleType::list_baryon_resonances(), nuc_or_anti_nuc, params);
        return result;
    }
    else{
        return two_to_one(incoming_particles_,params);
    }
}



void initialize(std::string particles_path, std::string decaymodes_path) {
  const auto pd = load_particles_and_decaymodes(particles_path, decaymodes_path);
  ParticleType::create_type_list(pd.first);
  DecayModes::load_decaymodes(pd.second);
  ParticleType::check_consistency();
}


