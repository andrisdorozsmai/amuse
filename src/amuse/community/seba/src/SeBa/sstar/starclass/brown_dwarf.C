//
// brown_dwarf.C
//
//to be based on Dantona, F., Mazzitelli, I., 1985, ApJ 296, 502
// and  2000astro.ph..5557, Chabrier, G.; Baraffe, I.; Allard, F.; 
// Hauschildt, P.

#include "brown_dwarf.h"
#include "main_sequence.h"
#include "proto_star.h"

brown_dwarf::brown_dwarf(proto_star & p) : single_star(p) {
    
      delete &p;

      real m_tot = get_total_mass();
      core_mass = brown_dwarf_core_mass();
      envelope_mass = m_tot - core_mass;

      last_update_age = 0;
      relative_age = 0;
	  
      instantaneous_element();
      update();

      post_constructor();
}

brown_dwarf::brown_dwarf(main_sequence & m) : single_star(m) {
    
      delete &m;

      real m_tot = get_total_mass();
      core_mass = brown_dwarf_core_mass();
      envelope_mass = m_tot - core_mass;

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;
      instantaneous_element();
      update();

      post_constructor();
}

void brown_dwarf::instantaneous_element() {

     luminosity = 1.e-4; 

     core_radius = radius = brown_dwarf_radius();
       
     update();
}

real brown_dwarf::get_evolve_timestep() {
    // (GN+SPZ Apr 28 1999) was a bit too small
    //  return max(next_update_age - relative_age
    //	     -0.5*cnsts.safety(minimum_timestep),
    //	     cnsts.safety(minimum_timestep));
    
    // (GN+SPZ May  5 1999) type change time must be small because of rapid
    // growth of giants at end phase 0.0001 seems to be OK (?)
    // return max(next_update_age - relative_age - (0.5*0.001), 0.001);

    return max(next_update_age - relative_age, 0.0001);
}



void brown_dwarf::evolve_element(const real end_time) {

    real dt = end_time - current_time;
    current_time = end_time;
    relative_age += dt;

    next_update_age = relative_age + cnsts.safety(maximum_timestep);

	//Burrows & Libert 1993, J. Rev. Mod. Phys. 65, 301
	luminosity = 938 * pow(relative_mass, 2.64); 
	if(relative_age>1){
		luminosity = 938 * pow(relative_mass, 2.64) / pow(relative_age, 1.3);
	}

    core_radius = radius = brown_dwarf_radius();
       
    update();
    this->determine_beta(get_mzams(), current_time);
    //if (identity == 0){
    //  	cout<<current_time<<","<<get_mzams()<<","<<beta_interpolated<<","<<k_interpolated<<endl;
    //}
     }

void brown_dwarf::update() {

     detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
     effective_radius = radius;

     }



// (SilT May 26 2021) 
// fit to Chen & Kipping (2017)
real brown_dwarf::brown_dwarf_radius() {

    real m_earth = cnsts.parameters(Mearth);
    real m_jupiter = cnsts.parameters(Mjupiter);
    real r_jupiter = cnsts.parameters(Rjupiter);           
       
    if (get_total_mass() > 0.414 * m_jupiter) // Jovian planets & brown dwarfs
        return (0.027095 + 1.205219* pow(get_total_mass()/m_jupiter,-0.044))*r_jupiter;
    else if (get_total_mass() >2.04 * m_earth) // Neptunian planets
        return  2.151774*r_jupiter * pow(get_total_mass()/m_jupiter, 0.589); 
    else // Terran planets
        return 0.31832 * pow(get_total_mass(), 0.28); 
    
     }


real brown_dwarf::brown_dwarf_core_mass() {
    
        return 0.01 * get_total_mass();
     }

star* brown_dwarf::subtrac_mass_from_donor(const real dt, real& mdot) {

    mdot = mdot_limit(dt, mdot);
    
        if (mdot<=envelope_mass)
	  envelope_mass -= mdot;
        else if (mdot>envelope_mass) 
	  envelope_mass = 0;

        return this;
     }


real brown_dwarf::add_mass_to_accretor(real mdot, bool, const real dt) {

        if (mdot<0) {
           cerr << "brown_dwarf::add_mass_to_accretor(mdot="
                 << mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   mdot = 0;
        }

        mdot = accretion_limit(mdot, dt);
 
        envelope_mass += mdot;
	relative_mass = max(relative_mass, get_total_mass());

	set_spec_type(Accreting);
	
        return mdot;
     }

real brown_dwarf::accretion_limit(const real mdot, const real dt) {

  if (dt < 0) return mdot;

        real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;

        if(mdot>=eddington)
	  return eddington;

        return mdot;
     }


real brown_dwarf::zeta_thermal() {

     return 0;
}

star* brown_dwarf::merge_elements(star* str) {

     real merger_core = str->get_core_mass();

     add_mass_to_accretor(str->get_envelope_mass(), 
			  cnsts.parameters(spiral_in_time), str->hydrogen_envelope_star());

     if (relative_mass < get_total_mass() + merger_core)
       relative_mass=get_total_mass() + merger_core;
     core_mass += merger_core;

     spec_type[Merger]=Merger;
     instantaneous_element();

     return this;
}

star* brown_dwarf::reduce_mass(const real mdot) {

      if (envelope_mass < mdot)
	envelope_mass = 0;
      else
	envelope_mass -= mdot;

      return this;
}

real brown_dwarf::gyration_radius_sq() {

    return cnsts.parameters(convective_star_gyration_radius_sq); 
}


stellar_type brown_dwarf::get_element_type() {

    // (SilT January 18 2021) 
    if (get_total_mass() < cnsts.parameters(brown_dwarf_mass_limit))
//  if (get_total_mass() < 0.1*cnsts.parameters(minimum_main_sequence))
      return Planet;
    else
      return Brown_Dwarf;
}


// Function reading in a csv table with predefined number of entries into a n*m array
void brown_dwarf::read_in_csv_table(real beta_array_SSO[14][15], char* filename){
   char buffer[1024] ;
   char *record,*line;
   int i=0,j=0;
   FILE *fstream = fopen(filename,"r");
   if(fstream == NULL){
      cout<<"\n file opening failed "<<endl;
      exit(-1);
   }
   while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL){
     record = strtok(line,",");
     j = 0;
     while(record != NULL){
     	beta_array_SSO[i][j++] = atof(record) ;
     	record = strtok(NULL,",");
     }
     i++;
   }
}

// Function for determining beta for given mass and time
void brown_dwarf::determine_beta(const real mass, const real time){
	//cout<<"do we ever get here? "<<endl;
	// convert time to yrs from Myrs
	real time_yrs = pow(10,6) * time;
	const int size_mbins = 14;
	const int size_tbins = 15;
	real m_bins [size_mbins] = {5.00e-04, 1.00e-03, 3.00e-03, 5.00e-03, 0.01,
						0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
	
	// check whether we are in relevant parameter space
	real max_time; 
	real max_time_for_a_given_mass [size_mbins];
	real max_mass_difference = 100.0;
	
	for(int i  = 0; i < size_mbins; i++){
		max_time_for_a_given_mass[i] = *max_element(timestamp_array[i], timestamp_array[i] + size_tbins); 
	}
	max_time =  *max_element(max_time_for_a_given_mass, max_time_for_a_given_mass + size_mbins); 
	
	real max_mass = *max_element(m_bins, m_bins + size_mbins);
	real min_mass = *min_element(m_bins, m_bins + size_mbins);
	
	if (mass >= min_mass && mass <= max_mass && time_yrs < max_time){
			
		real delta_temp;	
		real m_upper;
		real m_lower;    
		int im_upper;
		int im_lower;
		real delta_lower = -max_mass_difference ; 
		real delta_upper = max_mass_difference ; 
	
		for(int i=0;i<size_mbins;i++){	
   			delta_temp = m_bins[i] - mass;
   			if (delta_temp > 0) {
   				if (delta_temp < delta_upper) {
   					m_upper = m_bins[i];
   					im_upper = i;
   					delta_upper = delta_temp;
   				}
   			}else if (delta_temp < 0) {		
				if (delta_temp > delta_lower) {
					m_lower = m_bins[i];
					im_lower = i;
					delta_lower = delta_temp;
				}
   			}else{
   				m_lower = m_bins[i];
   				im_lower = i;
   				im_upper = i+1;
   			}
		}
					    
		//cout<<"m_lower: "<<m_lower<<", m actual: "<<mass<<", m_upper: "<<m_upper<<endl;
	
		//Second thing, for each m_lower, m_upper get the t_lower and t_upper
		// m_lower: t_ll, t_lu and m_upper: t_ul and t_uu
		real delta_t_upper = max_time;
		real delta_t_lower = -max_time; // some large number
		real delta_temp_m_lower;
		real delta_temp_m_upper;
	
		real t_ll;
		real t_lu;
		real t_ul;
		real t_uu;
	
		int it_ll;
		int it_lu;
		int it_ul;
		int it_uu;
		
		// are we in the correct time range?
		if (time_yrs > max_time_for_a_given_mass[im_lower] || time_yrs > max_time_for_a_given_mass[im_upper]){
			//do nothing and skip	
		}else{
			for(int i=0;i<size_tbins;i++){
				delta_temp_m_lower = timestamp_array_SSO[im_lower][i] - time_yrs;
				if (delta_temp_m_lower > 0) {
					if (delta_temp_m_lower < delta_t_upper) {
						t_lu = timestamp_array_SSO[im_lower][i];
						it_lu = i;
						delta_t_upper = delta_temp_m_lower;
					}
				}else if (delta_temp_m_lower <= 0) {		
					if (time_yrs == max_time_for_a_given_mass[im_lower]){
						t_lu = timestamp_array_SSO[im_lower][i];
						t_ll = timestamp_array_SSO[im_lower][i-1];
						it_lu = i;
						it_ll = i-1;
					}else if (delta_temp_m_lower > delta_t_lower) {
						t_ll = timestamp_array_SSO[im_lower][i];
						it_ll = i;
						delta_t_lower = delta_temp_m_lower;
					}
				}
			}	
	
			delta_t_upper = max_time;
			delta_t_lower = -max_time;
	
			for(int i=0;i<size_tbins;i++){
				delta_temp_m_upper = timestamp_array_SSO[im_upper][i] - time_yrs;
				if (delta_temp_m_upper > 0) {
					if (delta_temp_m_upper < delta_t_upper) {
						t_uu = timestamp_array_SSO[im_upper][i];
						it_uu = i;
						delta_t_upper = delta_temp_m_upper;
					}
				}else if (delta_temp_m_upper <= 0) {	
				    if (time_yrs == max_time_for_a_given_mass[im_upper]){
						t_lu = timestamp_array[im_upper][i];
						t_ll = timestamp_array[im_upper][i-1];
						it_lu = i;
						it_ll = i-1;
						break;
				 	}else if (delta_temp_m_upper > delta_t_lower){
						t_ul = timestamp_array_SSO[im_upper][i];
						it_ul = i;
						delta_t_lower = delta_temp_m_upper;
				 	}
				 }	
			}  

		//cout<<"m_lower: "<<m_lower<<", t_ll: "<<t_ll<<", t_lu: "<<t_lu<<endl;
		//cout<<"m_upper: "<<m_upper<<", t_ul: "<<t_ul<<", t_uu: "<<t_uu<<endl;
	
		// step three, interpolation in mass and time...
		//1st obtain beta for mlower at t and m_upper at t
		//then interpolate fore m
		real beta_m_lower_at_t = beta_array_SSO[im_lower][it_ll] + (time_yrs - t_ll) * (beta_array_SSO[im_lower][it_lu] - beta_array_SSO[im_lower][it_ll]) / (t_lu - t_ll);
		real beta_m_upper_at_t = beta_array_SSO[im_upper][it_ul] + (time_yrs - t_ul) * (beta_array_SSO[im_upper][it_uu] - beta_array_SSO[im_upper][it_ul]) / (t_uu - t_ul);

		beta_interpolated = beta_m_lower_at_t + (mass - m_lower) * (beta_m_upper_at_t - beta_m_lower_at_t) / (m_upper - m_lower);
	
		real log_k_interpolated;
	
		real k_m_lower_at_t = k_array_SSO[im_lower][it_ll] + (time_yrs - t_ll) * (k_array_SSO[im_lower][it_lu] - k_array_SSO[im_lower][it_ll]) / (t_lu - t_ll);
		real k_m_upper_at_t = k_array_SSO[im_upper][it_ul] + (time_yrs - t_ul) * (k_array_SSO[im_upper][it_uu] - k_array_SSO[im_upper][it_ul]) / (t_uu - t_ul);

		log_k_interpolated = k_m_lower_at_t + (mass - m_lower) * (k_m_upper_at_t - k_m_lower_at_t) / (m_upper - m_lower);
		k_interpolated = pow(10, log_k_interpolated);
	}
		
		//cout<<"k_interpolated bd: "<<k_interpolated<<", beta interpolated bd: "<<beta_interpolated<<endl;
	
	}else{
		//cout<<"Mass or time outside of interval!"<<endl;
	}

}
	
