#include "Decays.h"
#include "MonteCarloEmulator.h"


void runDecays(char* argv[]){

std::map<std::string,std::vector<double>> energies;  


/* energies["K⁺+p+Λ"] = {2.54916997, 2.55016081, 2.55118559, 2.55217598, 2.55320031,
       2.55419026, 2.55518   , 2.55620364, 2.5590004 , 2.56142057,
       2.56240872, 2.5644184 , 2.56870721, 2.57859634, 2.58819203,
       2.60779083, 2.7525205 , 2.83352778, 2.86453897, 2.7525205 ,
       2.78750718, 2.83251068, 2.60222457, 2.68575601, 2.63249552,
       2.66254264, 2.71867848, 2.63584521, 2.66918907, 2.68244825,
       2.70422738, 2.75066823, 2.70455642, 2.97810521, 3.34887618,
       3.47071231, 3.48400608, 3.4972528 , 3.62722724, 3.62722724,
       3.77761719, 3.8555388 , 4.07843806, 4.93436016, 6.84339572};
*/

/*
energies["K⁺+p+Λ"] = {2.97810521, 2.66254264, 2.55118559, 2.55620364, 2.55518   ,
       2.83251068, 2.56240872, 2.70455642, 2.55419026, 2.5644184 ,
       2.5590004 , 3.62722724, 6.84339572, 3.47071231, 3.77761719,
       3.34887618, 3.4972528 , 3.8555388 , 4.07843806, 4.93436016,
       3.48400608, 3.62722724};
*/

energies["K⁺+p+Λ"] = {2.7525205 , 2.55419026, 2.86453897, 2.5590004 , 2.58819203,
       2.55518   , 2.55016081, 2.56142057, 2.7525205 , 2.71867848,
       2.68244825, 2.63584521, 3.8555388 , 4.93436016, 6.84339572,
       3.4972528 , 3.77761719, 3.62722724, 3.47071231, 2.97810521,
       3.48400608, 3.62722724, 4.07843806, 3.34887618};


energies["K⁰+p+Σ⁺"] = {2.75131825, 2.78657377, 2.83168408, 2.97810521, 3.8555388 ,
       4.07843806, 4.53966087};

energies["K⁺+n+Σ⁺"] = {2.63584521, 2.66918907, 2.68244825, 2.70422738, 2.75066823,
       2.97810521, 3.34887618, 3.47071231, 3.62722724, 3.77761719};

energies["K̅⁰+K⁰+p+p"] = {3.34887618, 3.4972528 , 3.62722724, 3.8555388 , 4.07843806,
       4.53966087};

energies["K⁺+Σ⁺"] ={1.69440708, 1.72919162, 1.73241856, 1.7532553 , 1.75697019,
       1.76384936, 1.77700665, 1.78276629, 1.78902917, 1.79111211,
       1.81335795, 1.82105638, 1.82208046, 1.84548243, 1.85154014,
       1.87009813, 1.87657616, 1.8865001 , 1.89094933, 1.89637347,
       1.90913445, 1.92569738, 1.93440993, 1.9392337 , 1.95411339,
       1.95459153, 1.96888302, 1.96983217, 1.98495816, 1.99247894,
       2.01580389, 2.01904827, 2.02506   , 2.03059377, 2.05258214,
       2.05713395, 2.05940614, 2.07434115, 2.08872315, 2.10612034,
       2.11808201, 2.13743358, 2.14660884, 2.15791532, 2.18806892,
       2.19362376, 2.20171744, 2.22414056, 2.24300491, 2.26129939,
       2.28191221, 2.30356224, 2.3177495 , 2.3410721 , 2.35543296,
       2.47074911, 2.63948208, 2.80479395, 2.90012184, 3.20723056,
       3.35026323, 3.74677962, 4.49778812, 4.83941764, 5.56109162};

energies["K⁺+Σ⁻"] = {1.72163965, 1.74153056, 1.76384936, 1.79007093, 1.7978656 ,
       1.82310398, 1.84447094, 1.93054243, 1.97362432, 1.98542903,
       2.02043716, 2.02043716, 2.09319755, 2.13305071, 2.13743358,
       2.15052927, 2.18078379, 2.21908326, 2.22329847, 2.26502342,
       2.30600048, 2.36377014, 2.40930715, 2.44405711, 2.46315216,
       2.50465192, 2.55650099, 2.56016452, 2.60372681, 2.63236971,
       2.86433538, 2.88391017, 2.90012184, 2.95140245};

 energies["K⁰+Λ"] =  {1.61758006, 1.62101595, 1.62615715, 1.62786752, 1.63014543,
       1.63242038, 1.63298866, 1.63922757, 1.64092525, 1.64826298,
       1.65163933, 1.65388664, 1.66061146, 1.66061146, 1.66228867,
       1.66786793, 1.67009472, 1.67231871, 1.67287428, 1.67564949,
       1.67842038, 1.68118696, 1.68284484, 1.68284484, 1.68339713,
       1.68394924, 1.68450119, 1.68670725, 1.6889106 , 1.69111122,
       1.69385819, 1.69385819, 1.69385819, 1.69550436, 1.6982446 ,
       1.70152736, 1.70425842, 1.70753021, 1.71459874, 1.71622601,
       1.72163965, 1.72163965, 1.72326059, 1.72380057, 1.72488007,
       1.73510331, 1.74206514, 1.74259956, 1.74526933, 1.7500651 ,
       1.7580302 , 1.77648216, 1.79163248, 1.7978656 , 1.81951921,
       1.82463823, 1.84447094, 1.84699867, 1.84750381, 1.87906193,
       1.89785016, 1.90864518, 1.93054243, 1.93440993, 1.9377878 ,
       1.96650817, 1.97362432, 1.97835449, 1.98542903, 1.99950447,
       2.02043716, 2.02043716, 2.02736754, 2.05940614, 2.09319755,
       2.10434254, 2.13305071, 2.13743358, 2.15052927, 2.15921615,
       2.18078379, 2.18078379, 2.18292897, 2.20851042, 2.21908326,
       2.22329847, 2.25881333, 2.26502342, 2.30600048, 2.31613246,
       2.36377014, 2.40736199, 2.44405711, 2.46315216, 2.46315216,
       2.50465192, 2.56016452, 2.6019262 , 2.63236971, 2.86269813,
       2.8676071 , 2.88391017, 2.95140245, 3.05752617, 3.4874425 ,
       3.98930685, 4.43477051};

energies["K⁺+p+Σ⁰"] = {2.62451216, 2.62649268, 2.62850579, 2.62921015, 2.63152359,
       2.63199283, 2.63440518, 2.63530944, 2.6374186 , 2.64180079,
       2.65141815, 2.66121201, 2.68059467, 2.74852235, 2.78351436,
       2.8285349 , 2.75066823, 2.66918907, 2.68244825, 2.70422738,
       2.97810521, 3.34887618, 3.47071231, 3.4972528 , 3.62722724,
       3.62722724, 3.77761719, 3.8555388 , 4.07843806};




    using namespace smash;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    double sqrts =0; 
    initialize(argv[2],argv[3]);  
    



    PdgCode pdg_a("2212"), pdg_b("2212");
    const ParticleType &a = ParticleType::find(pdg_a);
    const ParticleType &b = ParticleType::find(pdg_b);
    
    ParticleData a_data(a), b_data(b);
    
    double s = sqrts*sqrts;
    double ma = a.mass();
    double mb = b.mass();
    double m2 = (ma + mb) * (ma+mb);
    double pz = pCM_from_s(s,ma,mb);
    
    a_data.set_4momentum(ma, 0.0, 0.0, 0);
    b_data.set_4momentum(mb, 0.0, 0.0, 0);
     
    ParticleList list = {a_data,b_data};

    double m_p = 0.938;
    double m_n = 0.939;

    Parameters param1("K⁺+p+Λ",sqrts,m_p,1,"3122","321");
    Parameters param2("K⁰+p+Σ⁺",sqrts,m_p,1,"3222","311"); 
    Parameters param3("K⁺+n+Σ⁺",sqrts,m_n,0,"3222","321");

    
    Parameters param4 ("K̅⁰+K⁰+p+p",sqrts,m_p,1,"333","2212");       
    Parameters param5 ("K̅⁰+K⁰+p+p",sqrts,m_p,1,"9010221","2212");        
    Parameters param6 ("K̅⁰+K⁰+p+p",sqrts,m_p,1,"9000111","2212");       
    
    
    Parameters param7("K⁺+Σ⁺",sqrts,0,0,"3222","321");
    param7.two_to_one = true;
   
    Parameters param8("K⁺+Σ⁻",sqrts,0,0,"3112","321");
    param8.two_to_one = true;
 
    Parameters param9("K⁰+Λ",sqrts,0,0,"3122","311");
    param9.two_to_one = true;
   
    Parameters param10("K⁺+p+Σ⁰",sqrts,m_p,1,"3212","321");




    Parameters paramsKK ("next",sqrts,m_p,1,"311","-311");    
    param4.next = &paramsKK;
    param5.next = &paramsKK;
    param6.next = &paramsKK;
    


    PdgCode pdg_pion("211");
    const ParticleType &pion= ParticleType::find(pdg_pion);
    ParticleData pion_data(pion);

    pion_data.set_4momentum(pion.mass(),0,0,0);
    ParticleList list_pion = {a_data,pion_data};


    PdgCode pdg_pion_neg("-211");
    const ParticleType &pion_neg= ParticleType::find(pdg_pion_neg);
    ParticleData pion_data_neg(pion_neg);
    std::cout<<pion_data_neg.type().name()<<std::endl;
    pion_data_neg.set_4momentum(pion_neg.mass(),0,0,0);
    ParticleList list_pion_neg = {a_data,pion_data_neg};

    //std::cout<<param1.name + " : " <<getXsection(list,param1)<<std::endl;  
    
    //std::cout<<param2.name + " : " <<getXsection(list,param2)<<std::endl;
    //std::cout<<param3.name + " : " <<getXsection(list,param3)<<std::endl;     
    //std::cout<<param4.name + " : " <<getXsection(list,param4)+getXsection(list,param5)+getXsection(list,param6)<<std::endl; 
    // Initialize the results map with zeros for each energy value
    


    std::vector<Parameters> parameters = {param1, param2, param3, param4, param5, param6,param7,param8,param9,param10};

    //std::vector<Parameters> parameters = {param8};
    std::map<std::string, std::vector<double>> results;
    std::string save_path = argv[1];
    for (const auto& energy_pair : energies) {
        const std::string& name = energy_pair.first;
        size_t size = energy_pair.second.size();
        results[name] = std::vector<double>(size, 0.0);
    }

    // Loop through each parameter set and compute the results
    for (auto& param : parameters) {
        const std::string& name = param.name;
        const std::vector<double>& energy_values = energies[name];
        std::vector<double>& result_vector = results[name];

        for (size_t i = 0; i < energy_values.size(); ++i) {
            param.changeEnergy(energy_values[i]);
            double res;

            if(!param.two_to_one ) res = getXsection(list, param); // Replace with your actual list
            else{
                res = getXsection((param.name == "K⁺+Σ⁻" or param.name == "K⁰+Λ") ? list_pion_neg : list_pion,param);
            }
            result_vector[i] += res;
        }
    }

    for(const auto& result_pair : results){
        
        std::ofstream outFile(save_path+"/"+result_pair.first + ".csv");
        //for(const auto& sigma : result_pair.second) outFile<<sigma<<"\n";
        auto energy_v = energies[result_pair.first];
        auto sigmas = result_pair.second;
        outFile << "sigma" <<"," <<"ssqrt"<<std::endl;
        for(int i = 0; i < energy_v.size(); ++i ){
            outFile << sigmas[i] << "," << energy_v[i]<<std::endl;
        }
    
        outFile.close();
    }

}


int main(int args, char* argv[]) {



    runDecays(argv);
    return 0;
}
