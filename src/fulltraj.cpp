/// code for evolving the set of tracers
/// 
/// after compiling run with ./fulltraj [tracer #] [reaction_vary_id_#]
/// if running for base case just run with ./fulltraj [tracer #] 0 
/// ../sh/run_jobs.sh can be used to send job scripts in ../run_sh/ to hpcc
/// 
/// khermans
//  30/10/2019

/// code based on 
/// \file r-process.cpp
/// \author jlippuner
/// \since Sep 3, 2014
///
/// \brief
///
///

#include <math.h>
#include <memory> 
#include <algorithm> 
#include <iostream> 
#include <string> 
#include <fstream> 
#include <sstream> 
#include <mpi.h>
#include <ctime>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/NSE.hpp"
#include "DensityProfiles/ExpTMinus3.hpp"
#include "DensityProfiles/PowerLawContinuation.hpp"
#include "DensityProfiles/ConstantFunction.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Reactions/NeutrinoReactionLibrary.hpp"
#include "Reactions/ReactionPostProcess.hpp" 
#include "Utilities/FunctionVsTime.hpp" 
#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

#include "PiecewiseProfile.hpp"
#include "PostProcessing.hpp"
#include "MultiplierReactionLibrary.hpp"

#include "EquationsOfState/NeutrinoHistoryBlackBody.hpp"

int main(int narg, char** args) {
    // Model strings
    std::string model = "traj_s12.0";
    std::string tracer_prefix = "stir_s12.0_tracer_";

    std::string PP_PATH = "../"; 
    int p;
    int my_rank;
  
    // for running parallel
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Finalize();
  
    // read the reaction vary ids to run from file
    std::ifstream reac_list("../inputs/reactions_base.txt"); 
    std::string line;
    std::vector<int> reacs;  
    while (std::getline(reac_list, line)) {
        int reac_temp;
        std::stringstream ss(line); 
        ss >> reac_temp;
        reacs.push_back(reac_temp); 
    }

    // name of output (trajectory, int used for determining NSE evolution)    
    std::string runIdent = "test";
    int runIdent_int = 0;
    if (narg > 1) { 
        runIdent = args[1];
        runIdent_int = std::stoi(runIdent); 
    }
    // run 1st, 2nd, 3rd, or 4th quarter of reaction ids
    std::string run_split = "0"; 
    if (narg > 2){
        run_split = args[2];
    }
  
    int N = reacs.size()/4+1; 
    if (N < 300)
        N = reacs.size();
    int offset = 0;
  
    if (run_split == "0"){
        N = 1;
    }else if (run_split == "2"){
        offset = N;
    }else if (run_split == "3"){
        offset = 2*N;
    }else if (run_split == "4")
    {
        offset = 3*N;
        N=reacs.size()-offset;
    }

    int count = N/p;
    int remainder = N % p;
    int start, stop;
    int nid;

    // setup all the inputs   
    std::string id;
    std::string outputPath;
    std::string LagPath = "../inputs/" + model + "/" + tracer_prefix + runIdent + ".dat";
    std::string nuclide_library = "../inputs/nuc_actual_z50_n70";
    std::string abundances_input = LagPath + ".input";
    std::string trajectory = LagPath;
    std::string reactionfile;

    std::cout << "Loading: " << LagPath << std::endl;

    // split reaction-varies amongs threads
    if (my_rank < remainder){
        start = my_rank * (count+1) + offset;
        stop = start + count;
    } else {
        start = my_rank * count + remainder + offset;
        stop = start + (count-1);
    }
    for ( int i = start; i<=stop; ++i){
        if ( i == start )
            std::cout<<start<<"\t"<<stop<<std::endl;

        // get input and output files
        nid = reacs[i];
        id = std::to_string(nid);
        if (narg > 3) {
            id += args[3];
        }
        reactionfile = PP_PATH + "variances_tot_100/json_output"+id+".json";
        outputPath = "../output/" + model + "/" + runIdent;
        std::cout << "Saving to: " << outputPath << std::endl;

        std::vector<std::string> nuclide_library_vec;
        nuclide_library_vec = PostProcessing::ReadNuclidesFromFile(nuclide_library);
        std::vector<double> abundances_vec;
        abundances_vec = PostProcessing::ReadAbundancesFromFile(abundances_input, nuclide_library_vec);

        // setup the nuclei library in skynet
        auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
            SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml", nuclide_library_vec );

        // setup skynet network options for the tracer evolution      
        NetworkOptions opts;
        opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
        opts.MassDeviationThreshold = 1.0E-10;
        opts.IsSelfHeating = false;
        opts.NSEEvolutionMinT9 = 9.5; 
        opts.MaxDt = 0.010;
        opts.DisableStdoutOutput = true; 

        // setup skynet network options for the post-tracer evolution     
        NetworkOptions opts2;
        opts2.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
        opts2.MassDeviationThreshold = 1.0E-10;
        opts2.IsSelfHeating = true;
        opts2.NSEEvolutionMinT9 = 9.5; 
        opts2.DisableStdoutOutput = true; 
      
        SkyNetScreening screen(nuclib);
        HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");
        
        bool doScreening = false;

        // setup skynet reaction libraries for the tracer evolution 
      /*
        REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
            ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
            "Strong reactions", nuclib, opts, doScreening);
      */
        MultiplierReactionLibrary<REACLIBReactionLibrary> strongReactionLibrary(
              REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
              ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
              "Strong reactions", nuclib, opts, false));
        REACLIBReactionLibrary symmetricFission(SkyNetRoot +
            "/data/netsu_panov_symmetric_0neut", ReactionType::Strong, false,
            LeptonMode::TreatAllAsDecayExceptLabelEC,
            "Symmetric neutron induced fission with 0 free neutrons", nuclib, opts, doScreening);
        REACLIBReactionLibrary spontaneousFission(SkyNetRoot +
            "/data/netsu_sfis_Roberts2010rates", ReactionType::Strong, false,
            LeptonMode::TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib,
            opts, doScreening);
      
        // vary the desired reaction rate(s)      
        strongReactionLibrary = PostProcessing::ReadReactionData(nuclib, strongReactionLibrary, reactionfile);
    //    strongReactionLibrary.SetDoOutput(true);
      
      
        // use only REACLIB weak rates
        REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
            ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
            "Weak reactions", nuclib, opts, doScreening);
        NeutrinoReactionLibrary nuLib("../inputs/neutrino_reactions.dat", 
            "weak capture rates", nuclib, opts);
//        nuLib.SetDoOutput(true);

        // output to verify reading the correct neutrino file
        { 
            std::ifstream f("../inputs/neutrino_reactions.dat");

            if (f.is_open())
            std::cout << f.rdbuf();
        }

         
        // Put all of the reaction library pointers in a container for the reaction 
        // network
        ReactionLibs reactionLibraries { &strongReactionLibrary, &symmetricFission,
          &spontaneousFission, &weakReactionLibrary, &nuLib};
      
        // implement the reaction network for tracer evolution 
        ReactionNetwork net(nuclib, reactionLibraries, &helm, &screen, opts);
      
        // Setup constant background and neutrino conditions
        double tstart = 0.0;
        // don't start any runs from NSE
        int nse_runs = -1;
        int t_start_index = 0;
        std::vector<double> times = ReadTimesFromFile( trajectory, 0 );
        std::vector<std::vector<double>> ENU = ReadENUFromFile( trajectory, t_start_index );
        std::vector<std::vector<double>> LNU = ReadLNUFromFile( trajectory, t_start_index );
        std::vector<double> radii = ReadRadiiFromFile( trajectory, t_start_index );
        std::vector<double> sub_times = ReadTimesFromFile( trajectory, t_start_index );
        std::vector<std::vector<double>> etas( sub_times.size(), { 0.0 , 0.0 } ) ;

        // if starting runs from NSE
        if (runIdent_int <= nse_runs) {
            t_start_index = FindStartTFromFile(trajectory);
            if (runIdent_int == 0) {
                t_start_index++;
            }
        } 
        tstart = times[t_start_index];

        auto nuHist = NeutrinoHistoryBlackBody::CreateTimeDependent( sub_times,
            radii, { NeutrinoSpecies::NuE, NeutrinoSpecies::AntiNuE },
            ENU , etas , LNU , false);

        // only load the nuHist files for the tracers with neutrino data
        if (runIdent_int < 98)
            net.LoadNeutrinoHistory(nuHist.MakeSharedPtr());

        // setup reaction libraries for the post-tracer evolution

      /*
        REACLIBReactionLibrary strongReactionLibrary2(SkyNetRoot + "/data/reaclib",
            ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
            "Strong reactions", nuclib, opts2, doScreening);
      */
        MultiplierReactionLibrary<REACLIBReactionLibrary> strongReactionLibrary2(
              REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
              ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
              "Strong reactions", nuclib, opts2, false));
        REACLIBReactionLibrary symmetricFission2(SkyNetRoot +
            "/data/netsu_panov_symmetric_0neut", ReactionType::Strong, false,
            LeptonMode::TreatAllAsDecayExceptLabelEC,
            "Symmetric neutron induced fission with 0 free neutrons", nuclib, opts2, doScreening);
        REACLIBReactionLibrary spontaneousFission2(SkyNetRoot +
            "/data/netsu_sfis_Roberts2010rates", ReactionType::Strong, false,
            LeptonMode::TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib,
            opts2, doScreening);
        REACLIBReactionLibrary weakReactionLibrary2(SkyNetRoot + "/data/reaclib",
            ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
            "Weak reactions", nuclib, opts2, doScreening);
        NeutrinoReactionLibrary nuLib2("../inputs/neutrino_reactions.dat", 
            "weak capture rates", nuclib, opts2);
      
        // vary the reaction rate(s) 
        strongReactionLibrary2 = PostProcessing::ReadReactionData(nuclib, strongReactionLibrary2, reactionfile);
     //   strongReactionLibrary2.SetDoOutput(true);
      
        ReactionLibs reactionLibraries2 { &strongReactionLibrary2, &symmetricFission2,
          &spontaneousFission2, &weakReactionLibrary2, &nuLib2};
      
        // implement the reaction network for post-tracer evolution 
        ReactionNetwork net2(nuclib, reactionLibraries2, &helm, &screen, opts2);
     
        // dummy T0    
        double T0 = 1.0;
        // run to 1e15 s so all nuclei decay to stability
        double tfinal = 1.e15; 

        // read in density, temperature, and ye profiles    
        auto DensProfile = ReadDensityProfileFromFile(trajectory); 
        auto TempProfile = ReadTempProfileFromFile(trajectory); 
        auto YeProfile = ReadYeProfileFromFile(trajectory);
      
        // Read the density profile just into vectors for genlinfun
        std::vector<double> rhos = ReadDensityFromFile( trajectory );
        double tinter = times.back();
    
        NSE nse(net.GetNuclideLibrary(), &helm, &screen);
    
        // density profile for use with the power law continuation
        auto genlinfun = GeneralPiecewiseLinearFunction<double>(times, rhos,0); 
    
        // density profile for evolving after the tracer ends  
        auto DensProfile2 = PowerLawContinuation(genlinfun,-3,0.1);
      
        
        // run evolution from tracer 
        std::vector<double> TempVsTime, FinalY;
        if (runIdent_int <= nse_runs) {
            //  if evolving from nse, use this 
            auto nseResult = nse.CalcFromTemperatureAndDensity(TempProfile(tstart),DensProfile(tstart),YeProfile(tstart));
            auto output = net.Evolve(nseResult.Y(), tstart, tinter, &TempProfile, &DensProfile, outputPath );   
            TempVsTime = output.TemperatureVsTime();
            FinalY = output.FinalY();
        } else {
            auto output = net.Evolve(abundances_vec, tstart, tinter, &TempProfile, &DensProfile, outputPath );   
            TempVsTime = output.TemperatureVsTime();
            FinalY = output.FinalY();
        }
 
        // get final temp from the end of the tracer evolution
        double final_temp = TempVsTime.back();
    
        // run evolution from end of tracer using power law profile     
        auto output2 = net2.EvolveSelfHeatingWithInitialTemperature(FinalY,tinter,tfinal,final_temp,&DensProfile2,outputPath + "_2");

    }
    return 0;

}

