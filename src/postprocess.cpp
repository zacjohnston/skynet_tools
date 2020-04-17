/// \file r-process.cpp
/// \author jlippuner
/// \since Sep 3, 2014
/// 
/// \brief
///
///

#include <math.h>
#include <string> 
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <numeric>
 
#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/NSE.hpp"
#include "DensityProfiles/ExpTMinus3.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Reactions/ReactionPostProcess.hpp" 
#include "PostProcessing.hpp"
#include "EquationsOfState/NeutrinoHistory.hpp"
#include "Network/NetworkOutput.hpp"
#include "MultiplierReactionLibrary.hpp"

int main(int, char** command) {

    clock_t setup_begin = clock();
    
    std::string runIdent = command[1];
    std::cout << "Analyzing " << runIdent << "..." <<  std::endl; 
    std::string nuclide_library = "../inputs/nuc_actual_z50_n70";

    std::string reactionfile = "../variances_tot_100/json_output0.json";
  
    std::cout << "after nuclide_library set"  <<  std::endl; 
    
    std::vector<std::string> nuclide_library_vec;
    std::cout << "after nuclide_library_vec init"  <<  std::endl; 
    nuclide_library_vec = PostProcessing::ReadNuclidesFromFile(nuclide_library);
    std::cout << "after nuclide_library_vec set"  <<  std::endl; 
   
    auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
        SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml", nuclide_library_vec );
    NetworkOptions opts;
    opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
    opts.MassDeviationThreshold = 1.0E-10;
    opts.IsSelfHeating = false;
    opts.NSEEvolutionMinT9 = 9.5; 
    opts.MaxDt = 0.010;
  
    NetworkOptions opts2;
    opts2.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
    opts2.MassDeviationThreshold = 1.0E-10;
    opts2.IsSelfHeating = true;
    opts2.NSEEvolutionMinT9 = 9.5; 
    
    bool doScreening = false;
  
  /*
    REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
        ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
        "Strong reactions", nuclib, opts, doScreening);
  */
    clock_t setup_multireaclib_begin = clock();
  
    MultiplierReactionLibrary<REACLIBReactionLibrary> strongReactionLibrary(
          REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
          ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
          "Strong reactions", nuclib, opts, false));
  
    strongReactionLibrary = PostProcessing::ReadReactionData(nuclib, strongReactionLibrary, reactionfile);
  
    MultiplierReactionLibrary<REACLIBReactionLibrary> strongReactionLibrary2(
          REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
          ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
          "Strong reactions", nuclib, opts2, false));
  
    strongReactionLibrary2 = PostProcessing::ReadReactionData(nuclib, strongReactionLibrary2, reactionfile);
  
    REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
        ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
        "Weak reactions", nuclib, opts, doScreening);
  
  //  ReactionLibs reactionLibraries { &strongReactionLibrary, &weakReactionLibrary };
    
    clock_t setup_multireaclib_end = clock();
    clock_t setup_readoutput_begin = clock();
  
    std::string filename = runIdent + ".h5";
    std::string filename2 = runIdent + "_2.h5";
    
    auto output = NetworkOutput::ReadFromFile(filename);  
    auto output2 = NetworkOutput::ReadFromFile(filename2);  
  
    std::vector<double> finalYVsA = output.FinalYVsA();
    std::vector<double> tempvstime = output.TemperatureVsTime();
    std::vector<double> DtVsTime = output.DtVsTime();
    std::vector<double> DensVsTime = output.DensityVsTime();
  
    std::vector<double> finalYVsA2 = output2.FinalYVsA();
    std::vector<double> tempvstime2 = output2.TemperatureVsTime();
    std::vector<double> DtVsTime2 = output2.DtVsTime();
    std::vector<double> DensVsTime2 = output2.DensityVsTime();
  
    std::shared_ptr<NeutrinoHistory> nuHist = std::shared_ptr<NeutrinoHistory>(new DummyNeutrinoHistory());
   
  //  auto x = Screening();
    SkyNetScreening screenlib(nuclib);
  //  std::unique_ptr<Screening> screen = screenlib.MakeUniquePtr();
    const Screening* screen = &screenlib;
    clock_t setup_readoutput_makepp_begin = clock();
    ReactionPostProcess postProcess( &strongReactionLibrary, screen, nuclib, output, nuHist );
    ReactionPostProcess postProcess2( &strongReactionLibrary2, screen, nuclib, output2, nuHist );
  
    clock_t setup_readoutput_end = clock();
    auto times = postProcess.GetTimes(); 
    auto times2 = postProcess2.GetTimes();
    unsigned int timesi;
    timesi = times.size();
    clock_t setup_end = clock();
    clock_t for_loops_begin = clock();
  
  /*
  skeleton
  
  for nuclide in nuclib:  
    prod = postProcess.GetNucleusProductionRates(nuclide)
    dest = postProcess.GetNucleusDestructionRates(nuclide)
  
    for flux in prod:
      for i in time:
        prod_total += flux.rate[i]
      prod_vector.append(flux,prod_total)
      
    for flux in dest:
      for i in time:
        dest_total += flux.rate[i]
      dest_vector.append(flux,dest_total)
  *///
  ///
    FILE * flowf = fopen((runIdent + ".flows_for_hendrik").c_str(), "w");
    FILE * f = fopen((runIdent + ".reactions").c_str(), "w");
    FILE * g = fopen((runIdent + ".vary_input").c_str(), "w");
    FILE * h = fopen((runIdent + ".exceptions").c_str(), "w");
    FILE * g2 = fopen((runIdent + ".vary_input2").c_str(), "w");
    int idx = 2;
    typedef std::pair<std::string,std::string> StrPairs;
    typedef std::pair<StrPairs, std::vector<long double>> RatePairs;
  
    std::vector<RatePairs> prod_rates;
    std::vector<RatePairs> dest_rates;
    long double test_value, test_value_diff = 0; 
    std::vector<long double> test_value_prod, test_value_dest; 
    std::vector<std::string> reactions;
    unsigned int count = 0;
    unsigned int ti44_index; 
    std::vector<std::tuple<std::string, std::string, long double, long double>> reaction_tuple;
    std::vector<StrPairs> reaction_strpairs;
    for (unsigned int i = 0; i<nuclide_library_vec.size(); ++i){
    //    if ( nuclide_library_vec[i] == "ti44" )
    //      ti44_index = i; 
        auto prod = postProcess.GetNucleusProductionRates(nuclide_library_vec[i]);
        auto dest = postProcess.GetNucleusDestructionRates(nuclide_library_vec[i]);
        auto prod2 = postProcess2.GetNucleusProductionRates(nuclide_library_vec[i]);
        auto dest2 = postProcess2.GetNucleusDestructionRates(nuclide_library_vec[i]);
        long double prev_total, prev_total2 = 0;
        if (prod.size()!=prod2.size())
            std::cout<<"prod len != prod2 len"<<std::endl; 
        for ( unsigned int j = 0; j<prod.size(); ++j ) {
            long double total_rate = 0;
            std::vector<long double> prod_total;
       //     printf( ("#[%2i] " + prod[j].reac.String() + "\n").c_str(), idx + j);
            long double total_rate2 = 0;
            std::vector<long double> dest_total;
            if(prod[j].reac.String() != prod2[j].reac.String())
                std::cout<<"NOT THE SAME REACTIONS "<<prod[j].reac.String()<<" "<<prod2[j].reac.String()<<std::endl;
            for ( unsigned int t = 0; t<times.size(); ++t ) {
                if( t > 0 )//<times.size()-1)
                {
        //          prod_total.push_back( prod[j].rate[t]);// *(times[t]-times[t-1]));
        //          dest_total.push_back( dest[j].rate[t] );//* (times[t]-times[t-1]) );
                    total_rate+=prod[j].rate[t]*(times[t]-times[t-1]); 
                    total_rate2+=dest[j].rate[t]*(times[t]-times[t-1]); 
                }
            }
            for ( unsigned int t = 0; t<times2.size(); ++t ) {
                if( t > 0 )//<times.size()-1)
                {
        //          prod_total.push_back( prod2[j].rate[t]);// *(times[t]-times[t-1]));
        //          dest_total.push_back( dest2[j].rate[t] );//* (times[t]-times[t-1]) );
                    total_rate+=prod2[j].rate[t]*(times2[t]-times2[t-1]); 
                    total_rate2+=dest2[j].rate[t]*(times2[t]-times2[t-1]); 
//                    if ( t == 1 && j < 100 && i < 5 ){
//                        printf("%10.3E %10.3E %10.3E %10.3E %10.3E", times[t], prod[j].rate[t], dest[j].rate[t], prod2[j].rate[t], dest2[j].rate[t]);
//                        printf("\n");
//                    }
                }
            }
      //      StrPairs str_pairs = std::make_pair( nuclide_library_vec[i] , prod[j].reac.String() );
      //      RatePairs str_rate = std::make_pair( str_pairs , prod_total );
      //      prod_rates.push_back(str_rate);
      //      test_value = std::accumulate(prod_total.begin(), prod_total.end(),0.0);
      //      test_value_prod.push_back(test_value);
      //      std::cout<<"prod:\t"<<nuclide_library_vec[i]<<"\t"<<prod[j].reac.String()<<"\t"<<total_rate<<std::endl;
      
      //      std::cout<<prod[j].reac.String();
            std::string s = prod[j].reac.String();
            std::vector<std::string> result, product;
            std::istringstream iss(s);
            unsigned int equals = 0;
            for(std::string s2 ; iss >> s2; )
            {
                if ( s2 == "+" )
                    continue;
                else if ( s2 == "->" ){
                    equals = result.size();
                    continue;
                }
                else
                    result.push_back(s2);
                    if ( equals > 0 )
                        product.push_back(s2); 
        //        if (s2=="n"||s2=="p"||s2=="he4") 
        //        std::cout<<", "<<result.back();
            }
      
            if ( std::distance(result.begin(), std::find( result.begin(), result.end(), nuclide_library_vec[i] )) >= equals )
            {
//                if (j < 100 && i < 5)
//                    std::cout<<prod[j].reac.String()<< nuclide_library_vec[i]<<std::distance(result.begin(), std::find( result.begin(), result.end(), nuclide_library_vec[i] ))<<std::endl;
//                total_rate2 = -total_rate2;
            }
            else
                total_rate = total_rate;
      
            long double net_rate = total_rate - total_rate2;
      
            int temp_id, temp_z, temp_a, temp_n;
            int temp_prod_z, temp_prod_n;
      
            std::string preapn, postapn, apn = "";
      
            StrPairs str_pairs = std::make_pair(nuclide_library_vec[i], prod[j].reac.String());
      
            bool casecheck = false;
            
            if ( std::find(reactions.begin(), reactions.end(), prod[j].reac.String()) == reactions.end() )
            { 
      //        std::cout<<casecheck<<std::endl;
                casecheck = true;
                reactions.push_back(prod[j].reac.String());
                reaction_strpairs.push_back(str_pairs);
      //        auto tuple_entry = std::make_tuple(nuclide_library_vec[i],prod[j].reac.String(),total_rate,total_rate2);
      //        reaction_tuple.push_back(tuple_entry);
            }  
            else if(std::find(reaction_strpairs.begin(), reaction_strpairs.end(), str_pairs) != reaction_strpairs.end())
            {
                casecheck = true;
                count++;
                total_rate += prev_total;
                total_rate2 += prev_total2;
            }
      //      else
      //        casecheck = false;
            if (casecheck)
            {
                if ( std::find(result.begin(), result.end(), "n") != result.end()  ||std::find(result.begin(), result.end(), "p") != result.end() ||std::find(result.begin(), result.end(), "he4") != result.end())
                {
                    for(unsigned int k = 0; k<result.size(); ++k)
                    {
                        if ( ( result[k] == "n" || result[k] == "p" || result[k] == "he4" ) )
                            if(result[k] == "he4")
                                apn+="a";
                            else
                                apn+=result[k];
                        else if (( std::find(product.begin(), product.end(), "n") == product.end() && std::find(product.begin(), product.end(), "p") == product.end() && std::find(product.begin(), product.end(), "he4") == product.end()) && k>=equals) {
                            apn+="g";
                            temp_id = nuclib.NuclideIdsVsNames().at(result[k]);
                            temp_prod_z = nuclib.Zs()[temp_id];
                            temp_prod_n = nuclib.Ns()[temp_id];
            //              continue;
                        }
                        else if (k < equals) 
                        {
                            temp_id = nuclib.NuclideIdsVsNames().at(result[k]);
                            temp_z = nuclib.Zs()[temp_id];
                            temp_a = nuclib.As()[temp_id];
                            temp_n = nuclib.Ns()[temp_id];
                        }
                        else if (k >= equals) 
                        {
                            temp_id = nuclib.NuclideIdsVsNames().at(result[k]);
                            temp_prod_z = nuclib.Zs()[temp_id];
                            temp_prod_n = nuclib.Ns()[temp_id];
                        }
                        
                    }
          //          std::cout<<"print to g"<<std::endl;
                    fprintf(g, "%s,%s,%10.15LE, %s,%s,%10.15LE",nuclide_library_vec[i].c_str(),prod[j].reac.String().c_str(),total_rate,nuclide_library_vec[i].c_str(),dest[j].reac.String().c_str(),total_rate2);
                    fprintf(g, ", %d,  %d,  %s \n",temp_z,temp_a,apn.c_str());
                }
          //      std::cout<<std::endl;
          
          
        //        if ( std::find(result.begin(), result.end(), "n") != result.end() )// ||std::find(result.begin(), result.end(), "p") != result.end() ||std::find(result.begin(), result.end(), "he4") != result.end())
        //        if(total_rate>1e0 || total_rate2>1e0)
          //      if(result.size()>4)
        //        std::cout<<"print to f"<<std::endl;
                fprintf(f, "%s,%s,%10.15LE, %s,%s,%10.15LE",nuclide_library_vec[i].c_str(),prod[j].reac.String().c_str(),total_rate,nuclide_library_vec[i].c_str(),dest[j].reac.String().c_str(),total_rate2);
                fprintf(f, ",%d,%d,%s\n",temp_z,temp_a,apn.c_str());
                fprintf(flowf, "%d,%d,%d,%d,%10.15LE\n", temp_z, temp_n, temp_prod_z, temp_prod_n, net_rate);
        //          fprintf(f, "%s,%10.15LE, %s,%10.15LE\n",prod[j].reac.String().c_str(),total_rate,dest[j].reac.String().c_str(),total_rate2);
          //      if ( fabs(test_value) > 2*pow(10,9) ) 
          //      if ( prod[j].reac.String().find("ti44") != std::string::npos && fabs(test_value) > 0) 
          //        std::cout<<"prod\t"<<nuclide_library_vec[i]<<"\t"<<str_rate.first<<"\t"<<(test_value)<<"\t"<<j<<std::endl;
              //    std::cout<<str_rate.first<<"\t"<<str_rate.second<<std::endl;
                if(result.size()>4 || apn.length() > 2)
                {
                    fprintf(h, "%s,%s,%10.15LE, %s,%s,%10.15LE",nuclide_library_vec[i].c_str(),prod[j].reac.String().c_str(),total_rate,nuclide_library_vec[i].c_str(),dest[j].reac.String().c_str(),total_rate2);
                    fprintf(h, ", %d,  %d,  %s \n",temp_z,temp_a,apn.c_str());
                }
            }
            
            prev_total = total_rate;
            prev_total2 = total_rate2; 
        }
    }
    return 0;
}
