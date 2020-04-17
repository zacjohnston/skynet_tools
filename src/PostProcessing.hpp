/// 
/// 
/// \author khermans21 /// \since Jan 29, 2018
///
/// \brief
/// various functions used for setting up the evolutions
///     reads input nuclei to be used for the network
///     reads input abundances and confirms the mass fractions sum to 1
///     reads in reaction info from json files

/// #ifndef POST_PROCESSING_HPP_
/// #define POST_PROCESSING_HPP_

#include <math.h>
#include <memory> 
#include <algorithm> 
#include <iostream> 
#include <string> 
#include <fstream> 
#include <sstream> 
#include <vector>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include "Reactions/Reaction.hpp"
#include "MultiplierReactionLibrary.hpp"

#include "Utilities/FunctionVsTime.hpp" 
#include <iomanip>


namespace PostProcessing {


// function to convert scientific notation to double
double sciToDub(const std::string& str) {

    std::stringstream ss(str);
    double d = 0;
    ss >> d;

    if (ss.fail()) {
        std::string s = "Unable to format ";
        s += str;
        s += " as a number!";
        throw (s);
    }

    return (d);
}


// Read in the nuclei from the given file
std::vector<std::string> ReadNuclidesFromFile(std::string fname) { 
    std::ifstream fin(fname);

    int n = 0;
    std::string line;

    while (getline(fin, line))
        ++n;
    
    std::string name[n];

    int i = 0;
    std::string name_temp;

    std::ifstream fin2(fname);

    while (fin2 >> name_temp)
    {
        name[i] = name_temp;
        i++;

    }
    
    std::vector<std::string> nuclides;
    nuclides.assign (name,name+n);

    return nuclides;
}



// Read in abundances from given file and check the mass fraction = 1.0
std::vector<double> ReadAbundancesFromFile(std::string fname, std::vector<std::string> nuclides) { 
 
    std::ifstream fin(fname);

    int n = 0;
    std::string line;

    while (getline(fin, line))
        ++n;

    std::string name[n];
    double abundances[n];

    long double abund_total = 0;
    int i = 0;
    std::string name_temp;
    std::string abund_temp;
    std::string nuclide_mass;

    std::ifstream fin2(fname);
    int nuclides_size = nuclides.size();
    std::string abund_name[nuclides_size];
    int mass_no[nuclides_size];

    while (fin2 >> name_temp >> abund_temp )
    {
        if(name_temp == "h1")
            name_temp = "p";
        name[i] = name_temp;
        abundances[i] = sciToDub(abund_temp);
        nuclide_mass = name_temp;
        nuclide_mass.erase(std::remove_if(nuclide_mass.begin(), nuclide_mass.end(), (int(*)(int))std::isalpha), nuclide_mass.end());
        if ( nuclide_mass.length() == 0 ){
            if ( name_temp == "t" )
                nuclide_mass = "3";
            else if ( name_temp == "d" )
                nuclide_mass = "2";
            else
                nuclide_mass = "1"; 
        }
        mass_no[i] = std::stoi(nuclide_mass);
        i++;

    }
    std::vector<std::string> abund_nuclides;
    abund_nuclides.assign (name,name+n);
    


    std::vector<double> abundance_vec;
    i = 0;

    for ( int j = 0; j<nuclides_size; j++ ){
        auto iterator = std::find(abund_nuclides.begin(),abund_nuclides.end(),nuclides[j]);
  
        if ( iterator == abund_nuclides.end() ){ 
            abundance_vec.push_back(0);
        }
        else {
            auto index = std::distance(abund_nuclides.begin(),iterator);
            abundance_vec.push_back(abundances[index]); 
            i++; 
        }
    }
 
    for ( int j = 0; j<nuclides_size; j++ ){
        nuclide_mass = nuclides[j];
        nuclide_mass.erase(std::remove_if(nuclide_mass.begin(), nuclide_mass.end(), (int(*)(int))std::isalpha), nuclide_mass.end());
        if ( nuclide_mass.length() == 0 ){
            if ( nuclides[j] == "t" ){
  //          std::cout<<nuclides[j]<<std::endl;
                nuclide_mass = "3";}
            else if ( nuclides[j] == "d" )
                nuclide_mass = "2";
            else
                nuclide_mass = "1"; 
        }
        mass_no[j] = std::stoi(nuclide_mass);
        abund_total = abund_total + abundance_vec[j]*mass_no[j];
    }

    std::cout.precision(100);
    int stop_loop = 0;
    std::cout << "Sum of input compositions: " << abund_total << "\n"; 
    while( abund_total != 1  && stop_loop < 100 ){
        long double new_total = 0;
            for( i = 0; i < nuclides_size; i++ ){
                if ( abundance_vec[i] == 0 ) int dummy = 1;
                else{
  //                std::cout << abundance_vec[i] << std::endl;
                    abundance_vec[i] = abundance_vec[i] /abund_total;//* pow(abund_total,-1);
  //                std::cout << abundance_vec[i] << std::endl;
                    new_total = new_total + (long double)abundance_vec[i]*(long double)(mass_no[i]);
                }
            }
        abund_total = new_total;
        std::cout.precision(100);
        stop_loop++;
    }
    std::cout<< "Sum of final compositions: " << abund_total<<std::endl; 
    std::cout<<"Number of iterations: "<<stop_loop<<std::endl;
    long double new_total = 0;
        for( i = 0; i < nuclides_size; i++ ){
            abundance_vec[i] = abundance_vec[i]/abund_total;
            new_total = new_total + (long double)abundance_vec[i]*(long double)(mass_no[i]);
        }

    return abundance_vec;
}


// Read in reaction data from the given json file
template<class LibraryType>
MultiplierReactionLibrary<LibraryType> ReadReactionData(
        NuclideLibrary nuclib, MultiplierReactionLibrary<LibraryType> reaclib, 
        std::string filename){
  
    std::vector<double> multipliers;
    std::vector<Reaction> reactions;
    boost::property_tree::ptree root;
    boost::property_tree::read_json(filename,root);
    BOOST_FOREACH( boost::property_tree::ptree::value_type const &v, root ){
  
        std::vector<std::string> reactants;
        std::vector<std::string> products;
        std::vector<int> reactantsN;
        std::vector<int> productsN;
        bool isWeak;
        bool isInverse;
        bool isDecay;
        std::string label;
    
        boost::property_tree::ptree root_child = root.get_child(v.first);
        for (boost::property_tree::ptree::value_type &reactant: root_child.get_child("reactants"))
            reactants.push_back(reactant.second.data());
        for (boost::property_tree::ptree::value_type &product: root_child.get_child("products"))
            products.push_back(product.second.data());
        for (boost::property_tree::ptree::value_type &reactantN: root_child.get_child("reactantsN"))
            reactantsN.push_back(std::stoi(reactantN.second.data()));
        for (boost::property_tree::ptree::value_type &productN: root_child.get_child("productsN"))
            productsN.push_back(std::stoi(productN.second.data()));
        isWeak = root_child.get<bool>("isWeak"); 
        isInverse = root_child.get<bool>("isInverse"); 
        isDecay = root_child.get<bool>("isDecay"); 
        label = root_child.get<std::string>("label");
        multipliers.push_back(root_child.get<double>("multiplier"));
        reactions.push_back( Reaction( reactants , products , reactantsN , productsN ,
             isWeak , isInverse , isDecay , label , nuclib ));
    }
  
    reaclib.SetRateMultipliers(reactions,multipliers);
    return reaclib;

}
}
/// #endif // POST_PROCESSING_HPP_
