#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>

#include "PostProcessing.hpp"

int main(int, char** command){

  std::string infile = command[1];
  std::ifstream ifs(infile.c_str(), std::ifstream::in);

  std::string line;
  std::getline(ifs, line);

  std::stringstream ss(line);
  std::string item;
  std::vector<std::string> init_vector, check_vector;
  std::vector<double> prod_vector, dest_vector;

  char delim = ',';
  while (std::getline(ss, item, delim))
    init_vector.push_back(item);

  FILE * f = fopen((infile + ".reduced").c_str(), "w");
  FILE * g = fopen((infile + ".input").c_str(), "w");
  FILE * h = fopen((infile + ".reduced1e0").c_str(), "w");

  bool check = false;
  unsigned int count = 0;
  double flow_thresh = 1e-10;

  while (std::getline(ifs, line)){
    std::stringstream ss2(line);
    std::vector<std::string> check_vector;
    while (std::getline(ss2, item, delim))
      check_vector.push_back(item);
    bool temp_check = false;
    if (init_vector[0] == check_vector[0] && init_vector[1] == check_vector[1])
    {
//      fprintf(f, "%s, %s\n", check_vector[1].c_str(), check_vector[2].c_str());
      temp_check = true; 
    }
//    else if (check)
//      continue; 
    else
    {
      fprintf(f, "%s, %s, %s\n", init_vector[1].c_str(), init_vector[2].c_str(), init_vector[5].c_str());
      temp_check = false;
      double prod = PostProcessing::sciToDub(init_vector[2]);
      double dest = PostProcessing::sciToDub(init_vector[5]);
      if ( prod > flow_thresh || dest > flow_thresh ){
//        std::string s = init_vector[1];
//        std::vector<std::string> result, product;
//        std::istringstream iss(s);
//        unsigned int equals = 0;
//        for(std::string s2 ; iss >> s2; )
//        {
//          if ( s2 == "+" )
//            continue;
//          else if ( s2 == "->" ){
//            equals = result.size();
//            continue;
//          }
//          else
//            result.push_back(s2);
//            if ( equals > 0 )
//              product.push_back(s2); 
////          if (s2=="n"||s2=="p"||s2=="he4") 
////          std::cout<<", "<<result.back();
//        }
//
//      if ( std::distance(result.begin(), std::find( result.begin(), result.end(), smalllibvec[i] )) >= equals )
//        total_rate2 = -total_rate2;
//      else
//        total_rate = total_rate;
//
//      long double net_rate = total_rate + total_rate2;
//
//      int temp_id, temp_z, temp_a, temp_n;
//      int temp_prod_z, temp_prod_n;
        fprintf(h, "%s, %s, %s, %s, %s, %s, %s, %s, %s\n", init_vector[0].c_str(), init_vector[1].c_str(), init_vector[2].c_str(), init_vector[3].c_str(), init_vector[4].c_str(), init_vector[5].c_str(),init_vector[6].c_str(), init_vector[7].c_str(), init_vector[8].c_str());
        fprintf(g, "%s %s %s\n", init_vector[6].c_str(), init_vector[7].c_str(),init_vector[8].c_str());
        count++;
      }
//      std::cout<<init_vector[8][0]<<std::endl;
      if(init_vector[8].size()==2){
        if( init_vector[8][0] == 'n' || init_vector[8][1] == 'n'){
          prod_vector.push_back(prod);
          dest_vector.push_back(dest); 
        }
      }
    }
    check = temp_check;
    init_vector = check_vector;
  }
{
  fprintf(f, "%s, %s, %s\n", init_vector[1].c_str(), init_vector[2].c_str(), init_vector[5].c_str());
  double prod = PostProcessing::sciToDub(init_vector[2]);
  double dest = PostProcessing::sciToDub(init_vector[5]);
  if ( prod > flow_thresh || dest > flow_thresh ){
    fprintf(h, "%s, %s, %s, %s, %s, %s, %s, %s, %s\n", init_vector[0].c_str(), init_vector[1].c_str(), init_vector[2].c_str(), init_vector[3].c_str(), init_vector[4].c_str(), init_vector[5].c_str(),init_vector[6].c_str(), init_vector[7].c_str(), init_vector[8].c_str());
//    fprintf(h, "%s, %s, %s\n", init_vector[1].c_str(), init_vector[2].c_str(), init_vector[5].c_str());
    fprintf(g, "%s %s %s\n", init_vector[6].c_str(), init_vector[7].c_str(),init_vector[8].c_str());
    count++;
  }
//    std::cout<<init_vector[8][0]<<std::endl;
  if(init_vector[8].size()==2){
    if( init_vector[8][0] == 'n' || init_vector[8][1] == 'n'){
      prod_vector.push_back(prod);
      dest_vector.push_back(dest); 
    }
  }
}
  std::sort(prod_vector.rbegin(),prod_vector.rend());
  std::sort(dest_vector.rbegin(),dest_vector.rend());
  std::cout<<prod_vector.size()<<std::endl;
  std::cout<<prod_vector[prod_vector.size()/8]<<std::endl;
  std::cout<<count<<std::endl;}
