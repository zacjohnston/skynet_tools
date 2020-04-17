/// \file MultiplierReactionLibrary.hpp
/// \author lroberts
/// \since May 4, 2018
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_MULTIPLIERREACTIONLIBRARY_HPP_
#define SKYNET_REACTIONS_MULTIPLIERREACTIONLIBRARY_HPP_

#include <iostream> 
#include <stdexcept>

#include "Reactions/ReactionLibraryBase.hpp"

// Template class for varying large numbers of rates in any type of library:
//
// This class takes another type of reaction library as a template parameter 
// and inherits from it. It overrides the DoCalculate rates method of the base 
// library by first calling the base library DoCalculate rates method. It then 
// gets the rates and inverse rates from the base library and multiplies them
// by the specified factors and stores them in its own arrays. It also overrides 
// Rates() and InverseRates() to return its own scaled arrays of rates. Inverse 
// rates are scaled by the same factor as forward rates to preserve detailed 
// balance. I ended up using templates instead of storing the base library with 
// a pointer to get around having to make the the base library classes friends 
// with MultiplierReactionLibrary. 

template<class LibraryType> 
class MultiplierReactionLibrary: public LibraryType {
public:
  MultiplierReactionLibrary(const LibraryType& BaseLibrary) 
      : LibraryType(BaseLibrary), 
      mMultipliers(std::vector<double>(BaseLibrary.NumAllReactions(), 1.0)),
      mRates(BaseLibrary.Rates().size()), 
      mInverseRates(BaseLibrary.InverseRates().size()) {}
  
  std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const {
    return std::unique_ptr<ReactionLibraryBase>(
        new MultiplierReactionLibrary(*this));
  }
  
  std::vector<bool> SetRateMultipliers(const std::vector<Reaction>& reacs, 
      const std::vector<double>& mults) {
    if (mults.size() != reacs.size()) 
      throw std::length_error("MultiplierReactionLibrary: reacs and mults don't have the same size.");
    
    std::vector<bool> found(reacs.size(), false);
    std::vector<double> allMults = mMultipliers.AllData();
    
    for (unsigned int j=0; j<reacs.size(); ++j) { 
      for (unsigned int i=0; i<LibraryType::Reactions().size(); ++i) {
        if (LibraryType::Reactions()[i].IsEquivalent(reacs[j]) || LibraryType::Reactions()[i].Inverse().IsEquivalent(reacs[j])) {
            allMults[i] = mults[j];  
            found[j] = true;
        }
      }
    }  
    std::vector<bool> active = mMultipliers.ActiveFlags(); 
    mMultipliers = allMults; 
    mMultipliers.SetActive(active); 
    return found;
  }
  
  //bool SetSingleRateMultiplier(const Reaction& reacIn, double mult) { 
  //  bool found = false;
  //  std::vector<double> mults = 
  //  for (unsigned int i=0; i<LibraryType::Reactions().size(); ++i) {
  //      if (LibraryType::Reactions()[i].IsEquivalent(reacIn)) {
  //          //mMultipliers[i] = mult;  
  //          found = true;
  //      }
  //  }
  //  return found;
  //}

  std::string Name() const {
    return "Multiply Rates Reaction Library (" + LibraryType::Name() +")";
  }

  const std::vector<double>& Rates() const {
    for (unsigned int i=0; i<LibraryType::Rates().size(); ++i) { 
      //if (mRates[i] != LibraryType::Rates()[i]) 
      //     std::cout << "Rate Difference " << mRates[i] << " " << LibraryType::Rates()[i] << " " << mMultipliers[i] << std::endl;
    }
    //return LibraryType::Rates();
    return mRates;
  }

  const std::vector<double>& InverseRates() const {
    //return LibraryType::InverseRates();
    return mInverseRates;
  }
  
  void CalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const
          pScreeningChemicalPotentialCorrection) { 
    DoCalculateRates(thermoState, partitionFunctionsWithoutSpinTerms, 
      expArgumentCap, pZs, pScreeningChemicalPotentialCorrection);
  }

protected:
  ReactionData<double> mMultipliers; 
  std::vector<double> mRates, mInverseRates;
  
  void DoLoopOverReactionData(
      const std::function<void(GeneralReactionData * const)>& func) {
    
    LibraryType::DoLoopOverReactionData(func); 
    
    auto numActive = mMultipliers.ActiveData().size();
    mRates = std::vector<double>(numActive);
    mInverseRates = std::vector<double>(numActive);

    func(&mMultipliers);
  }
  
  // Have to have both the old and the new call signatures of DoCalculateRates 
  // to allow for backward compatibility. This works because C++ templates never 
  // instantiate a function if it is not called, so 
  // LibraryType::DoCalculateRates(...) not existing does not matter as long as 
  // DoCalculateRates with the same call signature is never called.   
  void DoCalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const
          pScreeningChemicalPotentialCorrection) { 
    LibraryType::DoCalculateRates(thermoState, partitionFunctionsWithoutSpinTerms, 
        expArgumentCap, pZs, pScreeningChemicalPotentialCorrection); 
    MultiplyRates();
  }
  
  void DoCalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap) {
    LibraryType::DoCalculateRates(thermoState, partitionFunctionsWithoutSpinTerms, 
        expArgumentCap);
    MultiplyRates();
  }
 
  void MultiplyRates() {    
    if (mMultipliers.size() != LibraryType::Rates().size() || 
      mRates.size() != LibraryType::Rates().size() ||
      mInverseRates.size() != LibraryType::InverseRates().size()) 
        std::cout << mMultipliers.size() << " " 
            << LibraryType::Rates().size() << " " << mRates.size() << " "  
            << LibraryType::InverseRates().size() << " " << mInverseRates.size() 
            << std::endl;
 
    for (unsigned int i=0; i<mRates.size(); ++i) 
        mRates[i] = mMultipliers[i]*LibraryType::Rates()[i];
    for (unsigned int i=0; i<mInverseRates.size(); ++i) 
        mInverseRates[i] = mMultipliers[i]*LibraryType::InverseRates()[i];
  }

};

#endif // SKYNET_REACTIONS_MULTIPLIERREACTIONLIBRARY_HPP_
