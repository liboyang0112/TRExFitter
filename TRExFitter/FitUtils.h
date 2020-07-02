#ifndef FITUTILS_H
#define FITUTILS_H

#include <map>
#include <memory>
#include <vector>

class FittingTool;
class NormFactor;
class RooWorkspace;
class RooSimultaneous;

namespace FitUtils {

/**
 * Implementing external constraints on the workspace
 * @param workspace
 * @param Fitting tool
 * @param the PDF
 * @param vector of norma factors
 */
void ApplyExternalConstraints(RooWorkspace* ws,
                              FittingTool* fitTool,
                              RooSimultaneous* simPdf,
                              const std::vector<std::shared_ptr<NormFactor> >& normFactors);

/**
 * A helper function to set BinnedLikelihoodOptimisation
 * @param workspace
 */
void SetBinnedLikelihoodOptimisation(RooWorkspace* ws);

/**
 * A helper function to injects values of NPs
 * @param workspace
 * @param map of NP values
 */
void InjectGlobalObservables(RooWorkspace* ws, const std::map< std::string, double >& npValues);

/**
 * A helper function to set saturated model shapefactors to const
 * @param workspace
 */
void DisableSaturatedModel(RooWorkspace* ws);

/**
 * A helper function to set POI in a file
 * @param path to the file
 * @param POI
 */
void SetPOIinFile(const std::string& path, const std::string& poi);
}

#endif
