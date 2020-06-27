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
}
