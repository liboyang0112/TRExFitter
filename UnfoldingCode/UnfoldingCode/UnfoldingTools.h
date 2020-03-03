#ifndef UNFOLDINGTOOLS_H_
#define UNFOLDINGTOOLS_H_

class TH2;

namespace UnfoldingTools {

    void NormalizeMatrix(TH2* matrix, const bool byRow);

    void Correct2DMatrix(TH2* matrix);

}

#endif
