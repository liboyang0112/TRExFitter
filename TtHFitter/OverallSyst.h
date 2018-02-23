#ifndef OVERALLSYST_H
#define OVERALLSYST_H

class OverallSyst {
public:
    OverallSyst(std::string name, float up, float down);
    ~OverallSyst();

    std::string fName;
    float fUp;
    float fDown;
};

#endif
