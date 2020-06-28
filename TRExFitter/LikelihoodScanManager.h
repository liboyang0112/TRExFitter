#ifndef LIKELIHOOHSVANMANAGER_H
#define LIKELIHOOHSVANMANAGER_H

class LikelihoodScanManager {
public:

    explicit LikelihoodScanManager();
    ~LikelihoodScanManager();

    LikelihoodScanManager(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager(LikelihoodScanManager&& m) = delete;
    LikelihoodScanManager& operator=(const LikelihoodScanManager& m) = delete;
    LikelihoodScanManager& operator=(LikelihoodScanManager&& m) = delete;

private:

};

#endif
