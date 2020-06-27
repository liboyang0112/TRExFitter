#ifndef RANKINGMANAGER_H
#define RANKINGMANAGER_H

#include <string>

class RankingManager {
public:

    explicit RankingManager();
    ~RankingManager();
    RankingManager(const RankingManager& r) = delete;
    RankingManager(RankingManager&& r) = delete;
    RankingManager& operator=(const RankingManager& r) = delete;
    RankingManager& operator=(RankingManager&& r) = delete;

    inline void AddOutputPath(const std::string& path){fOutputPath = path;}

private:

    std::string fOutputPath;
};

#endif
