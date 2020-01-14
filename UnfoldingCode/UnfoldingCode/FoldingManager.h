#ifndef FOLDINGMANAGER_H_
#define FOLDINGMANAGER_H_

class FoldingManager {
public:
    explicit FoldingManager() = default;

    ~FoldingManager();

    FoldingManager(const FoldingManager& rhs) = delete;
    FoldingManager(FoldingManager&& rhs) = delete;
    FoldingManager& operator=(const FoldingManager& rhs) = delete;
    FoldingManager& operator=(FoldingManager&& rhs) = delete;
        
};

#endif
