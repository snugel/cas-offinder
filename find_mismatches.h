#include<vector>
#include<string>


struct match{
    size_t loc;
    int mismatches;
    size_t pattern_idx;
    char direction;
};
std::vector<match> find_matches(std::string genome, std::vector<std::string> patterns, int max_mismatches);
