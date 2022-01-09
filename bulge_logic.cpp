#include "bulge_logic.h"
#include "RangeIterator.h"


void get_bulge_info(bulge_info & result, std::string & genome, int bulge_num, int num_orig_patterns, int dna_bulges, int rna_bulges){

}
std::vector<std::string> augment_patterns_with_bulges(std::vector<std::string> patterns, int dna_bulges, int rna_bulges){
    std::vector<std::string> augmented;
    int pad_size = dna_bulges;
    for(std::string orig : patterns){
        augmented.push_back(std::string(pad_size, 'N') + orig);
        for(int i : range(1, dna_bulges + 1)){
            for(int j : range(1, orig.size() - 1)){
                augmented.push_back(std::string(pad_size-i, 'N') + std::string(orig.begin(), orig.begin() + j) + std::string(i, 'N') + std::string(orig.begin() + j, orig.end()));
            }
        }
        for(int i : range(1, rna_bulges + 1)){
            for(int j : range(i, orig.size())){
                augmented.push_back(std::string(pad_size+i, 'N') + std::string(orig.begin(), orig.begin() + j - i) + std::string(orig.begin() + j, orig.end()));
            }
        }
    }
}

TEST(test_bulge_logic){
    std::string pattern = "ACT";
    std::vector<std::string> augmented = augment_patterns_with_bulges({pattern}, 2, 2);

}