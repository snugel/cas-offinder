#include "test/test_framework.h"
#include "async_search.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>

using namespace std;

static std::vector<GenomeMatch> async_search_results;
static int intcmp(int i1, int i2){
    return i1 < i2 ? -1 : (i1 > i2 ? 1 : 0);
}
TEST(test_async_search){
    async_search_results.clear();

    constexpr size_t pattern_size = 8;
    constexpr size_t num_patterns = 2;
    constexpr size_t mismatches = 1;
    const char compares[num_patterns*pattern_size+1] =
        "ATRGCNGC"
        "CGNCGRTA"
    ;
    BlockConfig config{
        .OUT_CHUNK_SIZE=1<<16,
        .MIN_CMPS_PER_OUT=1<<16,
        .DEFAULT_CHUNK_BYTES=16,
    };
    auto callback = [](const GenomeMatch*gm){
        async_search_results.push_back(GenomeMatch{
                                           .dna_match=strdup(gm->dna_match),
                                           .chrom_name=strdup(gm->chrom_name),
                                           .chrom_loc=gm->chrom_loc,
                                           .pattern_idx=gm->pattern_idx,
                                           .mismatches=gm->mismatches,
                                       });
    };
    async_search(
                "../cas-offinder/test/t_data/",
                CPU,
                compares,
                pattern_size,
                num_patterns,
                mismatches,
                &config,
                callback
                );
    sort(async_search_results.begin(),async_search_results.end(),[](GenomeMatch gm1, GenomeMatch gm2){
        int res1 = strcmp(gm1.chrom_name, gm2.chrom_name);
        int res2 = intcmp(gm1.chrom_loc, gm2.chrom_loc);
        int res3 = intcmp(gm1.pattern_idx, gm2.pattern_idx);
        return res1 < 0 ||
                (res1 == 0 && res2 < 0) ||
                (res2 == 0 && res3 < 0);
    });
    ifstream ref_file("../cas-offinder/test/t_data/async_search_out.txt");
    t_assert(ref_file);
    stringstream ref_data;
    ref_data << ref_file.rdbuf();
    stringstream act_data;
    for(GenomeMatch gm : async_search_results){
        act_data << gm.pattern_idx <<'\t'  << gm.chrom_name<<'\t' << gm.chrom_loc<<'\t'<< gm.dna_match<<'\t'<< gm.mismatches << '\n';
    }
    //ref_file.flush();
    string act_data_str = act_data.str();
    string ref_data_str = ref_data.str();
    cerr << act_data_str << "\n\n";
    cerr << ref_data_str << "\n\n";
    return act_data_str == ref_data_str;
}
