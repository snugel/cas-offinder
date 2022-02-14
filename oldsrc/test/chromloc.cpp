#include "chromloc.h"
#include "test/test_framework.h"
#include <sstream>
#include <algorithm>

static std::string valid_chrom_pos_file = 
"Name\tLoc\n"
"chr1\t54\n"
"chr2\t89\n"
"chr51\t1099511627776\n"
;

static std::vector<chromloc> valid_file_result = {
    chromloc{.name="chr1",.loc=54},
    chromloc{.name="chr2",.loc=89},
    chromloc{.name="chr51",.loc=1099511627776},
};

static std::string chrom_pos_file_out_of_order = 
"Name\tLoc\n"
"chr1\t54\n"
"chr2\t89\n"
"chr51\t29\n"
;

static std::string chrom_pos_file_additional_data = 
"Name\tLoc\n"
"chr1\t54\n"
"chr2\t89\textra\n"
"chr51\t129\n"
;

TEST(test_parse_chrom_loc_valid){
    std::stringstream strstm(valid_chrom_pos_file);
    std::vector<chromloc> result = parse_chromloc_file(strstm);
    return result.size() == valid_file_result.size() &&
        std::equal(valid_file_result.begin(), valid_file_result.end(), result.begin(), [](chromloc a, chromloc b){
            return a.loc == b.loc && a.name == b.name;
        });
}

TEST(test_write_chrom_loc_valid){
    std::stringstream strstm;
    write_chromloc_file(valid_file_result, strstm);
    return strstm.str() == valid_chrom_pos_file;
}


TEST(test_parse_chrom_loc_out_of_order){
    try{
        std::stringstream strstm(chrom_pos_file_out_of_order);
        std::vector<chromloc> result = parse_chromloc_file(strstm);
        return false;
    }
    catch(std::runtime_error){
        return true;
    }
}

TEST(test_parse_chrom_loc_additional_data){
    try{
        std::stringstream strstm(chrom_pos_file_additional_data);
        std::vector<chromloc> result = parse_chromloc_file(strstm);
        return false;
    }
    catch(std::runtime_error){
        return true;
    }
}