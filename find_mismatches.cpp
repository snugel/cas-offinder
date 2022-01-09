#include "find_mismatches.h"
#include "test/test_framework.h"
#include "reverse_compliment.h"
#include "RangeIterator.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <ctime>




char to4bit(char c){
    constexpr char G = 0x1;
    constexpr char C = 0x2;
    constexpr char A = 0x4;
    constexpr char T = 0x8;
    switch (c)
    {
    case 'G':
        return G;
    case 'C':
        return C;
    case 'A':
        return A;
    case 'T':
        return T;
    case 'R':
        return A | G;
    case 'Y':
        return C | T;
    case 'S':
        return G | C;
    case 'W':
        return A | T;
    case 'K':
        return G | T;
    case 'M':
        return A | C;
    case 'B':
        return C | G | T;
    case 'D':
        return A | G | T;
    case 'H':
        return A | C | T;
    case 'V':
        return A | C | G;
    case 'N':
        return A | C | G | T;
    default:
        throw std::runtime_error("got unexpected letter, input must be mixed base");
    }
}

char complimentb4(char c){
    char left = (c & 0xb) >> 1;
    char right = (c & 0x5) << 1;
    return left | right;
}

void modifyto4bit(std::string & s){
    for(char & c : s){
        c = to4bit(c);
    }
}

struct local_match{
    size_t loc;
    int mismatches;
    size_t pattern_idx;
};
std::vector<local_match> find_matches_gold_local(std::string & genomeb4, std::vector<std::string> & patternsb4, int max_mismatches){
    size_t pattern_size = patternsb4.at(0).size();
    for(std::string & p : patternsb4){
        assert(p.size() == pattern_size);
    }
    std::vector<local_match> matches;
    for(size_t i : range(genomeb4.size())){
        for(size_t j : range(patternsb4.size())){
            std::string & pattern = patternsb4[j];
            int mismatches = 0;
            for(size_t j : range(pattern_size)){
                if(!(genomeb4[j+i] & pattern[j])){
                    mismatches++;
                    if(mismatches > max_mismatches){
                        break;
                    }
                }
            }
            if(mismatches <= max_mismatches){
                matches.push_back(local_match{.loc=i, .mismatches=mismatches, .pattern_idx=j});
            }
        }
    }
    return matches;
}
void add_matches(std::vector<match> & matches, std::vector<local_match> & local_matches, char dir){
    for(local_match lm : local_matches){
        matches.push_back(match{
            .loc=lm.loc,
            .mismatches=lm.mismatches,
            .pattern_idx=lm.pattern_idx,
            .direction=dir,
        });
    }
}
std::vector<match> find_matches_gold(std::string genome, std::vector<std::string> patterns, int max_mismatches){
    modifyto4bit(genome);
    std::vector<std::string> plus_patterns = patterns;
    for(std::string & pat : plus_patterns){
        modifyto4bit(pat);
    }
    std::vector<std::string> minus_patterns = patterns;
    for(std::string & pat : minus_patterns){
        pat = reverse_compliment(pat);
        modifyto4bit(pat);
    }
    std::vector<local_match> plus_local_matchs = find_matches_gold_local(genome, plus_patterns, max_mismatches);
    std::vector<local_match> minus_local_matchs = find_matches_gold_local(genome, minus_patterns, max_mismatches);
    std::vector<match> matches;
    add_matches(matches, plus_local_matchs, '+');
    add_matches(matches, minus_local_matchs, '-');
    return matches;
}
std::vector<match> find_matches(std::string genome, std::vector<std::string> patterns, int max_mismatches){
    return find_matches_gold(genome, patterns, max_mismatches);
}

void sort_matches(std::vector<match> & matches){
    std::sort(matches.begin(), matches.end(), [](match & a, match & b){
        return 
            a.loc < b.loc && 
            a.pattern_idx < b.pattern_idx && 
            a.direction < b.direction;
    });
}
TEST(test_find_matches_gold){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACG",
        "CGTAGCTAG",
        "GATCGACTG"
    };
    int mismatches = 3;
    std::vector<match> expected = {
        match{
            .loc=2,
            .mismatches=0,
            .pattern_idx=0,
            .direction='+',
        },
        match{
            .loc=21,
            .mismatches=0,
            .pattern_idx=1,
            .direction='+',
        },
        match{
            .loc=25,
            .mismatches=3,
            .pattern_idx=2,
            .direction='+',
        },
        match{
            .loc=3,
            .mismatches=3,
            .pattern_idx=0,
            .direction='-',
        },
        match{
            .loc=5,
            .mismatches=2,
            .pattern_idx=2,
            .direction='-',
        },
        match{
            .loc=13,
            .mismatches=0,
            .pattern_idx=2,
            .direction='-',
        },
        match{
            .loc=22,
            .mismatches=3,
            .pattern_idx=1,
            .direction='-',
        },
    };
    std::vector<match> actual = find_matches_gold(genome, patterns, mismatches);
    sort_matches(actual);
    sort_matches(expected);
    for(match m : actual){
        std::cout << m.loc << "\t" << m.mismatches << "\t" << m.pattern_idx << "\t" << m.direction << "\n";
    }

    if(actual.size() != expected.size()){
        return false;
    }
    return std::equal(actual.begin(), actual.end(), expected.begin(), [](match a, match b){
        return 
            a.direction == b.direction &&
            a.loc == b.loc &&
            a.mismatches == b.mismatches &&
            a.pattern_idx == b.pattern_idx;
    });
}
TEST(find_mismatches_perf){
    std::vector<std::string> patterns(25, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for(int i = 0; i < 10000; i++){
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 5;
    int start = clock();
    std::vector<match> actual = find_matches_gold(genome, patterns, mismatches);
    int end = clock();
    std::cout << "time: " << (end - start)/double(CLOCKS_PER_SEC) << "\n";
    return true;
}


TEST(test_complimentb4){
    return (
        (complimentb4(to4bit('T')) == to4bit('A')) &&
        (complimentb4(to4bit('A')) == to4bit('T')) &&
        (complimentb4(to4bit('G')) == to4bit('C')) &&
        (complimentb4(to4bit('C')) == to4bit('G')) &&
        (complimentb4(to4bit('N')) == to4bit('N')) 
    );
}