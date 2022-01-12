#include "find_mismatches.h"
#include "test/test_framework.h"
#include "RangeIterator.h"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cassert>
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

void modifyto4bit(std::string & s){
    for(char & c : s){
        c = to4bit(c);
    }
}

std::vector<match> find_matches_gold_local(std::string & genomeb4, std::vector<std::string> & patternsb4, int max_mismatches){
    size_t pattern_size = patternsb4.at(0).size();
    for(std::string & p : patternsb4){
        assert(p.size() == pattern_size);
    }
    std::vector<match> matches;
    for(size_t i : range(genomeb4.size() - pattern_size + 1)){
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
                matches.push_back(match{.loc=i, .mismatches=mismatches, .pattern_idx=j});
            }
        }
    }
    return matches;
}

std::vector<match> find_matches_gold(std::string genome, std::vector<std::string> patterns, int max_mismatches){
    modifyto4bit(genome);
    for(std::string & pat : patterns){
        modifyto4bit(pat);
    }
    return find_matches_gold_local(genome, patterns, max_mismatches);
}

std::vector<match> find_matches(std::string genome, std::vector<std::string> patterns, int max_mismatches){
    return find_matches_gold(genome, patterns, max_mismatches);
}

void sort_matches(std::vector<match> & matches){
    std::sort(matches.begin(), matches.end(), [](match & a, match & b){
        if(a.loc < b.loc) return true;
        if(a.loc > b.loc) return false;
        if(a.pattern_idx < b.pattern_idx) return true;
        if(a.pattern_idx > b.pattern_idx) return false;
        return false;
   });
}
bool matches_equal(std::vector<match> & m1, std::vector<match> & m2){
    return m1.size() == m2.size() && 
        std::equal(m1.begin(), m1.end(), m2.begin(), [](match a, match b){
            return 
                a.loc == b.loc &&
                a.mismatches == b.mismatches &&
                a.pattern_idx == b.pattern_idx;
        });
}

void atomic_print_match(match & m){
    // thread save printing: whole line guarenteed to print at once
    std::ostringstream oss;
    oss << m.loc << "\t" << m.mismatches << "\t" << m.pattern_idx << "\n";
    std::cout << oss.str();
}

TEST(test_find_matches_gold){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACG",
        "CGTAGCTAG",
        "CAGTCGATC"
    };
    int mismatches = 3;
    std::vector<match> expected = {
        match{
            .loc=2,
            .mismatches=0,
            .pattern_idx=0,
        },
        match{
            .loc=21,
            .mismatches=0,
            .pattern_idx=1,
        },
        match{
            .loc=5,
            .mismatches=2,
            .pattern_idx=2,
        },
        match{
            .loc=13,
            .mismatches=0,
            .pattern_idx=2,
        },
    };
    std::vector<match> actual = find_matches_gold(genome, patterns, mismatches);
    sort_matches(actual);
    sort_matches(expected);
    for(match m : actual){
        atomic_print_match(m);
    }
    return matches_equal(actual, expected);
}
TEST(find_mismatches_gold_perf){
    std::vector<std::string> patterns(25, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for(int i : range(100)){
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 5;
    int start = clock();
    std::vector<match> actual = find_matches_gold(genome, patterns, mismatches);
    int end = clock();
    std::cout << "time: " << (end - start)/double(CLOCKS_PER_SEC) << "\n";
    return true;
}
