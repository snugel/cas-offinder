#include "find_mismatches.h"
#include <algorithm>
#include <sstream>
#include <iostream>

void sort_matches(std::vector<match> & matches){
    std::sort(matches.begin(), matches.end(), [](match & a, match & b){
            if(a.pattern_idx < b.pattern_idx) return true;
            if(a.pattern_idx > b.pattern_idx) return false;
            if(a.loc < b.loc) return true;
            if(a.loc > b.loc) return false;
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