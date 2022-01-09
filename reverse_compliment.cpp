#include "reverse_compliment.h"
#include "test/test_framework.h"


std::string reverse_compliment(std::string s){
    std::string reversed(s.rbegin(), s.rend());
    for(char & c : reversed){
        switch (c)
        {
        case 'G':
            c = 'C';
            break;
        case 'C':
            c = 'G';
            break;
        case 'A':
            c = 'T';
            break;
        case 'T':
            c = 'A';
            break;
        }
    }
    return reversed;
}

TEST(test_reverse_compliment){
    std::string str = "ATTCGA";
    std::string reversed = "TCGAAT";

    return reverse_compliment(str) == reversed;
}