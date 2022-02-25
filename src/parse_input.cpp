#include "parse_input.h"
#include "RangeIterator.h"
#include "bit4ops.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

using namespace std;

static inline bool parseLine(istream& input,
                             string& line,
                             bool expect_line = true)
{
    if (expect_line && input.eof()) {
        throw runtime_error("Unexpected end of file.");
    } else {
        if (!getline(input, line))
            return false;

        // Remove possible Windows-only '\r'.
        if (line[line.size() - 1] == '\r')
            line.pop_back();
    }
    return true;
}

static char* to_chararr(string s)
{
    char* a = (char*)malloc(s.size() + 1);
    memcpy(a, s.data(), s.size());
    a[s.size()] = 0;
    return a;
}

static vector<string> split(string const& input)
{
    istringstream sbuffer(input);
    vector<string> ret((istream_iterator<string>(sbuffer)),
                       istream_iterator<string>());
    return ret;
}

static InFileInfo parse_input(istream& input)
{
    string genomedir;
    string pattern_line;
    string pattern;
    string line;
    vector<string> sline;
    string compares;
    int dna_bulges = 0;
    int rna_bulges = 0;
    vector<string> m_ids;
    vector<int> m_thresholds;
    parseLine(input, genomedir);
    parseLine(input, pattern_line);
    sline = split(pattern_line);
    pattern = sline.at(0);
    if(sline.size() > 1){
        dna_bulges = stoi(sline[1]);
    }
    if(sline.size() > 2){
        rna_bulges = stoi(sline[2]);
    }
    if(sline.size() > 3){
        throw runtime_error(
          "Second line of input file can have at most 3 elements: pattern, dna_bulges, rna_bulges");
    }

    try {
        size_t entrycnt = 0;
        while (getline(input, line)) {
            if (line.empty())
                break;
            if (line[line.length() - 1] == '\r')
                line = line.substr(0, line.length() - 1);
            sline = split(line);
            if (sline.size() != 2 && sline.size() != 3) {
                cerr << "Skipping malformed input guide line." << endl;
                break;
            }
            if (sline[0].length() != pattern.length()) {
                throw runtime_error(
                  "The length of target sequences should match with the length "
                  "of pattern sequence.");
            }
            compares += sline[0];
            m_thresholds.push_back(atoi(sline[1].c_str()));
            if (entrycnt == 0) {
                entrycnt = sline.size();
            } else if (entrycnt != sline.size()) {
                throw runtime_error(
                  "The number of entries below 2nd line should be consistent.");
            }
            if (sline.size() == 3)
                m_ids.push_back(sline[2]);
        }
        if (m_thresholds.size() == 0) {
            throw runtime_error("Must have at least one pattern in below the "
                                "2nd line to search.");
        }
        int threth = m_thresholds.back();
        for (int i : m_thresholds) {
            if (i != threth) {
                throw runtime_error(
                  "This version of cas-offinder requires that the number of "
                  "mismatches be conssitent across all lines.");
            }
            if (i < 0) {
                throw runtime_error("Number of mismatches must be positive.");
            }
        }
    } catch (const exception& e) {
        cerr << "Critical error! " << e.what() << endl;
        exit(1);
    }
    char** ids = nullptr;
    if (m_ids.size()) {
        ids = (char**)malloc(sizeof(ids[0]) * m_ids.size());
        for (int i : range(m_ids.size())) {
            ids[i] = to_chararr(m_ids[i]);
        }
    }
    return InFileInfo{
        .genome_path = to_chararr(genomedir),
        .pattern = to_chararr(pattern),
        .compares = to_chararr(compares),
        .pattern_size = pattern.size(),
        .num_patterns = m_thresholds.size(),
                .dna_bulges=dna_bulges,
                .rna_bulges=rna_bulges,
        .ids = ids,
        .mismatches = uint32_t(m_thresholds.back()),
    };
}

InFileInfo read_file(const char* inputfile)
{
    if (strlen(inputfile) == 1 && inputfile[0] == '-') {
        return parse_input(cin);
    } else {
        ifstream fi(inputfile);
        if (!fi) {
            cerr << "input file: " << inputfile << " not found\n";
            exit(1);
        }
        return parse_input(fi);
    }
}
