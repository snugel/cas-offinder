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

static void to_upper_s(string& s)
{
    for (char& c : s) {
        c = to_upper(c);
    }
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
    string pattern;
    string line;
    vector<string> sline;
    string compares;
    vector<string> m_ids;
    vector<int> m_thresholds;
    parseLine(input, genomedir);
    parseLine(input, pattern);
    to_upper_s(pattern);

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
            to_upper_s(sline[0]);
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
