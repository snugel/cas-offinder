#include "RangeIterator.h"
#include "bit4ops.h"
#include "ceildiv.h"
#include "read_genome.h"
#include "test/test_framework.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

TEST(test_read_twobit)
{
    TwoBitReader* reader =
      create_2bit_reader("../cas-offinder/test/t_data/upstream1000.2bit");
    assert(reader && "test file does not exist!");
    ChromData data;
    vector<uint64_t> chrompos;
    vector<string> chromnames;
    vector<uint8_t> all_data;
    chrompos.push_back(0);
    while (true) {
        data = read_next_2bit(reader);
        if (data.name == nullptr) {
            break;
        }
        chrompos.push_back(chrompos.back() + data.n_nucl);
        chromnames.push_back(data.name);
        all_data.insert(
          all_data.end(), data.bit4data, data.bit4data + cdiv(data.n_nucl, 2));
    }
    free_2bit_reader(&reader);
    // the last position is just the total size of the genome
    uint64_t genome_size = chrompos.back();
    chrompos.pop_back();
    char buf[20000] = { 0 };
    bit42str(buf, all_data.data(), 0, genome_size);
    std::stringstream buffer;
    ifstream file("../cas-offinder/test/t_data/expected.txt");
    buffer << file.rdbuf();
    string expected = buffer.str();
    bool equal = !memcmp(expected.data(), buf, expected.size());
    return equal && chromnames.size() == chrompos.size() &&
           chromnames.size() == 13 && chromnames.at(5).size() > 5;
}

TEST(test_read_fasta)
{
    FastaReader* reader =
      create_fasta_reader("../cas-offinder/test/t_data/upstream1000.fa");
    assert(reader && "test file does not exist!");
    ChromData data;
    vector<uint64_t> chrompos;
    vector<string> chromnames;
    vector<uint8_t> all_data;
    chrompos.push_back(0);
    while (true) {
        data = read_next_fasta(reader);
        if (data.name == nullptr) {
            break;
        }
        chrompos.push_back(chrompos.back() + data.n_nucl);
        chromnames.push_back(data.name);
        all_data.insert(
          all_data.end(), data.bit4data, data.bit4data + cdiv(data.n_nucl, 2));
    }
    free_fasta_reader(&reader);
    // the last position is just the total size of the genome
    uint64_t genome_size = chrompos.back();
    chrompos.pop_back();
    char buf[20000] = { 0 };
    bit42str(buf, all_data.data(), 0, genome_size);
    ifstream file("../cas-offinder/test/t_data/expected.txt");
    std::stringstream buffer;
    buffer << file.rdbuf();
    string s = buffer.str();
    bool equal = !memcmp(s.data(), buf, s.size());
    return equal && chromnames.size() == chrompos.size() &&
           chromnames.size() == 13 && chromnames.at(5).size() > 5;
}

TEST(test_read_folder)
{
    FolderReader* reader = create_folder_reader("../cas-offinder/test/t_data/");
    assert(reader && "test file does not exist!");
    ChromData data;
    vector<uint64_t> chrompos;
    vector<string> chromnames;
    chrompos.push_back(0);
    while (true) {
        data = read_next_folder(reader);
        if (data.name == nullptr) {
            break;
        }
        chrompos.push_back(chrompos.back() + data.n_nucl);
        chromnames.push_back(data.name);
    }
    free_folder_reader(&reader);
    // the last position is just the total size of the genome
    uint64_t genome_size = chrompos.back();
    chrompos.pop_back();
    return genome_size == 25014 && chromnames.size() == chrompos.size() &&
           chromnames.size() == 26 && chromnames.at(20).size() > 5;
}
