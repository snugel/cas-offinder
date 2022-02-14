
#include "RangeIterator.h"
#include "bit4ops.h"
#include "bulge_logic.h"
#include "chromloc.h"
#include "find_mismatches.h"
#include "postprocess.h"
#include "read_fasta.h"
#include "read_genome.h"
#include "search_genome.h"
#include "timing.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <cassert>
#include <cstdarg>
#include <cstdlib>
#include <cstring>

using namespace std;

int main(int argc, char* argv[])
{
    constexpr size_t num_pos_args = 5;
    if (argc < num_pos_args) {
        std::cerr << "Needs at least 4 command line arguments: data_folder, "
                     "mismatches, device, queries\n";
        exit(1);
    }
    double calc_time = time_spent([&]() {
        std::string data_folder(argv[1]);
        int mismatches = atoi(argv[2]);
        char device = *argv[3];
        char** pattern_chrs = &argv[4];
        int num_patterns = argc - num_pos_args + 1;

        std::cerr << "Number of patterns: " << num_patterns << std::endl;
        std::vector<std::string> patterns(num_patterns * 2);
        for (size_t i : range(num_patterns)) {
            std::string pattern(pattern_chrs[i]);
            while (pattern.back() == '\r')
                pattern.pop_back();
            patterns[i * 2] = pattern;
            patterns[i * 2 + 1] = reverse_compliment(pattern);
        }

        if (data_folder.back() != '/') {
            data_folder.push_back('/');
        }

        std::ifstream chrom_locs_file(data_folder + "chrom_locs.csv");
        std::vector<chromloc> chrom_locs = parse_chromloc_file(chrom_locs_file);

        Channel<GenomeInput> data_input_channel(6);
        Channel<GenomeMatch> data_output_channel;
        std::thread reader_thread(
          read_genome, data_folder, patterns.at(0).size(), &data_input_channel);
        std::thread searcher_thread(search_genome,
                                    patterns,
                                    chrom_locs,
                                    mismatches,
                                    &data_input_channel,
                                    &data_output_channel);

        // GGCCGACCTGTCGCTGACGCNNN chr8    517739       GcCCtgCaTGTgGCTGACGCAGG
        // +       5

        std::cout << "RNA\tChromosome\tLocation\tDNA\tDirection\tMismatches\n";

        GenomeMatch m;
        while (data_output_channel.receive(m)) {
            char dir;
            std::string pattern(patterns[m.pattern_idx]);
            if (m.pattern_idx % 2 == 1) {
                m.pattern_idx /= 2;
                m.dna_match = reverse_compliment(m.dna_match);
                pattern = reverse_compliment(pattern);
                dir = '-';
            } else {
                dir = '+';
            }
            indicate_mismatches_dna(m.dna_match, pattern);

            std::cout << pattern << "\t" << m.cromosome << "\t" << m.chrom_loc
                      << "\t" << m.dna_match << "\t" << dir << "\t"
                      << m.mismatches << "\n";
        }
        reader_thread.join();
        searcher_thread.join();
    });
    std::cerr << "Time spent: " << calc_time << std::endl;
    return 0;
}
