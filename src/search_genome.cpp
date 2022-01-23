
#include "search_genome.h"
#include "RangeIterator.h"
#include "chromloc.h"
#include "postprocess.h"
#include "find_mismatches.h"
#include "bit4ops.h"

#include <iostream>
#include <fstream>
#include <thread>
#include <cassert>
#include <vector>

using namespace std;

void search_genome(
    std::vector<std::string> patterns, 
    std::vector<chromloc> chrom_locs, 
    int mismatches, 
    Channel<GenomeInput> * inp_channel, 
    Channel<GenomeMatch> * out_channel
){  
    size_t pattern_size = patterns.at(0).size();

    // std::vector<uint32_t> data;
    GenomeInput inp;
    // while(inp_channel->receive(inp)){
    //     data.insert(data.end(), inp.data.get(), inp.data.get() + inp.size);
    // }
    Channel<WorkerOutput> match_channel;
    std::thread find_matches_thread(find_matches_worker,inp_channel,patterns,mismatches,&match_channel);

	// std::cerr << "Searching genome..." << data.size()*8 << std::endl;
    
	std::vector<uint64_t> chrom_poses(chrom_locs.size());
	for(size_t i : range(chrom_locs.size())){
		chrom_poses[i] = chrom_locs[i].loc;
	}

    WorkerOutput outs;
    constexpr int bit4_c = 8;
    while(match_channel.receive(outs)){
        // filter out and output matches
        for(size_t i : range(outs.num_matches)){
            match m = outs.matches.get()[i];
            if(m.loc < outs.input.size * bit4_c){
                size_t true_loc = m.loc + outs.input.idx * bit4_c;
                if(!crosses_chrom_wall(true_loc, chrom_poses, pattern_size)){
                    chrom_info chr_info = get_chrom_info(true_loc, chrom_poses);
                    chromloc chrloc = chrom_locs.at(chr_info.chrom_idx);
                    out_channel->send(GenomeMatch{
                        .dna_match=bit4tostr(outs.input.data.get(), m.loc, m.loc + pattern_size),
                        .cromosome=chrloc.name,
                        .chrom_loc=chr_info.rel_loc,
                        .pattern_idx=m.pattern_idx,
                        .mismatches=m.mismatches,
                    });
                }
            }
        }
    }
	std::cerr << "terminating finding\n"; 
    out_channel->terminate();
    find_matches_thread.join();
}
