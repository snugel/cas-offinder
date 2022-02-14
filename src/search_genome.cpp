#include "chunkify_data.h"
#include "search_genome.h"
#include "RangeIterator.h"
#include "chromloc.h"
#include "postprocess.h"
#include "find_mismatches.h"
#include "bit4ops.h"
#include "read_folder.h"

#include <iostream>
#include <fstream>
#include <thread>
#include <cassert>
#include <vector>

using namespace std;

void read_folder_t(std::string * filepath,  std::vector<std::string> * chrnames, Channel<std::string> * content, std::vector<uint64_t> * chrpos){
    
	if (read_folder(*filepath, *chrnames, *content, *chrpos)) {
        cerr << "Non-acceptable file/folder: " << *filepath << endl;
		exit(1);
	}
    content->terminate();
}

void search_genome(
    std::vector<std::string> patterns, 
    std::string path,
    int mismatches, 
    Channel<GenomeMatch> * out_channel
){  
    size_t pattern_size = patterns.at(0).size();
    
    Channel<string> file_read_data;
    std::vector<std::string> chrnames;
    std::vector<uint64_t>  chrpos;
    thread file_read_thread(read_folder_t, &path, &chrnames, &file_read_data, &chrpos); 
    
    int chunk_size = (1<<24) * 8;
    Channel<GenomeInput> genome_input;
    thread chunkify_thread(chunkify_data, &file_read_data, &genome_input, pattern_size, chunk_size);
     
    
    // std::vector<uint32_t> data;
    GenomeInput inp;
    // while(inp_channel->receive(inp)){
    //     data.insert(data.end(), inp.data.get(), inp.data.get() + inp.size);
    // }
    Channel<WorkerOutput> match_channel;
    std::thread find_matches_thread(find_matches_worker,&genome_input,patterns,mismatches,&match_channel);
       
    // file read needs to finish so that 
    // chrpos is defined...
    file_read_thread.join();
	std::cerr << "Searching genome..." << std::endl;
    
	std::vector<uint64_t> chrom_poses = chrpos;

    WorkerOutput outs;
    constexpr int bit4_c = 8;
    while(match_channel.receive(outs)){
        // filter out and output matches
        for(size_t i : range(outs.num_matches)){
            match m = outs.matches.get()[i];
            if(m.loc < outs.input.size * bit4_c){
                size_t true_loc = m.loc + outs.input.idx * bit4_c;
                if(!crosses_chrom_wall(true_loc, chrpos, pattern_size)){
                    chrom_info chr_info = get_chrom_info(true_loc, chrpos);
                    string chrname = chrnames.at(chr_info.chrom_idx);
                    out_channel->send(GenomeMatch{
                        .dna_match=bit4tostr(outs.input.data.get(), m.loc, m.loc + pattern_size),
                        .cromosome=chrname,
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
    chunkify_thread.join();
    find_matches_thread.join();
}
