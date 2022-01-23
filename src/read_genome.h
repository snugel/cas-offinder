#pragma once
#include "channel.h"
#include "find_mismatches.h"

void read_genome(std::string data_folder, size_t pattern_size, Channel<GenomeInput> * data_output_p);
