#pragma once
#include <string>
#include <vector>
#include "channel.h"
#include "find_mismatches.h"

void chunkify_data(Channel<std::string> * input_channel, Channel<GenomeInput> * genome_input, int pattern_size, int chunk_size);
