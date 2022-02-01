#include "read_twobit.h"
#include "read_fasta.h"

#include <iostream>
#include <sys/stat.h>
#include <cstdlib>
#include <dirent.h>

using namespace std;

int check_file(const char* filename) {
	struct stat file_stat;
	return (stat(filename, &file_stat) == 0);
}

int read_file(std::string &filepath, std::vector<std::string> &chrnames, Channel<FileChunk> &content, std::vector<uint64_t> &chrpos){
    if(read_twobit(filepath, chrnames, content, chrpos)){
        if(read_fasta(filepath, chrnames, content, chrpos)){
            std::cerr << "skipping invalid file: " << filepath << std::endl;
            return 1;
        }
    }
    return 0;
}

int read_folder(std::string &filepath, std::vector<std::string> &chrnames, Channel<FileChunk> &content, std::vector<uint64_t> &chrpos){
	DIR* dir;
	dirent *ent;
	if ((dir = opendir(filepath.c_str())) == NULL) {
		if (check_file(filepath.c_str())) {
			return read_file(filepath, chrnames, content, chrpos);
		} else {
            std::cerr << "An error has occured while opening directory: " << filepath << std::endl; 
			exit(1);
		}
	} else {
        int read_files_count = 0;
		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG || ent->d_type == DT_LNK) {
                if(filepath.back() == '/'){
                    filepath.pop_back();
                }
				filepath = filepath + "/" + ent->d_name;
			    read_files_count += read_file(filepath, chrnames, content, chrpos);
			}
		}
        return !read_files_count;
    }
}
