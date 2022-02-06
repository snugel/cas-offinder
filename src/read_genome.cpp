#include "read_genome.h"

#include <cstdlib>
#include <iostream>
#include <sys/stat.h>
#ifdef WIN32
#include "ext/dirent.h"
#else
#include <dirent.h>
#endif

struct FolderReader
{
    TwoBitReader* r2bit = nullptr;
    FastaReader* rfasta = nullptr;
    std::string filepath;
    DIR* dir = nullptr;
};

static int check_file(const char* filename)
{
    struct stat file_stat;
    return (stat(filename, &file_stat) == 0);
}
static void get_new_file(FolderReader* reader, const char* filepath)
{
    TwoBitReader* twobit_file = create_2bit_reader(filepath);
    if (twobit_file) {
        reader->r2bit = twobit_file;
        reader->rfasta = nullptr;
        return;
    }
    FastaReader* fasta_file = create_fasta_reader(filepath);
    if (fasta_file) {
        reader->r2bit = nullptr;
        reader->rfasta = fasta_file;
        return;
    }
    reader->r2bit = nullptr;
    reader->rfasta = nullptr;
}
FolderReader* create_folder_reader(const char* path)
{
    FolderReader* reader = new FolderReader();
    reader->filepath = std::string(path);
    if ((reader->dir = opendir(path)) == NULL) {
        if (check_file(path)) {
            get_new_file(reader, path);
        } else {
            return nullptr;
        }
    }
    return reader;
}
ChromData read_next_folder(FolderReader* reader)
{
    ChromData data{
        .name = nullptr,
        .bit4data = nullptr,
        .n_nucl = 0,
    };
    if (reader->r2bit) {
        data = read_next_2bit(reader->r2bit);
        if (!data.name) {
            free_2bit_reader(&reader->r2bit);
            reader->r2bit = nullptr;
        }
    }
    if (!data.name && reader->rfasta) {
        data = read_next_fasta(reader->rfasta);
        if (!data.name) {
            free_fasta_reader(&reader->rfasta);
            reader->rfasta = nullptr;
        }
    }
    if (!data.name && reader->dir) {
        dirent* ent;
        std::string filepath = reader->filepath;
        while ((ent = readdir(reader->dir)) != NULL) {
            if (ent->d_type == DT_REG || ent->d_type == DT_LNK) {
                if (filepath.back() == '/') {
                    filepath.pop_back();
                }
                filepath = filepath + "/" + ent->d_name;
                get_new_file(reader, filepath.c_str());
                data = read_next_folder(reader);
                break;
            }
        }
    }
    return data;
}
void free_folder_reader(FolderReader** fptr)
{
    FolderReader* reader = *fptr;
    if (reader->r2bit) {
        free_2bit_reader(&reader->r2bit);
    }
    if (reader->rfasta) {
        free_fasta_reader(&reader->rfasta);
    }
    closedir(reader->dir);
    delete reader;
    *fptr = nullptr;
}
