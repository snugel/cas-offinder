#include <inttypes.h>
struct FolderReader;
struct ChromData
{
    char* name; // c-string
    uint8_t* bit4data;
    uint64_t n_nucl;
};
FolderReader* create_folder_reader(const char* path);
ChromData read_next_folder(FolderReader* reader);
void free_folder_reader(FolderReader**);

struct TwoBitReader;
TwoBitReader* create_2bit_reader(const char* path);
ChromData read_next_2bit(TwoBitReader* reader);
void free_2bit_reader(TwoBitReader**);

struct FastaReader;
FastaReader* create_fasta_reader(const char* path);
ChromData read_next_fasta(FastaReader* reader);
void free_fasta_reader(FastaReader**);
