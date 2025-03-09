#ifndef SECOND_LOOKUP_TABLE_H_INCLUDED
#define SECOND_LOOKUP_TABLE_H_INCLUDED

#include "split_database.h"

typedef struct
{
	int*  kmer_counts;
	int** kmer_positions;
	//int*  kmer_offsets;
} read_index;

read_index*
destroy_read_index(read_index* ridx, uint32_t index_count, int column_number);

read_index*
create_read_index(uint32_t index_count , int kmer_size, int positions_per_kmer);

void fill_read_index(read_index* index, uint32_t index_count, int kmer_size, int positions_per_kmer, const char* read, int read_size);

void earase_read_index(read_index* ridx, uint32_t index_count);

#endif // SECOND_LOOKUP_TABLE_H_INCLUDED
