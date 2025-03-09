#include "second_lookup_table.h"
#include "packed_db.h"

#include <cstdio>

read_index*
destroy_read_index(read_index* ridx, uint32_t index_count, int column_number)
{
	for (uint32_t i=0;i<index_count;++i){
        safe_free(ridx->kmer_positions[i]);
	}
	safe_free(ridx->kmer_counts);
	//safe_free(ridx->kmer_starts);
	safe_free(ridx->kmer_positions);
	safe_free(ridx);
	return NULL;
}

void earase_read_index(read_index* ridx, uint32_t index_count){
    for (uint32_t i=0; i<index_count; ++i){
        ridx->kmer_counts[i]=0;
    }
}

read_index*
create_read_index(uint32_t index_count, int kmer_size, int positions_per_kmer)
{
	DynamicTimer dtimer(__func__);
	//uint32_t index_count = 1 << (kmer_size * 2);
	//uint32_t leftnum = 34 - 2 * kmer_size;
	read_index* index = (read_index*)malloc(sizeof(read_index));
	safe_malloc(index->kmer_positions, int*, index_count);
	safe_malloc(index->kmer_counts, int, index_count);
	for (uint32_t i=0; i<index_count;++i){
        safe_malloc(index->kmer_positions[i], int, positions_per_kmer);
        assert(index->kmer_counts[i] == 0);
	}
    return index;
}

void
fill_read_index(read_index* index, uint32_t index_count, int kmer_size, int positions_per_kmer, const char* read, int read_size)
{
    uint32_t eit = 0;
    uint32_t leftnum = 34 - 2 * kmer_size;
    for (int j = 0; j < read_size; ++j)
    {
        eit = (eit << 2) | read[j];
        if (eit>index_count){
            int y=10;
        }
        assert(eit < index_count);
        if (j >= kmer_size - 1)
        {
            //if (index->kmer_starts[eit] && eit >= riti->min_key && eit <= riti->max_key)
            //{
            if (index->kmer_counts[eit]<positions_per_kmer){
                index->kmer_positions[eit][index->kmer_counts[eit]] = j + 1 - kmer_size;
                ++index->kmer_counts[eit];
            }
            //}
            eit <<= leftnum;
            eit >>= leftnum;
        }
    }
}
