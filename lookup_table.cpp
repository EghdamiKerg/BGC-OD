#include "lookup_table.h"
#include "packed_db.h"

#include <cstdio>

ref_index*
destroy_ref_index(ref_index* ridx)//need to change??????????
{
	safe_free(ridx->kmer_counts);
	safe_free(ridx->kmer_starts);
	safe_free(ridx->kmer_offsets);
	safe_free(ridx);
	return NULL;
}

typedef struct
{
	ref_index* ridx;
	uint32_t min_key;
	uint32_t max_key;
	volume_t* v;
	int kmer_size;
} ref_index_thread_info;

int quick_partition(float *a,int *b,int start_index,int end_index)
{
    float pivot=a[end_index];
    //P-index indicates the pivot value index

    int P_index=start_index;
    int i,t1; //t is temporary variable
    float t2;
    //Here we will check if array value is
    //less than pivot
    //then we will place it at left side
    //by swapping

    for(i=start_index;i<end_index;i++)
    {
    	if(a[i]<=pivot)
        {
            t2=a[i];
            a[i]=a[P_index];
            a[P_index]=t2;
            t1=b[i];
            b[i]=b[P_index];
            b[P_index]=t1;
            P_index++;
        }
     }
     //Now exchanging value of
     //pivot and P-index
      t2=a[end_index];
      a[end_index]=a[P_index];
      a[P_index]=t2;
      t1=b[end_index];
      b[end_index]=b[P_index];
      b[P_index]=t1;

     //at last returning the pivot value index
     return P_index;
 }
 void Quicksort(float *a, int *b, int start_index,int end_index)
 {
    if(start_index<end_index)
    {
         int P_index=quick_partition(a,b,start_index,end_index);
             Quicksort(a,b,start_index,P_index-1);
             Quicksort(a,b,P_index+1,end_index);
    }
}


void*
fill_ref_index_offsets_func(void* arg)
{
	ref_index_thread_info* riti = (ref_index_thread_info*)(arg);
	volume_t* v = riti->v;
	int num_reads = v->num_reads;
	ref_index* index = riti->ridx;
	int kmer_size = riti->kmer_size;
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	int i, j;
	for (i = 0; i < num_reads; ++i)
	//for (i = 0; i < 1000; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		//int read_end=read_start+read_size-1;

		uint32_t eit = 0;
		for (j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			eit = (eit << 2) | c;
			assert(eit < index_count);

			if (j >= kmer_size - 1)
			{
				if (index->kmer_starts[eit] && eit >= riti->min_key && eit <= riti->max_key)
				{

                    index->kmer_starts[eit][index->kmer_counts[eit]] = k + 1 - kmer_size;
                    ++index->kmer_counts[eit];
                }
                eit <<= leftnum;
				eit >>= leftnum;
			}
		}
	}
	/*for (uint32_t i =riti->min_key; i<=riti->max_key;++i){
        if(riti->ridx->kmer_starts[i]){
            Quicksort(index->kmer_gc_starts[i], index->kmer_starts[i],0,index->kmer_counts[i]-1);
        }
	}*/


	return NULL;
}

ref_index*
create_ref_index(volume_t* v, int kmer_size, int num_threads)
{
	DynamicTimer dtimer(__func__);
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	ref_index* index = (ref_index*)malloc(sizeof(ref_index));
	safe_calloc(index->kmer_counts, int, index_count);
	int num_reads = v->num_reads;
	for (uint32_t i = 0; i != index_count; ++i) {
        assert(index->kmer_counts [i] == 0);
	}
	for (int i = 0; i != num_reads; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		uint32_t eit = 0;
		for (int j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			assert(c>= 0 && c < 4);
			eit = (eit << 2) | c;
			if (j >= kmer_size - 1)
			{
				assert(eit < index_count);
                ++index->kmer_counts[eit];
                eit = eit << leftnum;
                eit = eit >> leftnum;
            }
		}
	}

	int num_kmers = 0;
	for (uint32_t i = 0; i != index_count; ++i)
	{
		if (index->kmer_counts[i] > 128) index->kmer_counts[i] = 0;
		num_kmers += index->kmer_counts[i];
	}

	printf("number of kmers: %d\n", num_kmers);
	safe_malloc(index->kmer_offsets, int, num_kmers);
	safe_malloc(index->kmer_starts, int*, index_count);

	if (v->curr < 10 * 1000000) num_threads = 1;
	int kmers_per_thread = (num_kmers + num_threads - 1) / num_threads;
	fprintf(stderr, "%d threads are used for filling offset lists.\n", num_threads);
	uint32_t hash_boundaries[2 * num_threads];
	uint32_t L = 0;
	num_kmers = 0;
	int kmer_cnt = 0;
	int tid = 0;
	for (uint32_t i = 0; i != index_count; ++i)
	{

        if (index->kmer_counts[i]!=0){
            index->kmer_starts[i] = index->kmer_offsets + num_kmers;
            num_kmers += index->kmer_counts[i];
            kmer_cnt += index->kmer_counts[i];
            index->kmer_counts[i] = 0;// why???
        }
        else{
            index->kmer_starts[i]=NULL;
            }

        if (kmer_cnt >= kmers_per_thread)
        {
            printf("thread %d: %d\t%d\n", tid, L, i);
            hash_boundaries[2 * tid] = L;
            hash_boundaries[2 * tid + 1] = i;
            ++tid;
            L = i + 1;
            kmer_cnt = 0;
        }
    }


	if (kmer_cnt)
	{
		printf("thread %d: %d\t%d\n", tid, L, index_count - 1);
		hash_boundaries[2 * tid] = L;
		hash_boundaries[2 * tid + 1] = index_count - 1;
	}

	ref_index_thread_info ritis[num_threads];
	for (int i = 0; i != num_threads; ++i)
	{
		ritis[i].ridx = index;
		ritis[i].min_key = hash_boundaries[2 * i];
		ritis[i].max_key = hash_boundaries[2 * i + 1];
		ritis[i].v = v;
		ritis[i].kmer_size = kmer_size;
	}

	pthread_t tids[num_threads];
	for (int j = 0; j < num_threads; ++j)
		pthread_create(tids + j, NULL, fill_ref_index_offsets_func, (void*)(ritis + j));
	for (int j = 0; j < num_threads; ++j)
		pthread_join(tids[j], NULL);

	return index;
}



/*void*
fill_ref_index_offsets_func(void* arg)
{
	ref_index_thread_info* riti = (ref_index_thread_info*)(arg);
	volume_t* v = riti->v;
	int num_reads = v->num_reads;
	ref_index* index = riti->ridx;
	int kmer_size = riti->kmer_size;
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	int max_window_size = 30;
	int i, j;
	for (i = 0; i < num_reads; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		int read_end = read_start + read_size - 1;
		int window_start_point = v->offset_list->offset_list[i].offset;
		int window_end_point = window_start_point + max_window_size/2 -1;
		int current_window_size = max_window_size/2;
		uint8_t window_start_char = PackedDB::get_char(v->data, window_start_point);
		int window_gc_count = 0;
		uint8_t temp_char;
		for(int j=window_start_point;j<=window_end_point;++j){
            temp_char = PackedDB::get_char(v->data, j);
            if (temp_char==1 || temp_char==2)
                ++window_gc_count;
		}
		uint32_t eit = 0;
		for (j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			eit = (eit << 2) | c;
			assert(eit < index_count);
			if (j >= kmer_size - 1)
			{
				if (eit >= riti->min_key && eit <= riti->max_key)
				{
					if (index->kmer_in_gc_rich_regions_starts[eit] && window_gc_count>(current_window_size/2)*1.2){
                        index->kmer_in_gc_rich_regions_starts[eit][index->kmer_in_gc_rich_regions_counts[eit]] = k + 1 - kmer_size;
                        ++index->kmer_in_gc_rich_regions_counts[eit];
                    }
                    else if (index->kmer_in_at_rich_regions_starts[eit] && window_gc_count<(current_window_size/2)*0.8){
                        index->kmer_in_at_rich_regions_starts[eit][index->kmer_in_at_rich_regions_counts[eit]] = k + 1 - kmer_size;
                        ++index->kmer_in_at_rich_regions_counts[eit];
                    }
                    else if (index->kmer_in_other_regions_starts[eit]{
                        index->kmer_in_other_regions_starts[eit][index->kmer_in_other_regions_counts[eit]] = k + 1 - kmer_size;
                        ++index->kmer_in_other_regions_counts[eit];
                    }
				}
                eit <<= leftnum;
				eit >>= leftnum;
			}
            if (window_start_point!=read_start && window_end_point!=read_end){
                ++window_end_point;
                temp_char = PackedDB::get_char(v->data, window_end_point);
                if (temp_char==1 || temp_char==2)
                    ++window_gc_count;
                if (window_start_char==1 || window_start_char==2)
                    --window_gc_count;
                window_start_char= PackedDB::get_char(v->data, ++window_start_point);
            }
            else if (window_start_point == read_start && current_window_size < max_window_size){
                ++window_end_point;
                ++current_window_size;
                temp_char = PackedDB::get_char(v->data, window_end_point);
                if (temp_char==1 || temp_char==2)
                    ++window_gc_count;
            }
            else if (window_end_point == read_end && current_window_size<=max_window_size/2){
                --current_window_size;
                if (window_start_char==1 || window_start_char==2)
                    --window_gc_count;
                window_start_char= PackedDB::get_char(v->data, ++window_start_point);
            }
		}
	}

	return NULL;
}

ref_index*
create_ref_index(volume_t* v, int kmer_size, int num_threads)
{
	DynamicTimer dtimer(__func__);
	uint32_t index_count = 1 << (kmer_size * 2);
	uint32_t leftnum = 34 - 2 * kmer_size;
	int max_window_size = 30;
	ref_index* index = (ref_index*)malloc(sizeof(ref_index));
	safe_calloc(index->kmer_in_at_rich_regions_counts, int, index_count);
	safe_calloc(index->kmer_in_gc_rich_regions_counts, int, index_count);
	safe_calloc(index->kmer_in_other_regions_counts, int, index_count);
	int num_reads = v->num_reads;
	for (uint32_t i = 0; i != index_count; ++i) {
        assert(index->kmer_in_at_rich_regions_counts [i] == 0);
        assert(index->kmer_in_gc_rich_regions_counts [i] == 0);
        assert(index->kmer_in_other_regions_counts [i] == 0);
	}
	for (int i = 0; i != num_reads; ++i)
	{
		int read_start = v->offset_list->offset_list[i].offset;
		int read_size = v->offset_list->offset_list[i].size;
		int read_end = read_start + read_size - 1;
		int window_start_point = v->offset_list->offset_list[i].offset;
		int window_end_point = window_start_point + max_window_size/2 -1;
		int current_window_size = max_window_size/2;
		uint8_t window_start_char = PackedDB::get_char(v->data, window_start_point);
		int window_gc_count = 0;
		uint8_t temp_char;
		for(int j=window_start_point;j<=window_end_point;++j){
            temp_char = PackedDB::get_char(v->data, j);
            if (temp_char==1 || temp_char==2)
                ++window_gc_count;
		}
		uint32_t eit = 0;
		for (int j = 0; j < read_size; ++j)
		{
			int k = read_start + j;
			uint8_t c = PackedDB::get_char(v->data, k);
			assert(c>= 0 && c < 4);
			eit = (eit << 2) | c;
			if (j >= kmer_size - 1)
			{
				assert(eit < index_count);
				if (window_gc_count>(current_window_size/2)*1.2){
                  ++index->kmer_in_gc_rich_regions_counts[eit];
				}
				else if (window_gc_count<(current_window_size/2)*0.8){
                    ++index->kmer_in_at_rich_regions_counts[eit];
				}
				else{
                    ++index->kmer_in_other_regions_counts[eit];
				}
				eit = eit << leftnum;
				eit = eit >> leftnum;
            }
            if (window_start_point!=read_start && window_end_point!=read_end){
                ++window_end_point;
                temp_char = PackedDB::get_char(v->data, window_end_point);
                if (temp_char==1 || temp_char==2)
                    ++window_gc_count;
                if (window_start_char==1 || window_start_char==2)
                    --window_gc_count;
                window_start_char= PackedDB::get_char(v->data, ++window_start_point);
            }
            else if (window_start_point == read_start && current_window_size < max_window_size){
                ++window_end_point;
                ++current_window_size;
                temp_char = PackedDB::get_char(v->data, window_end_point);
                if (temp_char==1 || temp_char==2)
                    ++window_gc_count;
            }
            else if (window_end_point == read_end && current_window_size<=max_window_size/2){
                --current_window_size;
                if (window_start_char==1 || window_start_char==2)
                    --window_gc_count;
                window_start_char= PackedDB::get_char(v->data, ++window_start_point);
            }
		}
	}


	int num_kmers = 0;
	for (uint32_t i = 0; i != index_count; ++i)
	{
		if (index->kmer_in_at_rich_regions_counts[i] > 128) index->kmer_in_at_rich_regions_counts[i] = 0;
		if (index->kmer_in_gc_rich_regions_counts[i] > 128) index->kmer_in_gc_rich_regions_counts[i] = 0;
		if (index->kmer_in_other_regions_counts[i] > 128) index->kmer_in_other_regions_counts[i] = 0;
		num_kmers += (index->kmer_in_at_rich_regions_counts[i] + index->kmer_in_gc_rich_regions_counts[i] + index->kmer_in_other_regions_counts[i]);
	}

	printf("number of kmers: %d\n", num_kmers);
	safe_malloc(index->kmer_offsets, int, num_kmers);
	safe_malloc(index->kmer_in_at_rich_regions_starts, int*, index_count);
	safe_malloc(index->kmer_in_gc_rich_regions_starts, int*, index_count);
	safe_malloc(index->kmer_in_other_regions_starts, int*, index_count);

	if (v->curr < 10 * 1000000) num_threads = 1;
	int kmers_per_thread = (num_kmers + num_threads - 1) / num_threads;
	fprintf(stderr, "%d threads are used for filling offset lists.\n", num_threads);
	uint32_t hash_boundaries[2 * num_threads];
	uint32_t L = 0;
	num_kmers = 0;
	int kmer_cnt = 0;
	int tid = 0;
	for (uint32_t i = 0; i != index_count; ++i)
	{

        if (index->kmer_in_at_rich_regions_counts!=0){
            index->kmer_in_at_rich_regions_starts[i] = index->kmer_offsets + num_kmers;
            num_kmers += index->kmer_in_at_rich_regions_counts[i];
            kmer_cnt += index->kmer_in_at_rich_regions_counts[i];
            //index->kmer_counts[i] = 0;// why???
            index->kmer_in_at_rich_regions_counts[i] = 0;
        }
        else{
            index->kmer_in_at_rich_regions_starts=NULL;
        }
        if (index->kmer_in_gc_rich_regions_counts!=0){
            index->kmer_in_gc_rich_regions_starts[i] = index->kmer_offsets + num_kmers;
            num_kmers += index->kmer_in_gc_rich_regions_counts[i];
            kmer_cnt += index->kmer_in_gc_rich_regions_counts[i];
            //index->kmer_counts[i] = 0;// why???
            index->kmer_in_gc_rich_regions_counts[i] = 0;
        }
        else{
            index->kmer_in_gc_rich_regions_starts=NULL;
        }
        if (index->kmer_in_other_regions_counts!=0){
            index->kmer_in_other_regions_starts[i] = index->kmer_offsets + num_kmers;
            num_kmers += index->kmer_in_other_regions_counts[i];
            kmer_cnt += index->kmer_in_other_regions_counts[i];
            //index->kmer_counts[i] = 0;// why???
            index->kmer_in_other_regions_counts[i] = 0;
        }
        else{
            index->kmer_in_other_regions_starts=NULL;
        }

        if (kmer_cnt >= kmers_per_thread)
        {
            printf("thread %d: %d\t%d\n", tid, L, i);
            hash_boundaries[2 * tid] = L;
            hash_boundaries[2 * tid + 1] = i;
            ++tid;
            L = i + 1;
            kmer_cnt = 0;
        }
    }


	if (kmer_cnt)
	{
		printf("thread %d: %d\t%d\n", tid, L, index_count - 1);
		hash_boundaries[2 * tid] = L;
		hash_boundaries[2 * tid + 1] = index_count - 1;
	}

	ref_index_thread_info ritis[num_threads];
	for (int i = 0; i != num_threads; ++i)
	{
		ritis[i].ridx = index;
		ritis[i].min_key = hash_boundaries[2 * i];
		ritis[i].max_key = hash_boundaries[2 * i + 1];
		ritis[i].v = v;
		ritis[i].kmer_size = kmer_size;
	}

	pthread_t tids[num_threads];
	for (int j = 0; j < num_threads; ++j)
		pthread_create(tids + j, NULL, fill_ref_index_offsets_func, (void*)(ritis + j));
	for (int j = 0; j < num_threads; ++j)
		pthread_join(tids[j], NULL);

	return index;
}*/
