#include "split_database.h"
#include "pw_options.h"
#include "diff_gapalign.h"
#include "xdrop_gapalign.h"
#include "packed_db.h"
#include "lookup_table.h"
#include "second_lookup_table.h"
#include "pw_impl.h"

#include <algorithm>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <tgmath.h>

#define MSS MAX_SEQ_SIZE

static int MAXC = 100;
static int output_gapped_start_point = 1;
static int kmer_size = 13;
static const double ddfs_cutoff_pacbio = 0.25;
static const double ddfs_cutoff_nanopore = 0.25;
static double ddfs_cutoff = ddfs_cutoff_pacbio;
static int min_align_size = 0;
static int min_kmer_match = 0;
static int min_kmer_dist = 0;
int all_anchors=0;
int ignored_anchors=0;

using namespace std;

PWThreadData::PWThreadData(options_t* opt, volume_t* ref, volume_t* rd, ref_index* idx, std::ostream* o)
	: options(opt), used_thread_id(0), reference(ref), reads(rd), ridx(idx), out(o), m4_results(NULL), ec_results(NULL), next_processed_id(0), self_detection(false)
{
	pthread_mutex_init(&id_lock, NULL);
	if (options->task == TASK_SEED)
	{
		safe_malloc(ec_results, ExtensionCandidate*, options->num_threads);
		for (int i = 0; i < options->num_threads; ++i)
			safe_malloc(ec_results[i], ExtensionCandidate, kResultListSize);
	}
	else if (options->task == TASK_ALN)
	{
		safe_malloc(m4_results, M4Record*, options->num_threads);
		for (int i = 0; i < options->num_threads; ++i)
			safe_malloc(m4_results[i], M4Record, kResultListSize);
	}
	else
	{
		LOG(stderr, "Task must be either %d or %d, not %d!", TASK_SEED, TASK_ALN, options->task);
		abort();
	}
	pthread_mutex_init(&result_write_lock, NULL);
	pthread_mutex_init(&read_retrieve_lock, NULL);
}

PWThreadData::~PWThreadData()
{
	if (ec_results)
	{
		for (int i = 0; i < options->num_threads; ++i) safe_free(ec_results[i]);
		safe_free(ec_results);
	}
	if (m4_results)
	{
		for (int i = 0; i < options->num_threads; ++i) safe_free(m4_results[i]);
		safe_free(m4_results);
	}
}

void
reverse_complement(char* dst, const char* src, const int size)
{
	const uint8_t* rc_table = get_dna_complement_table();
	int i;
	for (i = 0; i < size; ++i)
	{
		uint8_t c = src[i];
		assert(c >= 0 && c < 4);
		c = rc_table[c];
		dst[size - 1 - i] = c;
	}
}

int
extract_kmers(const char* s, const int sstart, const int ssize, const uint8_t* mataData, int* kmer_ids, short* kmer_gc_percents)
{
	uint32_t max_id = 1 << (2 * kmer_size);
	//int num_kmers = ssize / kmer_size;
	int read_end = ssize - 1;
	int i=0, j;
    uint32_t eit = 0;
    //uint32_t leftnum = 34 - 2 * kmer_size;
    int p=kmer_size-1;
    for (j = 0; j < ssize; ++j)
    {
        eit = (eit << 2) | s[j];
			//assert(eit < index_count);

        if (j==p)
        {
            assert(eit < max_id);
            kmer_ids[i] = eit;
            kmer_gc_percents[i]=mataData[sstart+j-kmer_size+1];
            eit= 0;
            ++i;
            p=(i+1)*kmer_size-1;
            //p+=(kmer_size/2);
        }
        //if (j>=kmer_size-1){
            //eit <<= leftnum;
            //eit >>= leftnum;
        //}
    }
	return i-1;
}

void
extract_randomized_kmers(volume_t* v, const int start_pos, int* kmer_ids, uint32_t index_count)
{
	//int max_id = 1 << (2 * kmer_size);
	//int num_kmers = ssize / kmer_size;
	int i,j;
	bool bad_pos=false;
	//int end_pos = ssize - 1;
    uint32_t eit = 0;
    srand((unsigned) time(NULL));
    int random_pos;
    for (j = 0; j < NUM_RANDOMIZED_KMERS; ++j)
    {
        random_pos = start_pos + (rand() % (BLOCK_SIZE-READ_KMER_SIZE));
        for (i=random_pos; i<random_pos+READ_KMER_SIZE; ++i){
            uint8_t c = PackedDB::get_char(v->data, i);
            if (c>3 || c<0){
                bad_pos=true;
                break;
            }
			eit = (eit << 2) | c;
			assert(eit < index_count);//***************************
        }
        if(bad_pos){
            j--;
            eit=0;
            bad_pos=false;
            continue;
        }
        kmer_ids[j] = eit;
        eit= 0;
    }
}

/*int find_relevant_kmers(int kmer_id ,float kmer_gc_percent, ref_index* ridx, int** seed_arr){
    float start_threshold = kmer_gc_percent - 20;
    float end_threshold = kmer_gc_percent + 20;
    int num_kmers_all = ridx->kmer_counts[kmer_id];
    bool sw=false;
    int i,j;
    for(i=0; i<num_kmers_all; ++i){
        if (ridx->kmer_gc_starts[kmer_id][i] >= start_threshold){
            sw=true;
            *seed_arr= ridx->kmer_starts[kmer_id]+i;
            break;
        }
    }
    if (sw){
        for(j=i+1; j<num_kmers_all;++j){
            if(ridx->kmer_gc_starts[kmer_id][i]>end_threshold){
                return j-i;
            }
        }
        return j-i;
    }
    else{
        return 0;
    }
    //*seed_arr= ridx->kmer_starts[kmer_id];
    //return num_kmers_all;

}*/

SeedingBK::SeedingBK(const int ref_size)
{
	const int num_rows = MAX_SEQ_SIZE / BLOCK_SIZE + 1;
	const int num_columns = ref_size / BLOCK_SIZE + 1;
	safe_malloc(marked_blocks_list, marked_blocks, (num_columns + num_rows)*3);
	safe_malloc(database, int*, num_rows);
	for (int i=0;i<num_rows;i++){
        safe_malloc(database[i],int,num_columns);
	}
	safe_malloc(kmer_ids, int, MAX_SEQ_SIZE);
	safe_malloc(kmer_gc_percents, short, MAX_SEQ_SIZE);
	for (int i = 0; i < num_rows; ++i)
	{
	    for (int j=0; j<num_columns; ++j){
            //database[i][j].score = 0;
            database[i][j] = -1;
	    }
	}
}

SeedingBK::~SeedingBK()
{
	safe_free(marked_blocks_list);
	const int num_rows = MAX_SEQ_SIZE / BLOCK_SIZE + 1;
	for(int i=0;i<num_rows;i++){
        safe_free(database[i]);
	}
	safe_free(database);
	safe_free(kmer_ids);
	safe_free(kmer_gc_percents);
}



int
seeding(const char* read, const int read_start, const int read_size, const uint8_t* metaData, ref_index* ridx, SeedingBK* sbk)
{
	int* kmer_ids = sbk->kmer_ids;
	short* kmer_gc_percents = sbk->kmer_gc_percents;
	marked_blocks* marked_blocks_list = sbk->marked_blocks_list;
	marked_blocks* marked_blocks_list_spr = sbk->marked_blocks_list;
	//short* index_score = sbk->index_score;
	//short* index_ss = index_score;
	int** blocks = sbk->database;
    int* current_row;
    int* current_block;
	int num_kmers = extract_kmers(read, read_start, read_size, metaData, kmer_ids, kmer_gc_percents);
	int km;
	int used_segs = 0;
	for (km = 0; km < num_kmers; ++km)
	{
		int num_seeds = ridx->kmer_counts[kmer_ids[km]];
		int* seed_arr = ridx->kmer_starts[kmer_ids[km]];
		//int* seed_arr= NULL;
		//int num_seeds = find_relevant_kmers(kmer_ids[km] ,kmer_gc_percents[km], ridx, &seed_arr);
		if (num_seeds >= 0){
            int sid;
            int block_row = km*(kmer_size)/BLOCK_SIZE;
            current_row = *(blocks + block_row);
            for (sid = 0; sid < num_seeds; ++sid)
            {
                all_anchors++;
                if (metaData[seed_arr[sid]]<kmer_gc_percents[km]-GC_Dif || metaData[seed_arr[sid]]>kmer_gc_percents[km]+GC_Dif){
                    ignored_anchors++;
                    continue;
                }
		    //************************May change*****************************
                //if(seed_arr[sid]>=read_start && seed_arr[sid]<read_start+read_size){
                    //continue;
                //}
            //**************************************************************
                int block_column = seed_arr[sid] / BLOCK_SIZE;
                //int block_off = seed_arr[sid] % BLOCK_SIZE;
                current_block = current_row + block_column;
                //if (spr->score == 0 || spr->seednum < km + 1)
                //{
				//++current_block->score;
				//if (loc <= SM) { spr->loczhi[loc - 1] = seg_off; spr->seedno[loc - 1] = km + 1; }
				//else insert_loc(spr, seg_off, km + 1, BC);
				//int s_k;
				//if (seg_id > 0) s_k = spr->score + (spr - 1)->score;
				//else s_k = spr->score;
				//if (endnum < s_k) endnum = s_k;
				if (*current_block == -1)
				//if (spr->score==2)
				{
					//*(index_spr++) = seg_id;
					//*(index_ss++) = s_k;
					marked_blocks_list_spr->row_index=block_row;
					marked_blocks_list_spr->column_index=block_column;
					marked_blocks_list_spr->block_score=1;
					marked_blocks_list_spr++;
					//used_segs++;
					(*current_block) = used_segs++;
				}
				else marked_blocks_list[*current_block].block_score++;
			//}
			//spr->seednum = km + 1;
            }
		}
	}
	return used_segs;
}

int
self_seeding(const char* read, const int read_start, const int read_size, const uint8_t* metaData, ref_index* ridx, SeedingBK* sbk)
{
	int* kmer_ids = sbk->kmer_ids;
	short* kmer_gc_percents = sbk->kmer_gc_percents;
	marked_blocks* marked_blocks_list = sbk->marked_blocks_list;
	marked_blocks* marked_blocks_list_spr = sbk->marked_blocks_list;
	//short* index_score = sbk->index_score;
	//short* index_ss = index_score;
	int** blocks = sbk->database;
    int* current_row;
    int* current_block;
	int num_kmers = extract_kmers(read, read_start, read_size, metaData, kmer_ids, kmer_gc_percents);
	int km;
	int used_segs = 0;
	for (km = 0; km < num_kmers; ++km)
	{
		int num_seeds = ridx->kmer_counts[kmer_ids[km]];
		int* seed_arr = ridx->kmer_starts[kmer_ids[km]];
		//int* seed_arr= NULL;
		//int num_seeds = find_relevant_kmers(kmer_ids[km] ,kmer_gc_percents[km], ridx, &seed_arr);
		if (num_seeds > 0){
            int sid;
            int block_row = km*(kmer_size)/BLOCK_SIZE;
            current_row = *(blocks + block_row);
            for (sid = 0; sid < num_seeds; ++sid)
            {
		    //************************May change*****************************
                if(seed_arr[sid] >= read_start){
                    continue;
                }
            //**************************************************************
                all_anchors++;
                if (metaData[seed_arr[sid]]<kmer_gc_percents[km]-GC_Dif || metaData[seed_arr[sid]]>kmer_gc_percents[km]+GC_Dif){
                    ignored_anchors++;
                    continue;
                }

                int block_column = seed_arr[sid] / BLOCK_SIZE;
                //int block_off = seed_arr[sid] % BLOCK_SIZE;
                current_block = current_row + block_column;
                //if (spr->score == 0 || spr->seednum < km + 1)
                //{
				//++current_block->score;
				//if (loc <= SM) { spr->loczhi[loc - 1] = seg_off; spr->seedno[loc - 1] = km + 1; }
				//else insert_loc(spr, seg_off, km + 1, BC);
				//int s_k;
				//if (seg_id > 0) s_k = spr->score + (spr - 1)->score;
				//else s_k = spr->score;
				//if (endnum < s_k) endnum = s_k;
				if (*current_block == -1)
				//if (spr->score==2)
				{
					//*(index_spr++) = seg_id;
					//*(index_ss++) = s_k;
					marked_blocks_list_spr->row_index=block_row;
					marked_blocks_list_spr->column_index=block_column;
					marked_blocks_list_spr->block_score=1;
					marked_blocks_list_spr++;
					//used_segs++;
					(*current_block) = used_segs++;
				}
				else marked_blocks_list[*current_block].block_score++;
			//}
			//spr->seednum = km + 1;
            }
		}
	}
	return used_segs;
}

bool compare_marked_blocks (marked_blocks i, marked_blocks j){
    if (i.column_index<j.column_index)
        return true;
    if(i.column_index==j.column_index && i.row_index<j.row_index)
        return true;
    return false;
}

bool compare_temp_arr_elements(temp_arr_element i, temp_arr_element j){
    return i.score>j.score;
}

void chaining (marked_blocks* marked_blocks_list, temp_arr_element* temp_arr, int num_segs, int read_size){
    marked_blocks* marked_blocks_list_spr1 = marked_blocks_list;
    marked_blocks* marked_blocks_list_spr2;
    temp_arr_element* temp_arr_spr = temp_arr;
    int i, j;
    int temp_score;
    int temp_gap_penalty;
    int num_read_blocks = DIV_BLOCK_SIZE(read_size)+1;
    /*for (i=0; i<num_segs; ++i, ++marked_blocks_list_spr1, ++temp_arr_spr){
        marked_blocks_list_spr1->visited=false;
        temp_arr_spr->score=marked_blocks_list_spr1->chaining_score=1;
        marked_blocks_list_spr1->prev=-1;
        temp_arr_spr->index=i;
        //temp_arr_spr->score=2;
        for (j=0 , marked_blocks_list_spr2 = marked_blocks_list_spr1-1; j<20 && j<i; ++j , --marked_blocks_list_spr2){
            if (marked_blocks_list_spr2->column_index - marked_blocks_list_spr1->column_index > num_read_blocks)
                break;
            if(marked_blocks_list_spr1->row_index<marked_blocks_list_spr2->row_index){
                continue;
            }
            temp_gap_penalty = max(marked_blocks_list_spr1->row_index - marked_blocks_list_spr2->row_index,
                                   marked_blocks_list_spr1->column_index - marked_blocks_list_spr2->column_index);
            if(temp_gap_penalty>MAX_GAP)
                continue;
            temp_gap_penalty = abs((marked_blocks_list_spr1->row_index - marked_blocks_list_spr2->row_index)-
                                   (marked_blocks_list_spr1->column_index - marked_blocks_list_spr2->column_index));
            temp_gap_penalty = 0.1 * temp_gap_penalty;// + 0.5 * log2 (temp_gap_penalty);
            temp_score=marked_blocks_list_spr2->chaining_score + 1;//min(marked_blocks_list_spr1->row_index - marked_blocks_list_spr2->row_index,marked_blocks_list_spr1->column_index - marked_blocks_list_spr2->column_index)+marked_blocks_list_spr1->block_score;
            temp_score-=temp_gap_penalty;
            if (temp_score > marked_blocks_list_spr1->chaining_score){
                marked_blocks_list_spr1->chaining_score=temp_score;
                marked_blocks_list_spr1->prev=i-j-1;
                temp_arr_spr->score=temp_score;
            }
        }
    }*/
    for (i=0; i<num_segs; ++i, ++marked_blocks_list_spr1, ++temp_arr_spr){
        marked_blocks_list_spr1->visited=false;
        marked_blocks_list_spr1->chaining_score=1;//was changed
        marked_blocks_list_spr1->prev=-1;
        temp_arr_spr->index=i;
        //temp_arr_spr->score=BLOCK_SIZE;
        temp_arr_spr->score=1;//was changed
        for (j=0 , marked_blocks_list_spr2 = marked_blocks_list_spr1-1; j<20 && j<i; ++j , --marked_blocks_list_spr2){
            if (marked_blocks_list_spr2->column_index - marked_blocks_list_spr1->column_index > num_read_blocks)
                break;
            if(marked_blocks_list_spr1->row_index<marked_blocks_list_spr2->row_index){
                continue;
            }
            temp_gap_penalty = max(marked_blocks_list_spr1->row_index - marked_blocks_list_spr2->row_index,
                                   marked_blocks_list_spr1->column_index - marked_blocks_list_spr2->column_index);
            if(temp_gap_penalty>MAX_GAP)
                continue;
            temp_gap_penalty = abs((marked_blocks_list_spr1->row_index - marked_blocks_list_spr2->row_index)-
                                   (marked_blocks_list_spr1->column_index - marked_blocks_list_spr2->column_index));
            temp_score=marked_blocks_list_spr2->chaining_score + min(marked_blocks_list_spr1->row_index - marked_blocks_list_spr2->row_index,
                                   marked_blocks_list_spr1->column_index - marked_blocks_list_spr2->column_index);
            temp_score-=(temp_gap_penalty*0.1);
            //temp_score=MUL_BLOCK_SIZE(temp_score);//was changed
            if (temp_score > marked_blocks_list_spr1->chaining_score){
                marked_blocks_list_spr1->chaining_score=temp_score;
                marked_blocks_list_spr1->prev=i-j-1;
                temp_arr_spr->score=temp_score;
            }
        }
    }
}

int
get_candidates(volume_t* ref,
			   SeedingBK* sbk,
			   read_index* read_index,
			   const int32_t read_index_count,
			   int num_segs,
			   const char* read,
			   const int read_id,
			   const int read_size,
			   const char chain,
			   candidate_save* candidates,
			   int candidatenum)
{
	marked_blocks* marked_blocks_list = sbk->marked_blocks_list;
	marked_blocks* marked_blocks_list_spr = marked_blocks_list;
	marked_blocks* marked_blocks_list_end = marked_blocks_list+num_segs;

	int** database = sbk->database;
	/*int current_row;
	int current_column;
	Back_List* current_block;
	int* new_seeds_positions;
	int* new_kmer_ids;
	int new_seeds_count;
	int index;
	safe_malloc(new_kmer_ids, int, NUM_RANDOMIZED_KMERS);*/

	int index;
	//float min_chain_score=2*BLOCK_SIZE;
	float min_chain_score=2;//was changed
	candidate_save *candidate_loc = candidates, candidate_temp;
	int location_loc[4];
	int i, j, k;
	bool sw=false;
	/*for (i = 0; i < num_segs; ++i, ++marked_blocks_list_spr){
		current_row = marked_blocks_list_spr->row_index;
		current_column = marked_blocks_list_spr->column_index;
		current_block = *(database + current_row) + current_column;
		//if (current_block->score < min_kmer_match){//*******************************may change
                //continue;
		//}
		if (database[current_row][current_column+1].index!=-1 || database[current_row+1][current_column+1].index!=-1 || database[current_row+2][current_column+1].index!=-1){
            continue;
		}
		if (sw==false){
            sw=true;
            fill_read_index(read_index,read_index_count,READ_KMER_SIZE,POSITIONS_PER_KMER, read, read_size);
		}
		int ref_start_pos = MUL_BLOCK_SIZE(current_column);
        extract_randomized_kmers(ref, ref_start_pos, new_kmer_ids, read_index_count);
        for (j=0 ; j<NUM_RANDOMIZED_KMERS; ++j){
            new_seeds_positions = read_index->kmer_positions[new_kmer_ids[j]];
            new_seeds_count = read_index->kmer_counts[new_kmer_ids[j]];
            for (k=0; k<new_seeds_count; ++k){
                int new_seed_pos_row = DIV_BLOCK_SIZE(*(new_seeds_positions+k));
                if (new_seed_pos_row==current_row || new_seed_pos_row==current_row+1 || new_seed_pos_row==current_row+2){
                    ++database[new_seed_pos_row][current_column+1].score;
                    if(database[new_seed_pos_row][current_column+1].score>=3 && database[new_seed_pos_row][current_column+1].index==-1){
                        marked_blocks_list_end->column_index=current_column+1;
                        marked_blocks_list_end->row_index=new_seed_pos_row;
                        marked_blocks_list_end++;
                        database[new_seed_pos_row][current_column+1].index=num_segs++;
                    }
                }
            }
        }
	}*/
    int last_sid=-1;
    sort(marked_blocks_list, marked_blocks_list+num_segs,compare_marked_blocks);
    temp_arr_element* temp_arr;
    safe_malloc(temp_arr,temp_arr_element,num_segs);
    if(read_id==34)
        int y=0;
    chaining (marked_blocks_list, temp_arr, num_segs, read_size);
    sort(temp_arr, temp_arr+num_segs,compare_temp_arr_elements);
    i=0;
    float chain_score;
    temp_arr_element* temp_arr_spr = temp_arr;
    index = temp_arr_spr->index;
    chain_score = marked_blocks_list[index].chaining_score;
    bool prev_visited_chain;
    if (read_id==47)
        int y=0;
    while (i<num_segs && chain_score>min_chain_score){


            //index = temp_arr_spr->index;
            //chain_score = marked_blocks_list[index].chaining_score;
            //if (chain_score<5){//??????????????????
                //break;
            //}
            //while (!marked_blocks_list[index].visited && i<num_segs){
                    //marked_blocks_list[index].visited=true;
                    //++temp_arr_spr;
                    //index = temp_arr_spr->index;
                    //++i;
            //}
            //if (i>num_segs){
                //int y=0;
            //}
            prev_visited_chain=false;
            candidate_temp.score=marked_blocks_list[index].chaining_score;
            //marked_blocks_list[index].visited=true;
            location_loc[2] = MUL_BLOCK_SIZE(marked_blocks_list[index].column_index + 1);
            location_loc[3] = MUL_BLOCK_SIZE(marked_blocks_list[index].row_index + 1);
            while (marked_blocks_list[index].prev!=-1){
                    if(marked_blocks_list[index].visited){
                        prev_visited_chain=true;
                        break;
                    }
                    marked_blocks_list[index].visited=true;
                    index=marked_blocks_list[index].prev;
                    ++i;
            }
            if(prev_visited_chain){
                do{
                    temp_arr_spr++;
                    index=temp_arr_spr->index;
                    chain_score=marked_blocks_list[index].chaining_score;
                }while(i<num_segs && chain_score>min_chain_score && marked_blocks_list[index].visited);
                continue;
            }
            marked_blocks_list[index].visited=true;
            ++i;
            location_loc[0] = MUL_BLOCK_SIZE(marked_blocks_list[index].column_index);
            location_loc[1] = MUL_BLOCK_SIZE(marked_blocks_list[index].row_index);
            if (location_loc[2]-location_loc[0]<2000){
                //break;
                do{
                    temp_arr_spr++;
                    index=temp_arr_spr->index;
                    chain_score=marked_blocks_list[index].chaining_score;
                }while(i<num_segs && chain_score>min_chain_score && marked_blocks_list[index].visited);
                continue;
            }
			candidate_temp.chain = chain;
			//int loc_list = location_loc[0];
			int sid = get_read_id_from_offset_list(ref->offset_list, (location_loc[0]+location_loc[2])/2);
			if (sid == last_sid) {
                do{
                    temp_arr_spr++;
                    index=temp_arr_spr->index;
                    chain_score=marked_blocks_list[index].chaining_score;
                }while(i<num_segs && chain_score>min_chain_score && marked_blocks_list[index].visited);
                continue;/////*****************why
            }
			int sstart = ref->offset_list->offset_list[sid].offset;
			int ssize = ref->offset_list->offset_list[sid].size;
			int send = sstart + ssize + 1;/////////////why +1
			sid += ref->start_read_id;//////important
			candidate_temp.readno = sid;
            candidate_temp.readstart = sstart;
            //***********************************************
            if (location_loc[0] < sstart)
                location_loc[0]=sstart;
            //***********************************************
            int left_length1 = location_loc[0] - sstart + kmer_size - 1;
            int right_length1 = send - location_loc[0];
            int left_length2 = location_loc[1] + kmer_size - 1;
            int right_length2 = read_size - location_loc[1];
            int num1 = (left_length1 > left_length2) ? left_length2 : left_length1;////////////////what
            int num2 = (right_length1 > right_length2) ? right_length2 : right_length1;
            if (num1 + num2 < min_kmer_dist) {
                do{
                    temp_arr_spr++;
                    index=temp_arr_spr->index;
                    chain_score=marked_blocks_list[index].chaining_score;
                }while(i<num_segs && chain_score>min_chain_score && marked_blocks_list[index].visited);
                continue;/////*****************why
            }
            candidate_temp.loc1=location_loc[0] - sstart;
            if (candidate_temp.loc1<0)
                int y=0;
            candidate_temp.num1=num1;
            candidate_temp.loc2=location_loc[1];
            if (candidate_temp.loc2<0)
                int y=0;
            candidate_temp.num2=num2;
            candidate_temp.left1=left_length1;
            candidate_temp.left2=left_length2;
            candidate_temp.right1=right_length1;
            candidate_temp.right2=right_length2;
            candidate_temp.score=chain_score;
            //insert canidate position or delete this position
            int low=0;
            int high=candidatenum-1;
            int mid;
            while(low<=high)
            {
                mid=(low+high)/2;
                if(mid>=candidatenum||candidate_loc[mid].score<candidate_temp.score) high=mid-1;
                else low=mid+1;
            }
            if(candidatenum<MAXC)
                for(j=candidatenum-1; j>high; j--)
                    candidate_loc[j+1]=candidate_loc[j];
            else
                for(j=candidatenum-2; j>high; j--)
                    candidate_loc[j+1]=candidate_loc[j];
            if(high+1<MAXC){
                candidate_loc[high+1]=candidate_temp;
                last_sid=sid;
            }
            if(candidatenum<MAXC)
                candidatenum++;
            else candidatenum=MAXC;

            do{
                temp_arr_spr++;
                index=temp_arr_spr->index;
                chain_score=marked_blocks_list[index].chaining_score;
            }while(i<num_segs && chain_score>min_chain_score && marked_blocks_list[index].visited);
    }
	for(i=0,marked_blocks_list_spr=marked_blocks_list; i<num_segs; i++,marked_blocks_list_spr++)
	{
	    int r=marked_blocks_list_spr->row_index;
	    int c=marked_blocks_list_spr->column_index;
		//database[r][c].score=0;
		database[r][c]=-1;
	}
	//safe_free(new_kmer_ids);
	safe_free(temp_arr);
	if (sw)
        earase_read_index(read_index, read_index_count);
	return candidatenum;
}

void
fill_m4record(GapAligner* aligner, const int qid, const int sid,
			  const char qchain, int qsize, int ssize,
			  int qstart, int sstart, int vscore, M4Record* m)
{
	if (qchain == 'F')
	{
		m->qid = sid;
		m->sid = qid;
		m->ident = aligner->calc_ident();
		m->vscore = vscore;
		m->qdir = 0;
		m->qoff = aligner->target_start();
		m->qend = aligner->target_end();
		m->qsize = ssize;
		m->sdir = 0;
		m->soff = aligner->query_start();
		m->send = aligner->query_end();
		m->ssize = qsize;
		m->qext = sstart;
		m->sext = qstart;
	}
	else
	{
		m->qid = sid;
		m->sid = qid;
		m->ident = aligner->calc_ident();
		m->vscore = vscore;
		m->qdir = 0;
		m->qoff = aligner->target_start();
		m->qend = aligner->target_end();
		m->qsize = ssize;
		m->sdir = 1;
		m->soff = qsize - aligner->query_end();
		m->send = qsize - aligner->query_start();
		m->ssize = qsize;
		m->qext = sstart;
		m->sext = qsize - 1 - qstart;
	}
}


void output_m4record(ostream& out, const M4Record& m4)
{
	const char sep = '\t';

	out << m4qid(m4)    << sep
	    << m4sid(m4)    << sep
	    << m4ident(m4) << sep
	    << m4vscore(m4)   << sep
	    << m4qdir(m4)   << sep
	    << m4qoff(m4)   << sep
	    << m4qend(m4)   << sep
	    << m4qsize(m4)  << sep
	    << m4sdir(m4)   << sep
	    << m4soff(m4)   << sep
	    << m4send(m4)   << sep
	    << m4ssize(m4);

	if (output_gapped_start_point)
		out << sep
			<< m4qext(m4)	<< sep
			<< m4sext(m4);
	out << "\n";
}

void
print_m4record_list(ostream* out, M4Record* m4_list, int num_m4)
{
	for (int i = 0; i < num_m4; ++i) output_m4record(*out, m4_list[i]);
}

struct CmpM4RecordByQidAndOvlpSize
{
	bool operator()(const M4Record& a, const M4Record& b)
	{
		if (m4qid(a) != m4qid(b)) return m4qid(a) < m4qid(b);
		int o1 = M4RecordOverlapSize(a);
		int o2 = M4RecordOverlapSize(b);
		return o1 > o2;
	}
};

void
check_records_containment(M4Record* m4v, int s, int e, int* valid)
{
	const int soft = 100;

	for (int i = s; i < e; ++i)
	{
		if (!valid[i]) continue;
		int qb1 = m4qoff(m4v[i]);
		int qe1 = m4qend(m4v[i]);
		int sb1 = m4soff(m4v[i]);
		int se1 = m4send(m4v[i]);
		for (int j = i + 1; j < e; ++j)
		{
			if (!valid[j]) continue;
			if (m4sdir(m4v[i]) != m4sdir(m4v[j])) continue;
			int qb2 = m4qoff(m4v[j]);
			int qe2 = m4qend(m4v[j]);
			int sb2 = m4soff(m4v[j]);
			int se2 = m4send(m4v[j]);

			if (qb2 + soft >= qb1 && qe2 - soft <= qe1 && sb2 + soft >= sb1 && se2 - soft <= se1) valid[j] = 0;
		}
	}
}

void
append_m4v(M4Record* glist, int* glist_size,
		   M4Record* llist, int* llist_size,
		   ostream* out, pthread_mutex_t* results_write_lock)
{
	sort(llist, llist + *llist_size, CmpM4RecordByQidAndOvlpSize());
	int i = 0, j;
	int valid[*llist_size];
	fill(valid, valid + *llist_size, 1);
	while (i < *llist_size)
	{
		idx_t qid = m4qid(llist[i]);
		j = i + 1;
		while (j < *llist_size && m4qid(llist[j]) == qid) ++j;
		if (j - i > 1) check_records_containment(llist, i, j, valid);
		i = j;
	}

	if ((*glist_size) + (*llist_size) > PWThreadData::kResultListSize)
	{
		pthread_mutex_lock(results_write_lock);
		print_m4record_list(out, glist, *glist_size);
		*glist_size = 0;
		pthread_mutex_unlock(results_write_lock);
	}

	for (i = 0; i < *llist_size; ++i)
		if (valid[i])
		{
			glist[*glist_size] = llist[i];
			++(*glist_size);
		}

	*llist_size = 0;
}

inline void
get_next_chunk_reads(PWThreadData* data, int& Lid, int& Rid)
{
	pthread_mutex_lock(&data->read_retrieve_lock);
	Lid = data->next_processed_id;
	Rid = Lid + CHUNK_SIZE;
	if (Rid > data->reads->num_reads) Rid = data->reads->num_reads;
	data->next_processed_id += CHUNK_SIZE;
	pthread_mutex_unlock(&data->read_retrieve_lock);
}

void
pairwise_mapping(PWThreadData* data, int tid)
{
}

/*void
pairwise_mapping(PWThreadData* data, int tid)
{
	char *read, *read1, *read2, *subject;
	safe_malloc(read1, char, MSS);
	safe_malloc(read2, char, MSS);
	safe_malloc(subject, char, MSS);
	SeedingBK* sbk = new SeedingBK(data->reference->curr);
	candidate_save candidates[MAXC];
	int num_candidates = 0;
	M4Record* m4_list = data->m4_results[tid];
	int m4_list_size = 0;
	M4Record* m4v = new M4Record[MAXC];
	int num_m4 = 0;
	GapAligner* aligner = NULL;
	if (data->options->tech == TECH_PACBIO) {
		aligner = new DiffAligner(0);
	} else if (data->options->tech == TECH_NANOPORE) {
		aligner = new XdropAligner(0);
	} else {
		ERROR("TECH must be either %d or %d", TECH_PACBIO, TECH_NANOPORE);
	}

	int rid, Lid, Rid;
	while (1)
	{
		get_next_chunk_reads(data, Lid, Rid);
		if (Lid >= data->reads->num_reads) break;
		for (rid = Lid; rid < Rid; ++rid)
		{
			int rsize = data->reads->offset_list->offset_list[rid].size;
			extract_one_seq(data->reads, rid, read1);
			reverse_complement(read2, read1, rsize);
			int s;
			char chain;
			num_candidates = 0;
			for (s = 0; s < 2; ++s)
			{
				if (s%2) { chain = 'R'; read = read2; }
				else { chain = 'F'; read = read1; }
				int num_segs = seeding(read, rsize, data->ridx, sbk);
				num_candidates = get_candidates(data->reference,
												sbk,
												num_segs,
												rid + data->reads->start_read_id,
												rsize,
												chain,
												candidates,
												num_candidates);
			}

			for (s = 0; s < num_candidates; ++s)
			{
				if (candidates[s].chain == 'F') read = read1;
				else read = read2;
				extract_one_seq(data->reference, candidates[s].readno - data->reference->start_read_id, subject);
				int sstart = candidates[s].loc1;
				int qstart = candidates[s].loc2;
				if (qstart && sstart)
				{
					qstart += kmer_size / 2;
					sstart += kmer_size / 2;
				}
				int ssize = data->reference->offset_list->offset_list[candidates[s].readno - data->reference->start_read_id].size;

				int flag = aligner->go(read, qstart, rsize, subject, sstart, ssize, min_align_size);

				if (flag)
				{
					fill_m4record(aligner, rid + data->reads->start_read_id,
								  candidates[s].readno, candidates[s].chain,
								  rsize, ssize, qstart, sstart, candidates[s].score,
								  m4v + num_m4);
					++num_m4;
				}
			}

			append_m4v(m4_list, &m4_list_size, m4v, &num_m4, data->out, &data->result_write_lock);
		}
	}

		if (m4_list_size)
		{
			pthread_mutex_lock(&data->result_write_lock);
			print_m4record_list(data->out, m4_list, m4_list_size);
			m4_list_size = 0;
			pthread_mutex_unlock(&data->result_write_lock);
		}

		safe_free(read1);
		safe_free(read2);
		safe_free(subject);
		delete sbk;
		delete aligner;
		delete[] m4v;
}*/

void
candidate_detect(PWThreadData* data, int tid)
{
	char *read, *read1, *read2, *subject;
	safe_malloc(read1, char, MAX_SEQ_SIZE);
	safe_malloc(read2, char, MAX_SEQ_SIZE);
	safe_malloc(subject, char, MAX_SEQ_SIZE);
	SeedingBK* sbk = new SeedingBK(data->reference->curr);
	uint32_t read_index_count = 1 << (kmer_size * 2);
	read_index* read_index;//= create_read_index(read_index_count, READ_KMER_SIZE, POSITIONS_PER_KMER);
	Candidate candidates[MAXC];
	int num_candidates = 0;
	r_assert(data->ec_results);
	ExtensionCandidate* eclist = data->ec_results[tid];
	bool self_detection = data->self_detection;
	int nec = 0;
	ExtensionCandidate ec;

	int rid, Lid, Rid;
	while (1)
	{
		get_next_chunk_reads(data, Lid, Rid);
		if (Lid >= data->reads->num_reads) break;
		//if (Lid >=1000) break;
        for (rid = Lid; rid < Rid; ++rid)
        {
            int rsize = data->reads->offset_list->offset_list[rid].size;
            //int rstart= data->reads->offset_list->offset_list[rid].offset;
            if (rsize >= MAX_SEQ_SIZE) {
                cout << "rsize = " << rsize << "\t" << MAX_SEQ_SIZE << endl;
                abort();
            }
            extract_one_seq(data->reads, rid, read1);
            reverse_complement(read2, read1, rsize);
            int s;
            int chain;
            num_candidates = 0;
            if (rid==34)
                int y=0;
            for (s = 0; s < 2; ++s)
            {
                if (s%2) { chain = REV; read = read2; }
                else { chain = FWD; read = read1; }
                int num_segs;
                int rstart= data->reads->offset_list->offset_list[rid].offset;
                if (self_detection)
                    num_segs = self_seeding(read, rstart, rsize, data->reference->metaData, data->ridx, sbk);
                else
                    num_segs = seeding(read, rstart, rsize, data->reference->metaData, data->ridx, sbk);
                if (num_segs>1)
                    num_candidates = get_candidates(data->reference,
											sbk,
											read_index,
											read_index_count,
											num_segs,
											read,
											rid + data->reads->start_read_id,
											rsize,
											chain,
											candidates,
											num_candidates);
            }

            for (s = 0; s < num_candidates; ++s)
            {
                int qstart = candidates[s].loc2;
                int sstart = candidates[s].loc1;
                if (qstart && sstart)
                {
                    qstart += kmer_size / 2;
                    sstart += kmer_size / 2;
                }
                int qdir = candidates[s].chain;
                int sdir = FWD;
                int qid = rid + data->reads->start_read_id;
                int sid = candidates[s].readno;
                int score = candidates[s].score;

                ec.qid = qid;
                ec.qdir = qdir;
                ec.qext = qstart;
                ec.sid = sid;
                ec.sdir = sdir;
                ec.sext = sstart;
                ec.score = score;
                ec.qsize = rsize;
                ec.ssize = data->reference->offset_list->offset_list[sid - data->reference->start_read_id].size;
                if (ec.qdir == REV) ec.qext = ec.qsize - 1 - ec.qext;
                if (ec.sdir == REV) ec.sext = ec.ssize - 1 - ec.sext;
                eclist[nec] = ec;
                ++nec;
                if (nec == PWThreadData::kResultListSize)
                {
                    pthread_mutex_lock(&data->result_write_lock);
                    for (int i = 0; i < nec; ++i) (*data->out) << eclist[i];
                    nec = 0;
                    pthread_mutex_unlock(&data->result_write_lock);
                }
            }
        }
	}

	if (nec)
	{
		pthread_mutex_lock(&data->result_write_lock);
		for (int i = 0; i < nec; ++i) (*data->out) << eclist[i];
		nec = 0;
		pthread_mutex_unlock(&data->result_write_lock);
	}

	safe_free(read1);
	safe_free(read2);
	safe_free(subject);
	delete sbk;
	//destroy_read_index(read_index, read_index_count, POSITIONS_PER_KMER);
}

void*
multi_thread_func(void* p)
{
	PWThreadData* data = (PWThreadData*)p;
	int t = 0;
	pthread_mutex_lock(&data->id_lock);
	t = data->used_thread_id;
	++data->used_thread_id;
	pthread_mutex_unlock(&data->id_lock);
	r_assert(data->options->task == TASK_SEED || data->options->task == TASK_ALN);
	if (data->options->task == TASK_SEED) candidate_detect(data, t);
	else pairwise_mapping(data, t);
	return NULL;
}

void
process_one_volume(options_t* options, const int svid, const int evid, volume_names_t* vn, ostream* out)
{
	MAXC = options->num_candidates;
	output_gapped_start_point = options->output_gapped_start_point;
	min_align_size = options->min_align_size;
	min_kmer_match = options->min_kmer_match;

	if (options->tech == TECH_PACBIO) {
		ddfs_cutoff = ddfs_cutoff_pacbio;
		min_kmer_dist = 1800;
	} else if (options->tech == TECH_NANOPORE) {
		ddfs_cutoff = ddfs_cutoff_nanopore;
		min_kmer_dist = 400;
	} else {
		ERROR("TECH must be either %d or %d", TECH_PACBIO, TECH_NANOPORE);
	}

	const char* ref_name = get_vol_name(vn, svid);
	volume_t* ref = load_volume(ref_name);
	ref_index* ridx = create_ref_index(ref, kmer_size, options->num_threads);
	pthread_t tids[options->num_threads];
	char volume_process_info[1024];;
	int vid, tid;
	for(vid = svid; vid < evid; ++vid)
	{
		sprintf(volume_process_info, "process volume %d", vid);
		DynamicTimer dtimer(volume_process_info);
		const char* read_name = get_vol_name(vn, vid);
		LOG(stderr, "processing %s\n", read_name);
		volume_t* read = load_volume(read_name);
		PWThreadData* data = new PWThreadData(options, ref, read, ridx, out);
		if (vid == svid)
            data->self_detection=true;
		for (tid = 0; tid < options->num_threads; ++tid)
		{
			int err_code = pthread_create(tids + tid, NULL, multi_thread_func, (void*)data);
			if (err_code)
			{
				LOG(stderr, "Error: return code is %d\n", err_code);
				abort();
			}
		}
		for (tid = 0; tid < options->num_threads; ++tid) pthread_join(tids[tid], NULL);
		read = delete_volume_t(read);
		delete data;
	}
	ref = delete_volume_t(ref);
	ridx = destroy_ref_index(ridx);
}
