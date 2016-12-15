/*
 * Experiments with Wafe Function Collapse algorithms for procedural generation.
 *  - Allen Webster
 * Created 02.10.2016
 *
 * Thanks to "mxgmn" and the github page https://github.com/mxgmn/WaveFunctionCollapse
 *  for providing the initial concept.
 *
 */


#define DEBUGGING 5


#include <stdint.h>
#include <math.h>
#include <string.h>

#define ArrayCount(a) (sizeof(a)/sizeof(*(a)))

#if defined(DEBUGGING)
#include <stdio.h>
#include <stdlib.h>
#endif

typedef uint64_t TIME;
typedef struct{
    TIME start;
    TIME total;
    uint64_t start_count;
    uint64_t end_count;
} TIME_COUNTER;

static TIME_COUNTER null_time_counter = {0};
#define MAKE_TIME_COUNTER() (null_time_counter)

#if defined(DEBUGGING) && (DEBUGGING == 4 || DEBUGGING == 5)
# include <intrin.h>
# define BEGIN_TIME(v) TIME v = (__rdtsc())
# define END_TIME_MESSAGE(s,m) do { TIME t = __rdtsc() - s; printf("%20s: %12llu\n", m, t); } while(0)
# define END_TIME(s) END_TIME_MESSAGE(s,#s)

# if (DEBUGGING == 4)
#  define BEGIN_TIME_COUNTER(c) (c)->start_count++; (c)->start = __rdtsc()
#  define END_TIME_COUNTER(c) (c)->end_count++; (c)->total += (__rdtsc() - (c)->start)
#  define DISPLAY_TIME_COUNTER_MESSAGE(c,m) do {                  \
    if ((c)->start_count == (c)->end_count)                      \
    printf("%20s: %12llu\n%44s: %10llu\n%44s: %10llu\n",         \
    m, (c)->total,                                        \
    "average", (c)->total / (c)->start_count,             \
     "count", (c)->start_count);                           \
    else printf("%30s: COUNT MISMATCH ERROR %llu vs  %llu\n",    \
    m, (c)->start_count, (c)->end_count);            \
} while(0)
#  define DISPLAY_TIME_COUNTER(c) DISPLAY_TIME_COUNTER_MESSAGE(c,#c)
# endif

#endif

#if !defined(BEGIN_TIME)
# define BEGIN_TIME(v)
# define END_TIME_MESSAGE(s,m)
# define END_TIME(s)
#endif

#if !defined(BEGIN_TIME_COUNTER)
# define BEGIN_TIME_COUNTER(c)
# define END_TIME_COUNTER(c)
# define DISPLAY_TIME_COUNTER_MESSAGE(c,m)
# define DISPLAY_TIME_COUNTER(c)
#endif

// RANDOM
typedef struct Random{
    uint64_t state;
    uint64_t inc;
} Random;

// STUDY(allen): LEARN THIS
static uint32_t
pcg32_random(Random* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((~rot + 1) & 31));
}

static float
pcg32_random_float_01(Random *rng){
    float result = (pcg32_random(rng) % 1024)/1024.f;
    return(result);
}

// HARD CODED TEST CASE

#define SRC_IMG_W 11
#define SRC_IMG_H 11

uint32_t test_base_image[SRC_IMG_W][SRC_IMG_H] = {
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
};

// MEMORY
typedef struct Partition{
    uint8_t *data;
    int32_t pos;
    int32_t max;
} Partition;

static Partition
part_make(void *data, int32_t max){
    Partition part = {0};
    part.data = (uint8_t*)data;
    part.max = max;
    return(part);
}

static void*
part_push(Partition *part, int32_t size){
    void *result = 0;
    if (part->pos + size <= part->max){
        result = part->data + part->pos;
        part->pos += size;
    }
    return(result);
}

static void*
part_current(Partition *part){
    void *result = 0;
    if (part->pos < part->max){
        result = part->data + part->pos;
    }
    return(result);
}

static void
part_clear(Partition *part){
    part->pos = 0;
}

#define part_begin_array(p,T) (T*)(part_current(p))
#define push_struct(p,T) (T*)(part_push((p), sizeof(T)))
#define push_array(p,T,c) (T*)(part_push((p), sizeof(T)*(c)))

// 2D WAVE FUNCTION CODE
typedef struct Wave2D_Sample{
    uint32_t *data;
    float weight;
} Wave2D_Sample;

typedef struct Wave2D_Samples{
    Partition part;
    
    uint32_t w;
    uint32_t h;
    
    Wave2D_Sample *samples;
    int32_t sample_count;
} Wave2D_Samples;

static void
wave2d_samples_memory(Wave2D_Samples *samples, void *memory, int32_t size){
    samples->part = part_make(memory, size);
}

static void
wave2d_begin_samples(Wave2D_Samples *samples, uint32_t w, uint32_t h){
    samples->w = w;
    samples->h = h;
    
    samples->samples = part_begin_array(&samples->part, Wave2D_Sample);
}

static void
wave2d_add_sample(Wave2D_Samples *samples, uint32_t *data, float weight){
    Wave2D_Sample *sample = push_struct(&samples->part, Wave2D_Sample);
    sample->data = data;
    sample->weight = weight;
    ++samples->sample_count;
}

static void
wave2d_end_samples(Wave2D_Samples *samples){}

typedef struct Wave2D_State{
    Partition part;
    
    uint32_t output_w, output_h;
    uint32_t sample_count;
    uint32_t *feasible_states;
    
    float log_T;
    float total_weight;
    float *weight_by_log_weight;
    
    Wave2D_Samples *samples;
} Wave2D_State;

static void
wave2d_state_memory(Wave2D_State *state, void *memory, int32_t size){
    state->part = part_make(memory, size);
}

static void
wave2d_initialize_state(Wave2D_State *state, Wave2D_Samples *samples, uint32_t output_w, uint32_t output_h){
    uint32_t grid_cell_count = output_w*output_h;
    uint32_t sample_count = samples->sample_count;
    uint32_t feasible_state_count = grid_cell_count * (sample_count + 1);
    
    state->output_w = output_w;
    state->output_h = output_h;
    state->sample_count = sample_count;
    state->samples = samples;
    
    part_clear(&state->part);
    state->feasible_states = push_array(&state->part, uint32_t, feasible_state_count);
    
    state->log_T = logf((float)sample_count);
    state->weight_by_log_weight = push_array(&state->part, float, sample_count);
    
    float *weight_by_log_weight = state->weight_by_log_weight;
    float total_weight = 0.f;
        float *weight_ptr = &samples->samples[0].weight;
    int32_t stride = sizeof(Wave2D_Sample);
        for (uint32_t j = 0; j < sample_count; ++j){
            *(weight_by_log_weight++) = (*weight_ptr) * logf(*weight_ptr);
            total_weight += *weight_ptr;
            weight_ptr = (float*)((uint8_t*)weight_ptr + stride);
        }
        
        state->total_weight = total_weight;
    }

typedef struct Wave2D_Change{
    uint32_t x, y;
} Wave2D_Change;

enum Wave2D_Add_Change_Result{
    AddChange_Failed,
    AddChange_Success,
    AddChange_Success_Already_On_Queue,
};

static int32_t
wave2d_add_change(Wave2D_Change change, Wave2D_Change *changes, int32_t *change_count, int32_t change_max, uint8_t *on_change_queue, uint32_t w){
    int32_t result = AddChange_Failed;
    uint32_t i = change.x + change.y*w;
    if (on_change_queue[i] == 0){
        if (*change_count < change_max){
            result = AddChange_Success;
            changes[(*change_count)++] = change;
        }
    }
    else{
        result = AddChange_Success_Already_On_Queue;
    }
    return(result);
}

enum{
    Fail_ValidSampleCountZero,
    Fail_DieRollBug,
    Fail_InitialAddChangeBug,
    Fail_ImpossibleToResolve,
    Fail_ChangeQueueOutOfMemory,
};

typedef struct Wave2D_Failure_Info{
    uint32_t type;
    
    union{
        struct{
            uint32_t x, y;
        } contradiction;
    };
} Wave2D_Failure_Info;

#if defined(DEBUGGING) && (DEBUGGING != 0)
static void
wave2d_get_debug_output(uint32_t *map, uint32_t w, uint32_t h, uint8_t *cell_finished, uint32_t *feasible_states, uint32_t sample_count, Wave2D_Samples *samples){
    memset(map, 0xFF, 4*w*h);
    for (uint32_t y = 0; y < h; ++y){
        for (uint32_t x = 0; x < w; ++x){
            uint32_t *cell_base = &feasible_states[(x + y*w)*(sample_count+1)];
            if (*cell_base == 1){
                uint32_t *cell = cell_base + 1;
                uint32_t t = cell[0];
                
                uint32_t *sample = samples->samples[t].data;
                for (uint32_t yy = 0; yy < samples->h; ++yy){
                    for (uint32_t xx = 0; xx < samples->w; ++xx){
                        uint32_t xx_x = (xx + x) % w;
                        uint32_t yy_y = (yy + y) % h;
                        map[xx_x + yy_y*w] = sample[xx + yy*samples->w];
                    }
                }
            }
        }
    }
}
#endif

static void*
round_ptr(void *ptr, uint32_t boundary){
    uint64_t pos = (uint64_t)ptr;
    --boundary;
    pos = (pos+boundary)&(~boundary);
    ptr = (void*)pos;
    return(ptr);
}

static int32_t
wave2d_generate_output(Wave2D_State *state, Random *rng, uint32_t *out, void *scratch, int32_t scratch_size){
#if defined(DEBUGGING) && (DEBUGGING == 4)
    TIME_COUNTER observation = MAKE_TIME_COUNTER();
    TIME_COUNTER propogation = MAKE_TIME_COUNTER();
    TIME_COUNTER single_change = MAKE_TIME_COUNTER();
    TIME_COUNTER dxy_preproc = MAKE_TIME_COUNTER();
    TIME_COUNTER full_check = MAKE_TIME_COUNTER();
    TIME_COUNTER coexist_check = MAKE_TIME_COUNTER();
    TIME_COUNTER add_change = MAKE_TIME_COUNTER();
    #endif
    
    int32_t result = 0;
    
    uint32_t output_w = state->output_w;
    uint32_t output_h = state->output_h;
    uint32_t cell_count = output_w*output_h;
    uint32_t sample_count = state->sample_count;
    
    Wave2D_Samples *samples = state->samples;
    
    uint32_t sample_w = samples->w;
    uint32_t sample_h = samples->h;
    
    uint32_t *feasible_states = state->feasible_states;
    float *weight_by_log_weight = state->weight_by_log_weight;
    
    uint8_t *scratch_start = scratch;
    uint8_t *scratch_ptr = scratch_start;
    uint8_t *scratch_end = scratch_ptr + scratch_size;
    (void)scratch_end;
    
    uint8_t *on_change_queue = (uint8_t*)scratch_ptr;
    memset(on_change_queue, 0, cell_count);
    scratch_ptr += cell_count;
    
    uint8_t *cell_finished = (uint8_t*)scratch_ptr;
    memset(cell_finished, 0, cell_count);
    scratch_ptr += cell_count;
    
    scratch_ptr = round_ptr(scratch_ptr, 4);
    
    Wave2D_Change *changes = (Wave2D_Change*)scratch_ptr;
    memset(changes, 0, sizeof(Wave2D_Change)*cell_count);
    int32_t change_max = cell_count;
    int32_t change_count = 0;
    scratch_ptr += sizeof(Wave2D_Change)*cell_count;
    
#if defined(DEBUGGING) && (DEBUGGING != 0)
    if (scratch_ptr > scratch_end){
        *(int*)0 = 0xA11E;
    }
    #endif
    
    // CLEAR
    uint32_t *feasible_state_ptr = feasible_states;
    uint32_t cell_stride = sample_count + 1;
    for (uint32_t i = 0; i < cell_count; ++i){
        *(feasible_state_ptr++) = sample_count;
        for (uint32_t j = 0; j < sample_count; ++j){
            *(feasible_state_ptr++) = j;
        }
    }
    
    // COLLAPSE
    uint32_t failed = 0;
    Wave2D_Failure_Info fail_info = {0};
    for (int32_t step = 0; ; ++step){
        // OBSERVE
        //  find smallest entropy
        
        BEGIN_TIME_COUNTER(&observation);
        
        uint32_t best_x = output_w;
        uint32_t best_y = output_h;
        float best_total_weight = 0;
        float min_entropy = 100000.f;
        
        for (uint32_t x = 0; x < output_w; ++x){
            for (uint32_t y = 0; y < output_h; ++y){
                
                if (!cell_finished[(x + y*output_w)]){
                    uint32_t *cell_feasible_states_base = &feasible_states[(x + y*output_w) * cell_stride];
                    uint32_t *cell_feasible_states = cell_feasible_states_base+1;
                    uint32_t *cell_feasible_count = cell_feasible_states_base;
                    
                    uint32_t valid_sample_count = *cell_feasible_count;
                    float total_weight = 0.f;
                    float entropy = 0.f;
                    
                    // TODO(allen): Is it really faster this way?
                    // I guess the point is to avoid too many array lookups
                    // in common short-circuitable cases?  Is this actually
                    // speeding it up for us?
                    if (valid_sample_count == 0){
                        // TODO(allen): This can be an assert once we start using that.
                        failed = 1;
                        fail_info.type = Fail_ValidSampleCountZero;
                        goto finished;
                    }
                    else if (valid_sample_count == sample_count){
                        entropy = state->log_T;
                        total_weight = state->total_weight;
                    }
                    else{
                        for (uint32_t t_i = 0; t_i < valid_sample_count; ++t_i){
                            uint32_t t = cell_feasible_states[t_i];
                            total_weight += samples->samples[t].weight;
                        }
                        
                        float entropy_sum = 0;
                        for (uint32_t t_i = 0; t_i < valid_sample_count; ++t_i){
                            uint32_t t = cell_feasible_states[t_i];
                                entropy_sum += weight_by_log_weight[t];
                        }
                        
                        float log_total_weight = logf(total_weight);
                        entropy = log_total_weight - entropy_sum/total_weight;
                    }
                    
                    if (entropy > 0){
                        float noise = 0.0000001f*pcg32_random_float_01(rng);
                        entropy += noise;
                        if (entropy < min_entropy){
                            min_entropy = entropy;
                            best_x = x;
                            best_y = y;
                            best_total_weight = total_weight;
                        }
                    }
                }
            }
        }
        
        if (best_x == output_w && best_y == output_h){
            END_TIME_COUNTER(&observation);
            goto finished;
        }
        
        // randomly choose a value for the cell with lowest entropy
        cell_finished[(best_x + best_y*output_w)] = 1;
        
        uint32_t *cell_feasible_states_base = &feasible_states[(best_x + best_y*output_w) * cell_stride];
        uint32_t *cell_feasible_states = cell_feasible_states_base+1;
        uint32_t valid_sample_count = *cell_feasible_states_base;
        
        float die_roll = best_total_weight * pcg32_random_float_01(rng);
        
        uint32_t hit_element = sample_count;
        for (uint32_t t_i = 0; t_i < valid_sample_count; ++t_i){
            uint32_t t = cell_feasible_states[t_i];
            die_roll -= samples->samples[t].weight;
            if (die_roll < 0){
                hit_element = t;
                break;
            }
        }
        
        if (hit_element == sample_count){
            failed = 1;
            fail_info.type = Fail_DieRollBug;
            goto finished;
        }
        
        *cell_feasible_states_base = 1;
        cell_feasible_states[0] = hit_element;
        
        {
            Wave2D_Change new_change;
            new_change.x = best_x;
            new_change.y = best_y;
            int32_t added_change = wave2d_add_change(new_change, changes, &change_count, change_max, on_change_queue, output_w);
            
            if (added_change != AddChange_Success){
                failed = 1;
                fail_info.type = Fail_InitialAddChangeBug;
                goto finished;
            }
        }
        
        END_TIME_COUNTER(&observation);
        BEGIN_TIME_COUNTER(&propogation);
        
        // PROPOGATE
        for (int32_t i = 0; i < change_count; ++i){
            BEGIN_TIME_COUNTER(&single_change);
            Wave2D_Change change = changes[i];
            
            uint32_t *cell1_base = &feasible_states[(change.x + change.y*output_w)*cell_stride];
            uint32_t *cell1 = cell1_base + 1;
            
            for (int32_t delta_x = -(int32_t)(sample_w) + 1; delta_x < (int32_t)(sample_w); ++delta_x){
                for (int32_t delta_y = -(int32_t)(sample_h) + 1; delta_y < (int32_t)(sample_h); ++delta_y){
                    BEGIN_TIME_COUNTER(&dxy_preproc);
                    if (delta_x == 0 && delta_y == 0){
                        END_TIME_COUNTER(&dxy_preproc);
                        continue;
                    }
                    
                    uint32_t x2 = (uint32_t)((change.x + delta_x + output_w) % output_w);
                    uint32_t y2 = (uint32_t)((change.y + delta_y + output_h) % output_h);
                    
                    // cell1 changed, now propogate to cell2
                    int32_t minx = 0, maxx = (int32_t)sample_w;
                    int32_t miny = 0, maxy = (int32_t)sample_h;
                    
                    if (minx < delta_x){
                        minx = delta_x;
                    }
                    
                    if (maxx > delta_x + (int32_t)sample_w){
                        maxx = delta_x + (int32_t)sample_w;
                    }
                    
                    if (miny < delta_y){
                        miny = delta_y;
                    }
                    
                    if (maxy > delta_y + (int32_t)sample_h){
                        maxy = delta_y + (int32_t)sample_h;
                    }
                    
                    uint32_t cell_2_changed = 0;
                    uint8_t cell2_is_finished = cell_finished[x2 + y2*output_w];
                    uint32_t *cell2_base = &feasible_states[(x2 + y2*output_w)*cell_stride];
                    uint32_t *cell2 = cell2_base + 1;
                    
                    END_TIME_COUNTER(&dxy_preproc);
                    
                    BEGIN_TIME_COUNTER(&full_check);
                    END_TIME_COUNTER(&full_check);
                    
                    for (uint32_t t2_i = 0; t2_i < *cell2_base; ++t2_i){
                        uint32_t t2 = cell2[t2_i];
                        
                            uint32_t can_coexist = 0;
                        for (uint32_t t1_i = 0; t1_i < *cell1_base; ++t1_i){
                            BEGIN_TIME_COUNTER(&coexist_check);
                                uint32_t t1 = cell1[t1_i];
                                    uint32_t coexist_check_result = 1;
                                    
                                    // COEXIST CHECK
                                    uint32_t *sample_1 = samples->samples[t1].data;
                                    uint32_t *sample_2 = samples->samples[t2].data;
                                    
                                    for (int32_t sx = minx; sx < maxx; ++sx){
                                        for (int32_t sy = miny; sy < maxy; ++sy){
                                            if (sample_1[sx + sy*sample_w] != sample_2[(sx-delta_x) + (sy-delta_y)*sample_w]){
                                                coexist_check_result = 0;
                                                goto finish_coexist_check;
                                            }
                                        }
                                    }
                                    finish_coexist_check:;
                                    
                                    if (coexist_check_result){
                                        can_coexist = 1;
                                        END_TIME_COUNTER(&coexist_check);
                                        break;
                                    }
                                    END_TIME_COUNTER(&coexist_check);
                            }
                            
                            if (!can_coexist){
                                if (cell2_is_finished){
                                    failed = 1;
                                    fail_info.type = Fail_ImpossibleToResolve;
                                    fail_info.contradiction.x = change.x;
                                    fail_info.contradiction.y = change.y;
                                    goto finished;
                                }
                                
                                --(*cell2_base);
                                cell2[t2_i] = cell2[(*cell2_base)];
                                --t2_i;
                                
                                cell_2_changed = 1;
                            }
                    }
                    
                    if (!cell2_is_finished && *cell2_base == 1){
                        cell_finished[x2 + y2*output_w] = 1;
                    }
                    
                    if (cell_2_changed){
                        BEGIN_TIME_COUNTER(&add_change);
                        if (change_count == change_max){
                            int32_t n = i;
                            memmove(changes, changes + n, change_count - n);
                            change_count -= n;
                            i = 0;
                        }
                        
                        Wave2D_Change new_change;
                        new_change.x = x2;
                        new_change.y = y2;
                        
                        int32_t added_change = wave2d_add_change(new_change, changes, &change_count, change_max, on_change_queue, output_w);
                        if (added_change == AddChange_Failed){
                            failed = 1;
                            fail_info.type = Fail_ChangeQueueOutOfMemory;
                            END_TIME_COUNTER(&add_change);
                            goto finished;
                        }
                        END_TIME_COUNTER(&add_change);
                    }
                }
            }
            
            on_change_queue[change.x + change.y*output_w] = 0;
            END_TIME_COUNTER(&single_change);
        }
        
        change_count = 0;
        
        END_TIME_COUNTER(&propogation);
        
#if defined(DEBUGGING) && (DEBUGGING == 1)
        for (uint32_t y = 0; y < output_h; ++y){
            for (uint32_t x = 0; x < output_w; ++x){
                if (cell_finished[x + y*output_w]){
                    uint8_t *cell = &coefficients[(x + y*output_w)*sample_count];
                    uint32_t t = 0;
                    for (; t < sample_count; ++t){
                        if (cell[t]) break;
                    }
                    printf("%u", t);
                }
                else{
                    printf("?");
                }
            }
            printf("\n");
        }
        printf("\n");
#endif
        
#if defined(DEBUGGING) && (DEBUGGING == 3)
        {
            uint32_t *map = (uint32_t*)malloc(4*output_w*output_h);
            wave2d_get_debug_output(map, output_w, output_h, cell_finished, feasible_states, sample_count, samples);
            
        for (uint32_t y = 0; y < output_h; ++y){
            for (uint32_t x = 0; x < output_w; ++x){
                uint32_t v = map[x + y*output_w];
                if (v == 0xFFFFFFFF){
                    printf("?");
                }
                else{
                    printf("%u", v);
                }
            }
            printf("\n");
        }
        printf("\n");
        free(map);
    }
#endif
    }
    
    finished:;
    DISPLAY_TIME_COUNTER(&observation);
    DISPLAY_TIME_COUNTER(&propogation);
    DISPLAY_TIME_COUNTER(&single_change);
    DISPLAY_TIME_COUNTER(&dxy_preproc);
    DISPLAY_TIME_COUNTER(&full_check);
    DISPLAY_TIME_COUNTER(&coexist_check);
    DISPLAY_TIME_COUNTER(&add_change);
    DISPLAY_TIME_COUNTER(&update_out_grid);
    
    if (!failed){
        for (uint32_t x = 0; x < output_w; ++x){
            for (uint32_t y = 0; y < output_h; ++y){
                uint32_t *cell_base = &feasible_states[(x + y*output_w)*cell_stride];
                uint32_t *cell = cell_base + 1;
                uint32_t t = cell[0];
                uint32_t *sample = samples->samples[t].data;
                out[x + y*output_w] = sample[0];
            }
        }
        result = 1;
    }
    else{
#if defined(DEBUGGING) && (DEBUGGING != 0)
        switch (fail_info.type){
            case Fail_ValidSampleCountZero:
            {
                printf("Found a cell with a sample count of zero.\n");
            }break;
            
            case Fail_DieRollBug:
            {
                printf("Had an invalid die roll durring sample selection.\n");
            }break;
            
            case Fail_InitialAddChangeBug:
            {
                printf("Had a bug trying to seed the change queue after random resolution stage.\n");
            }break;
            
            case Fail_ImpossibleToResolve:
            {
                printf("Reached a contradictive state at some tile on the board.\n");
                
                uint32_t *map = (uint32_t*)malloc(4*output_w*output_h);
                wave2d_get_debug_output(map, output_w, output_h, cell_finished, feasible_states, sample_count, samples);
                
                for (uint32_t y = 0; y < sample_h; ++y){
                    for (uint32_t x = 0; x < sample_w; ++x){
                        uint32_t yy = fail_info.contradiction.y + y;
                        uint32_t xx = fail_info.contradiction.x + x;
                        
                        yy %= output_h;
                        xx %= output_w;
                        
                        map[xx + yy*output_w] = 0xEEEEEEEE;
                    }
                }
                
                for (uint32_t y = 0; y < output_h; ++y){
                    for (uint32_t x = 0; x < output_w; ++x){
                        uint32_t v = map[x + y*output_w];
                        if (v == 0xFFFFFFFF){
                            printf("?");
                        }
                        else if (v == 0xEEEEEEEE){
                            printf("X");
                        }
                        else{
                            printf("%u", v);
                        }
                    }
                    printf("\n");
                }
                printf("\n");
                free(map);
            }break;
            
            case Fail_ChangeQueueOutOfMemory:
            {
                printf("Could not push a new change onto the change queue because it was out of memory.\n");
            }break;
        }
        #endif
    }
    
    return(result);
}

// SET SAMPLE FROM INPUT IMAGE
typedef struct Wave2D_Image_Processing_Params{
    uint32_t w;
    uint32_t h;
    uint32_t window_w;
    uint32_t window_h;
    
    uint32_t sample_x_positions;
    uint32_t sample_y_positions;
    uint32_t max_sample_count;
    uint32_t scratch_size;
    
    uint32_t next_sample_stride;
    
    uint32_t wrapped_image_w;
    uint32_t wrapped_image_h;
    uint32_t wrapped_image_size;
    
    uint32_t apply_rotation;
    uint32_t apply_mirroring;
} Wave2D_Image_Processing_Params;

static Wave2D_Image_Processing_Params
wave2d_image_processing_params(uint32_t w, uint32_t h, uint32_t window_w, uint32_t window_h, int32_t wrapped, uint32_t rotation_level, uint32_t mirror_type){
    Wave2D_Image_Processing_Params params;
    params.w = w;
    params.h = h;
    params.window_w = window_w;
    params.window_h = window_h;
    
    if (wrapped){
        params.sample_x_positions = params.w;
        params.sample_y_positions = params.h;
        params.max_sample_count = params.sample_x_positions*params.sample_y_positions;
        
        params.wrapped_image_w = w + window_w - 1;
        params.wrapped_image_h = h + window_h - 1;
        params.wrapped_image_size = params.wrapped_image_w*params.wrapped_image_h;
        
        params.next_sample_stride = params.window_w*params.window_h;
    }
    else{
        params.sample_x_positions = params.w - params.window_w + 1;
        params.sample_y_positions = params.h - params.window_h + 1;
        params.max_sample_count = params.sample_x_positions*params.sample_y_positions;
        
        params.next_sample_stride = params.window_w*params.window_h;
        params.wrapped_image_w = 0;
        params.wrapped_image_h = 0;
        params.wrapped_image_size = 0;
    }
    
    if (rotation_level == 1){
        params.max_sample_count *= 2;
        params.apply_rotation = 1;
    }
    else if (rotation_level == 2){
        params.max_sample_count *= 4;
        params.apply_rotation = 2;
    }
    else{
        params.apply_rotation = 0;
    }
    
    if (mirror_type == 1){
        params.max_sample_count *= 2;
        params.apply_mirroring = 1;
    }
    else if (mirror_type == 2){
        params.max_sample_count *= 2;
        params.apply_mirroring = 2;
    }
    else{
        params.apply_mirroring = 0;
    }
    
    params.scratch_size = params.max_sample_count*(params.w*params.h*4 + 1) + 4*(params.wrapped_image_size);
    
    return(params);
}

static uint32_t*
wave2d_get_processed_samples_from_scratch(Wave2D_Image_Processing_Params params, void *scratch_memory){
    uint32_t *result = (uint32_t*)(scratch_memory);
    return(result);
}

static uint8_t*
wave2d_processed_counts_from_scratch(Wave2D_Image_Processing_Params params, void *scratch_memory){
    uint32_t *sample_memory = (uint32_t*)scratch_memory;
    uint32_t *wrapped_image = (uint32_t*)(sample_memory + params.max_sample_count*params.next_sample_stride);
    uint8_t *result = (uint8_t*)(wrapped_image + params.wrapped_image_size);
    return(result);
}

static void
wave2d_image_to_sample_blit(uint32_t *image, uint32_t image_stride, uint32_t x, uint32_t y, uint32_t *output, uint32_t output_w, uint32_t output_h){
    uint32_t *src_line = image + (x + y*image_stride);
    uint32_t *dst_line = output;
    for (uint32_t y = 0; y < output_h; ++y){
        uint32_t *src = src_line;
        uint32_t *dst = dst_line;
        for (uint32_t x = 0; x < output_w; ++x){
            *dst = *src;
            ++dst;
            ++src;
        }
        src_line += image_stride;
        dst_line += output_w;
    }
}

static void
wave2d_rotated_sample_to_sample_blit(uint32_t *src_sample, uint32_t rotation, uint32_t *output, uint32_t output_w, uint32_t output_h){
    switch (rotation % 4){
        case 0:
        {
            // TODO(allen): this could be memcpy
            uint32_t *src_line = src_sample;
            uint32_t *dst_line = output;
            for (uint32_t y = 0; y < output_h; ++y){
                uint32_t *src = src_line;
                uint32_t *dst = dst_line;
                for (uint32_t x = 0; x < output_w; ++x){
                    *dst = *src;
                    ++dst;
                    ++src;
                }
                src_line += output_w;
                dst_line += output_w;
            }
        }break;
        
        case 1:
        {
            uint32_t *src_line = src_sample;
            uint32_t *dst_col = output + output_w - 1;
            for (uint32_t y = 0; y < output_h; ++y){
                uint32_t *src = src_line;
                uint32_t *dst = dst_col;
                for (uint32_t x = 0; x < output_w; ++x){
                    *dst = *src;
                    dst += output_w;
                    ++src;
                }
                src_line += output_w;
                --dst_col;
            }
        }break;
        
        case 2:
        {
            // TODO(allen): this could be backcpy
            uint32_t *src_line = src_sample;
            uint32_t *dst_line = output + output_w*output_h - 1;
            for (uint32_t y = 0; y < output_h; ++y){
                uint32_t *src = src_line;
                uint32_t *dst = dst_line;
                for (uint32_t x = 0; x < output_w; ++x){
                    *dst = *src;
                    --dst;
                    ++src;
                }
                src_line += output_w;
                dst_line -= output_w;
            }
        }break;
        
        case 3:
        {
            uint32_t *src_line = src_sample;
            uint32_t *dst_col = output + (output_h-1)*output_w;
            for (uint32_t y = 0; y < output_h; ++y){
                uint32_t *src = src_line;
                uint32_t *dst = dst_col;
                for (uint32_t x = 0; x < output_w; ++x){
                    *dst = *src;
                    dst -= output_w;
                    ++src;
                }
                src_line += output_w;
                ++dst_col;
            }
        }break;
    }
}

static void
wave2d_mirrored_sample_to_sample_blit(uint32_t *src_sample, uint32_t vertical_mirror, uint32_t *output, uint32_t output_w, uint32_t output_h){
    switch (vertical_mirror % 2){
        case 0:
        {
            uint32_t *src_line = src_sample;
            uint32_t *dst_line = output + output_w - 1;
            for (uint32_t y = 0; y < output_h; ++y){
                uint32_t *src = src_line;
                uint32_t *dst = dst_line;
                for (uint32_t x = 0; x < output_w; ++x){
                    *dst = *src;
                    --dst;
                    ++src;
                }
                src_line += output_w;
                dst_line += output_w;
            }
        }break;
        
        case 1:
        {
            uint32_t *src_line = src_sample;
            uint32_t *dst_line = output + (output_h-1)*output_w;
            for (uint32_t y = 0; y < output_h; ++y){
                uint32_t *src = src_line;
                uint32_t *dst = dst_line;
                for (uint32_t x = 0; x < output_w; ++x){
                    *dst = *src;
                    ++dst;
                    ++src;
                }
                src_line += output_w;
                dst_line -= output_w;
            }
        }break;
    }
}

typedef struct {
    int32_t is_repeat_pattern;
    uint32_t repeat_index;
} Wave2D_Sample_Check_Result;

static Wave2D_Sample_Check_Result
wave2d_check_new_sample(Wave2D_Image_Processing_Params params, uint32_t *existing_samples, uint32_t sample_count, uint32_t *new_sample){
    Wave2D_Sample_Check_Result check_result = {0};
    
    uint32_t *end_sample = existing_samples + sample_count*params.next_sample_stride;
    
    for (uint32_t *prev_sample = existing_samples;
         prev_sample < end_sample;
         prev_sample += params.next_sample_stride, ++check_result.repeat_index){
        
        check_result.is_repeat_pattern = 1;
        uint32_t *a_ptr = prev_sample;
        uint32_t *b_ptr = new_sample;
        for (uint32_t i = 0; i < params.next_sample_stride; ++i){
            if ((*a_ptr) != (*b_ptr)){
                check_result.is_repeat_pattern = 0;
                break;
            }
            ++a_ptr;
            ++b_ptr;
        }
        
        if (check_result.is_repeat_pattern){
            break;
        }
    }
    
    return(check_result);
}

static uint32_t
wave2d_extract_samples_from_image(Wave2D_Image_Processing_Params params, uint32_t *image, void *scratch_memory){
    uint32_t *sample_memory = (uint32_t*)scratch_memory;
    uint32_t *wrapped_image = (uint32_t*)(sample_memory + params.max_sample_count*params.next_sample_stride);
    uint8_t *use_count = (uint8_t*)(wrapped_image + params.wrapped_image_size);
    
    // Make wrappable image
    uint32_t image_stride = params.w;
    if (params.wrapped_image_size > 0){
        uint32_t i = 0;
        for (uint32_t y = 0; y < params.wrapped_image_h; ++y){
            for (uint32_t x = 0; x < params.wrapped_image_w; ++x){
                uint32_t yy = y % params.h;
                uint32_t xx = x % params.w;
                wrapped_image[i] = image[xx + yy*params.w];
                ++i;
            }
        }
        image = wrapped_image;
        image_stride = params.wrapped_image_w;
    }
    
    // Extract samples from image and add to sample memory
    uint32_t *current_sample = sample_memory;
    uint8_t *current_count = use_count;
    for (uint32_t y = 0; y < params.sample_y_positions; ++y){
        for (uint32_t x = 0; x < params.sample_x_positions; ++x){
            wave2d_image_to_sample_blit(image, image_stride, x, y, current_sample, params.window_w, params.window_h);
            
            Wave2D_Sample_Check_Result check_result = 
                wave2d_check_new_sample(params, sample_memory, (uint32_t)(current_count - use_count), current_sample);
            
            if (check_result.is_repeat_pattern){
                ++use_count[check_result.repeat_index];
            }
            else{
                *current_count = 1;
                current_sample += params.next_sample_stride;
                ++current_count;
            }
        }
    }
    
    // Perform rotations
    if (params.apply_rotation){
    uint32_t total_count = (uint32_t)(current_count - use_count);
        
        uint32_t rot_level_1[1] = {2};
        uint32_t rot_level_2[3] = {1,2,3};
        
        uint32_t *rotation_array = 0;
        uint32_t rotation_count = 0;
        
        if (params.apply_rotation == 1){
            rotation_array = rot_level_1;
            rotation_count = ArrayCount(rot_level_1);
        }
        else if (params.apply_rotation == 2){
            rotation_array = rot_level_2;
            rotation_count = ArrayCount(rot_level_2);
        }
        
        if (rotation_count != 0){
    for (uint32_t rotation_index = 0; rotation_index < rotation_count; ++rotation_index){
        uint32_t rotation = rotation_array[rotation_index];
        uint32_t *end_sample = sample_memory + total_count*params.next_sample_stride;
        for (uint32_t *prev_sample = sample_memory;
             prev_sample < end_sample;
             prev_sample += params.next_sample_stride){
            wave2d_rotated_sample_to_sample_blit(prev_sample, rotation, current_sample, params.window_w, params.window_h);
            
            Wave2D_Sample_Check_Result check_result = 
                wave2d_check_new_sample(params, sample_memory, (uint32_t)(current_count - use_count), current_sample);
            
            if (check_result.is_repeat_pattern){
                ++use_count[check_result.repeat_index];
            }
            else{
                *current_count = 1;
                current_sample += params.next_sample_stride;
                ++current_count;
            }
    }
}
}
}
    
// Perform mirroring
if (params.apply_mirroring){
    uint32_t total_count = (uint32_t)(current_count - use_count);
    
        uint32_t *end_sample = sample_memory + total_count*params.next_sample_stride;
        for (uint32_t *prev_sample = sample_memory;
             prev_sample < end_sample;
             prev_sample += params.next_sample_stride){
            wave2d_mirrored_sample_to_sample_blit(prev_sample, params.apply_mirroring, current_sample, params.window_w, params.window_h);
            
            Wave2D_Sample_Check_Result check_result = 
                wave2d_check_new_sample(params, sample_memory, (uint32_t)(current_count - use_count), current_sample);
            
            if (check_result.is_repeat_pattern){
                ++use_count[check_result.repeat_index];
            }
            else{
                *current_count = 1;
                current_sample += params.next_sample_stride;
                ++current_count;
            }
        }
    }

    uint32_t total_count = (uint32_t)(current_count - use_count);
    return(total_count);
}

#if defined(DEBUGGING)
static void
wave2d_print_samples(Wave2D_Samples samples){
    Wave2D_Sample *sample_ptr = samples.samples;
    Wave2D_Sample *end_ptr = samples.samples + samples.sample_count;
    
    for (; sample_ptr < end_ptr; ++sample_ptr){
        printf("weight %f:\n", sample_ptr->weight);
        uint32_t i = 0;
        for (uint32_t y = 0; y < samples.h; ++y){
            for (uint32_t x = 0; x < samples.w; ++x){
                printf("%u", sample_ptr->data[i]);
                ++i;
            }
            printf("\n");
        }
        printf("\n");
    }
}
#endif

static void
print_dst(uint32_t *dst){
    uint32_t i = 0;
    for (uint32_t y = 0; y < 3; ++y){
        for (uint32_t x = 0; x < 3; ++x){
            printf("%u", dst[i]);
            ++i;
        }
        printf("\n");
    }
    printf("\n");
}

static void
rotation_test(){
    uint32_t src_[3][3] = {
        {1, 2, 3}, 
        {4, 5, 6}, 
        {7, 8, 9}, 
    };
    uint32_t *src = &src_[0][0];
    
    uint32_t dst_[3][3] = {
        {0, 0, 0}, 
        {0, 0, 0}, 
        {0, 0, 0}, 
    };
    uint32_t *dst = &dst_[0][0];
    
    print_dst(dst);
    wave2d_rotated_sample_to_sample_blit(src, 0, dst, 3, 3);
    print_dst(dst);
    wave2d_rotated_sample_to_sample_blit(src, 1, dst, 3, 3);
    print_dst(dst);
    wave2d_rotated_sample_to_sample_blit(src, 2, dst, 3, 3);
    print_dst(dst);
    wave2d_rotated_sample_to_sample_blit(src, 3, dst, 3, 3);
    print_dst(dst);
}

#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv){
    // INIT VARS
#define SAMPLE_WINDOW_W 3
#define SAMPLE_WINDOW_H 3
#define OUTPUT_W 30
#define OUTPUT_H 30
    
    // SETUP SAMPLES
    BEGIN_TIME(get_samples_time);
    
    Wave2D_Samples samples = {0};
    
    Wave2D_Image_Processing_Params img_proc_params = wave2d_image_processing_params(SRC_IMG_W, SRC_IMG_H, SAMPLE_WINDOW_W, SAMPLE_WINDOW_H, 1, 2, 1);
    void *img_proc_scratch_memory = malloc(img_proc_params.scratch_size);
    int32_t sample_count = wave2d_extract_samples_from_image(img_proc_params, &test_base_image[0][0], img_proc_scratch_memory);
    
    uint32_t *current_sample = wave2d_get_processed_samples_from_scratch(img_proc_params, img_proc_scratch_memory);
    uint8_t *current_count = wave2d_processed_counts_from_scratch(img_proc_params, img_proc_scratch_memory);
    
    int32_t sample_memory_size = sample_count*sizeof(Wave2D_Sample);
    void *sample_memory = malloc(sample_memory_size);
    wave2d_samples_memory(&samples, sample_memory, sample_memory_size);
    wave2d_begin_samples(&samples, SAMPLE_WINDOW_W, SAMPLE_WINDOW_H);
    for (int32_t i = 0; i < sample_count; ++i){
        wave2d_add_sample(&samples, current_sample, (float)(*current_count));
        current_sample += img_proc_params.next_sample_stride;
        ++current_count;
    }
    wave2d_end_samples(&samples);
    
    END_TIME(get_samples_time);
    
#if defined(DEBUGGING) && (DEBUGGING == 2)
    wave2d_print_samples(samples);
#endif
    
    // RUN THE WAVE COLLAPSE
    BEGIN_TIME(generator_init_time);
    
    Wave2D_State state = {0};
    uint32_t out[OUTPUT_W*OUTPUT_H];
    memset(out, 0, sizeof(out));
    
    int32_t state_memory_size = (1 << 19);
    void *state_memory = malloc(state_memory_size);
    wave2d_state_memory(&state, state_memory, state_memory_size);
    
    wave2d_initialize_state(&state, &samples, OUTPUT_W, OUTPUT_H);
    
    int32_t scratch_memory_size = (1 << 10)*64;
    void *scratch_memory = malloc(scratch_memory_size);
    
    END_TIME(generator_init_time);
    
    Random rng = {0};
    rng.inc = 17;
    
    int32_t print_result = 1;
    
#if defined(DEBUGGING) && (DEBUGGING == 4 || DEBUGGING == 5)
    print_result = 0;
    #endif
    
    for (int32_t i = 0; i < 10; ++i){
        rng.state = 12 + i;
        
        BEGIN_TIME(generator_run_time);
        int32_t result = wave2d_generate_output(&state, &rng, out, scratch_memory, scratch_memory_size);
        END_TIME(generator_run_time);
        
        // SHOW OUTPUT
        if (print_result){
        if (result){
            uint32_t *line = &out[0];
            for (int32_t y = 0; y < OUTPUT_H; ++y){
                uint32_t *src = line;
                for (int32_t x = 0; x < OUTPUT_W; ++x){
                    printf("%u", *src);
                    src += 1;
                }
                printf("\n");
                line += 30;
            }
            printf("\n");
        }
        else{
            printf("FAILED\n\n");
        }
    }
}

    return(0);
}


