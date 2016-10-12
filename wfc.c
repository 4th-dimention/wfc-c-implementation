/*
 * Experiments with Wafe Function Collapse algorithms for procedural generation.
 *  - Allen Webster
 * Created 02.10.2016
 *
 * Thanks to "mxgmn" and the github page https://github.com/mxgmn/WaveFunctionCollapse
 *  for providing the initial concept.
 *
 */


//#define DEBUGGING 1


#include <stdint.h>
#include <math.h>
#include <string.h>

#define ArrayCount(a) (sizeof(a)/sizeof(*(a)))

#if defined(DEBUGGING)
#include <stdio.h>
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

uint32_t test_samples[][2][2] = {
    {{0, 1}, {1, 0}},
    {{1, 0}, {0, 1}},
    
    {{0, 0}, {1, 1}},
    {{1, 1}, {0, 0}},
    
    {{1, 0}, {1, 0}},
    {{0, 1}, {0, 1}},
    
    {{1, 1}, {1, 0}},
    {{1, 1}, {0, 1}},
    {{1, 0}, {1, 1}},
    {{0, 1}, {1, 1}},
    
    {{0, 0}, {0, 1}},
    {{0, 0}, {1, 0}},
    {{0, 1}, {0, 0}},
    {{1, 0}, {0, 0}},
    
    {{0, 0}, {0, 0}},
    
    {{1, 1}, {1, 1}},
    
    {{1, 1}, {1, 2}},
    {{1, 1}, {2, 1}},
    {{1, 2}, {1, 1}},
    {{2, 1}, {1, 1}},
    
    {{2, 0}, {1, 1}},
    {{1, 1}, {2, 0}},
    {{0, 2}, {1, 1}},
    {{1, 1}, {0, 2}},
    
    {{1, 2}, {1, 0}},
    {{2, 1}, {0, 1}},
    {{1, 0}, {1, 2}},
    {{0, 1}, {2, 1}},
    
    {{0, 2}, {1, 0}},
    {{2, 0}, {0, 1}},
    {{0, 1}, {2, 0}},
    {{1, 0}, {0, 2}},
    
    {{3, 0}, {0, 0}},
    {{0, 3}, {0, 0}},
    {{0, 0}, {3, 0}},
    {{0, 0}, {0, 3}},
};

int32_t occurences[] = {
    // DIAGONALS
    1,1,
    
    // ROWS
    10,10,
    
    // COLUMNS
    10,10,
    
    // LS
    5,5,5,5,
    
    // CORNERS
    5,5,5,5,
    
    // BIG OPEN SPACE
    25,
    
    // BIG THICK WALL
    1,
    
    // TREASURE
    1,1,1,1,
    1,1,1,1,
    1,1,1,1,
    1,1,1,1,
    
    // STATUE
    1,1,1,1,
};

#if 0
uint32_t test_samples[][3][3] = {
    // OPEN SPACE
    {
        {0,0,0},
        {0,0,0},
        {0,0,0},
    },
    
    // ROWS
    {
        {1,1,1},
        {0,0,0},
        {0,0,0},
    },
    
    {
        {0,0,0},
        {1,1,1},
        {0,0,0},
    },
    
    {
        {0,0,0},
        {0,0,0},
        {1,1,1},
    },
    
    {
        {1,1,1},
        {0,0,0},
        {1,1,1},
    },
    
    // COLUMNS
    {
        {1,0,0},
        {1,0,0},
        {1,0,0},
    },
    
    {
        {0,1,0},
        {0,1,0},
        {0,1,0},
    },
    
    {
        {0,0,1},
        {0,0,1},
        {0,0,1},
    },
    
    {
        {1,0,1},
        {1,0,1},
        {1,0,1},
    },
    
    // BENDS
    {
        {1,1,1},
        {0,0,1},
        {0,0,1},
    },
    
    {
        {0,0,1},
        {0,0,1},
        {1,1,1},
    },
    
    {
        {1,0,0},
        {1,0,0},
        {1,1,1},
    },
    
    {
        {1,1,1},
        {1,0,0},
        {1,0,0},
    },
    
    {
        {0,0,0},
        {1,1,1},
        {0,0,1},
    },
    
    {
        {0,0,1},
        {1,1,1},
        {0,0,0},
    },
    
    {
        {1,0,0},
        {1,1,1},
        {0,0,0},
    },
    
    {
        {0,0,0},
        {1,1,1},
        {1,0,0},
    },
    
    {
        {1,1,0},
        {0,1,0},
        {0,1,0},
    },
    
    {
        {0,1,0},
        {0,1,0},
        {1,1,0},
    },
    
    {
        {0,1,0},
        {0,1,0},
        {0,1,1},
    },
    
    {
        {0,1,1},
        {0,1,0},
        {0,1,0},
    },
    
    {
        {0,0,0},
        {1,1,0},
        {0,1,0},
    },
    
    {
        {0,1,0},
        {1,1,0},
        {0,0,0},
    },
    
    {
        {0,1,0},
        {0,1,1},
        {0,0,0},
    },
    
    {
        {0,0,0},
        {0,1,1},
        {0,1,0},
    },
    
    {
        {1,1,1},
        {0,0,1},
        {1,0,1},
    },
    
    {
        {1,0,1},
        {0,0,1},
        {1,1,1},
    },
    
    {
        {1,0,1},
        {1,0,0},
        {1,1,1},
    },
    
    {
        {1,1,1},
        {1,0,0},
        {1,0,1},
    },
    
    // CORNERS
    {
        {0,0,0},
        {0,0,0},
        {1,0,0},
    },
    
    {
        {0,0,0},
        {0,0,0},
        {0,0,1},
    },
    
    {
        {0,0,1},
        {0,0,0},
        {0,0,0},
    },
    
    {
        {1,0,0},
        {0,0,0},
        {0,0,0},
    },
    
    // Ts
    {
        {0,0,1},
        {1,1,1},
        {0,0,1},
    },
    
    {
        {1,0,0},
        {1,1,1},
        {1,0,0},
    },
    
    {
        {1,1,1},
        {0,1,0},
        {0,1,0},
    },
    
    {
        {0,1,0},
        {0,1,0},
        {1,1,1},
    },
    
    {
        {0,0,0},
        {1,1,1},
        {0,1,0},
    },
    
    {
        {0,1,0},
        {1,1,0},
        {0,1,0},
    },
    
    {
        {0,1,0},
        {1,1,1},
        {0,0,0},
    },
    
    {
        {0,0,0},
        {1,1,1},
        {0,1,0},
    },
    
    {
        {0,1,0},
        {1,1,1},
        {0,1,0},
    },
};

int32_t occurences[] = {
    // OPEN SPACE
    1,
    
    // ROWS
    1,1,1,1,
    
    // COLUMNS
    1,1,1,1,
    
    // BENDS
    1,1,1,1,
    1,1,1,1,
    1,1,1,1,
    1,1,1,1,
    1,1,1,1,
    
    // CORNERS
    1,1,1,1,
    
    // Ts
    1,1,1,1,
    1,1,1,1,
    1,
};
#endif

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
    
    uint32_t w, h;
    uint32_t sample_count;
    uint32_t coefficient_count;
    uint8_t *coefficients;
    
    float log_T;
    float *weight_by_log_weight;
} Wave2D_State;

static void
wave2d_state_memory(Wave2D_State *state, void *memory, int32_t size){
    state->part = part_make(memory, size);
}

static void
wave2d_set_output_size(Wave2D_State *state, Wave2D_Samples *samples, uint32_t w, uint32_t h){
    uint32_t grid_cell_count = w*h;
    uint32_t sample_count = samples->sample_count;
    uint32_t coefficient_count = grid_cell_count * sample_count;
    
    state->w = w;
    state->h = h;
    state->sample_count = sample_count;
    state->coefficient_count = coefficient_count;
    
    part_clear(&state->part);
    state->coefficients = push_array(&state->part, uint8_t, coefficient_count);
    
    state->log_T = logf((float)sample_count);
    state->weight_by_log_weight = push_array(&state->part, float, sample_count);
    
    float *weight_by_log_weight = state->weight_by_log_weight;
    float *weight_ptr = &samples->samples[0].weight;
    int32_t stride = sizeof(samples->samples[0]);
    for (uint32_t i = 0; i < sample_count; ++i){
        weight_by_log_weight[i] = (*weight_ptr) * logf(*weight_ptr);
        weight_ptr = (float*)((uint8_t*)weight_ptr + stride);
    }
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
            on_change_queue[i] = 1;
        }
    }
    else{
        result = AddChange_Success_Already_On_Queue;
    }
    return(result);
}

static int32_t
wave2d_generate_output(Wave2D_State *state, Wave2D_Samples *samples, Random *rng, uint32_t *out, void *scratch, int32_t scratch_size){
    int32_t result = 0;
    
    uint32_t coefficient_count = state->coefficient_count;
    
    uint32_t w = state->w;
    uint32_t h = state->h;
    uint32_t cell_count = w*h;
    uint32_t sample_count = state->sample_count;
    
    uint32_t sample_w = samples->w;
    uint32_t sample_h = samples->h;
    
    uint8_t *coefficients = state->coefficients;
    Wave2D_Sample *sample_array = samples->samples;
    float *weight_by_log_weight = state->weight_by_log_weight;
    
    uint8_t *scratch_start = scratch;
    uint8_t *scratch_ptr = scratch_start;
    uint8_t *scratch_end = scratch_ptr + scratch_size;
    
    uint8_t *on_change_queue = (uint8_t*)scratch_ptr;
    memset(on_change_queue, 0, cell_count);
    scratch_ptr += cell_count;
    
    uint8_t *cell_finished = (uint8_t*)scratch_ptr;
    memset(cell_finished, 0, cell_count);
    scratch_ptr += cell_count;
    
    int32_t pos = (int32_t)(scratch_ptr - scratch_start);
    pos = (pos+3)&(~3);
    scratch_ptr = (uint8_t*)(scratch) + pos;
    
    Wave2D_Change *changes = (Wave2D_Change*)scratch_ptr;
    int32_t change_max = (int32_t)(scratch_end - scratch_ptr)/sizeof(Wave2D_Change);
    int32_t change_count = 0;
    
    // CLEAR
    for (uint32_t i = 0; i < coefficient_count; ++i){
        coefficients[i] = 1;
    }
    
    // COLLAPSE
    uint32_t failed = 0;
    for (int32_t step = 0; ; ++step){
        // OBSERVE
        //  find smallest entropy
        
        uint32_t best_x = w;
        uint32_t best_y = h;
        float best_total_weight = 0;
        float min_entropy = 100000.f;
        
        for (uint32_t x = 0; x < w; ++x){
            for (uint32_t y = 0; y < h; ++y){
                
                if (!cell_finished[(x + y*w)]){
                    uint8_t *cell_coefficients = &coefficients[(x + y*w) * sample_count];
                    
                    uint32_t valid_sample_count = 0;
                    float total_weight = 0.f;
                    float entropy = 0;
                    
                    float *weight_ptr = &sample_array[0].weight;
                    int32_t stride = sizeof(sample_array[0]);
                    for (uint32_t t = 0; t < sample_count; ++t){
                        if (cell_coefficients[t]){
                            ++valid_sample_count;
                            total_weight += *weight_ptr;
                        }
                        weight_ptr = (float*)((uint8_t*)weight_ptr + stride);
                    }
                    
                    // TODO(allen): Is it really faster this way?
                    // I guess the point is to avoid too many array lookups
                    // in common short-circuitable cases?  Is this actually
                    // speeding it up for us?
                    if (valid_sample_count == 0){
                        // TODO(allen): This can be an assert once we start using that.
                        failed = 1;
                        goto finished;
                    }
                    else if (valid_sample_count == sample_count){
                        entropy = state->log_T;
                    }
                    else{
                        float entropy_sum = 0;
                        float log_total_weight = logf(total_weight);
                        for (uint32_t t = 0; t < sample_count; ++t){
                            if (cell_coefficients[t]){
                                entropy_sum += weight_by_log_weight[t];
                            }
                        }
                        entropy = log_total_weight - entropy_sum/total_weight;
                    }
                    
                    float noise = 0.00001f*pcg32_random_float_01(rng);
                    entropy += noise;
                    if (entropy > 0){
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
        
        if (best_x == w && best_y == h){
            goto finished;
        }
        
        cell_finished[(best_x + best_y*w)] = 1;
        
        uint8_t *cell_coefficients = &coefficients[(best_x + best_y*w) * sample_count];
        
        float die_roll = best_total_weight * pcg32_random_float_01(rng);
        
        float *weight_ptr = &sample_array[0].weight;
        int32_t stride = sizeof(sample_array[0]);
        uint32_t hit_element = sample_count;
        for (uint32_t t = 0; t < sample_count; ++t){
            die_roll -= (*weight_ptr)*cell_coefficients[t];
            if (die_roll < 0){
                hit_element = t;
                break;
            }
            
            weight_ptr = (float*)((uint8_t*)weight_ptr + stride);
        }
        
        if (hit_element == sample_count){
            failed = 1;
            goto finished;
        }
        
        for (uint32_t t = 0; t < sample_count; ++t){
            cell_coefficients[t] = (t == hit_element);
        }
        
        {
            Wave2D_Change new_change;
            new_change.x = best_x;
            new_change.y = best_y;
            int32_t added_change = wave2d_add_change(new_change, changes, &change_count, change_max, on_change_queue, w);
            
            if (added_change != AddChange_Success){
                failed = 1;
                goto finished;
            }
        }
        
        // PROPOGATE
        for (int32_t i = 0; i < change_count; ++i){
            Wave2D_Change change = changes[i];
            
            int32_t delta_x = 0, delta_y = 0;
            
            for (delta_x = -(int32_t)(sample_w) + 1; delta_x < (int32_t)(sample_w); ++delta_x){
                for (delta_y = -(int32_t)(sample_h) + 1; delta_y < (int32_t)(sample_h); ++delta_y){
                    uint32_t x2 = (uint32_t)((change.x + delta_x + w) % w);
                    uint32_t y2 = (uint32_t)((change.y + delta_y + h) % h);
                    
                    uint8_t *cell1 = &coefficients[(change.x + change.y*w)*sample_count];
                    uint8_t *cell2 = &coefficients[(x2 + y2*w)*sample_count];
                    
                    uint32_t cell_2_changed = 0;
                    
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
                    
                    uint8_t cell2_is_finished = cell_finished[x2 + y2*w];
                    int32_t potential_state_count = 0;
                    
                    for (uint32_t t2 = 0; t2 < sample_count; ++t2){
                        if (cell2[t2]){
                            uint32_t can_coexist = 0;
                            for (uint32_t t1 = 0; t1 < sample_count; ++t1){
                                if (cell1[t1]){
                                    uint32_t coexist_check_result = 1;
                                    
                                    // COEXIST CHECK
                                    uint32_t *sample_1 = samples->samples[t1].data;
                                    uint32_t *sample_2 = samples->samples[t2].data;
                                    
                                    for (int32_t sx = minx; sx < maxx; ++sx){
                                        for (int32_t sy = miny; sy < maxy; ++sy){
                                            if (sample_1[sx + sy*sample_w] != sample_2[(sx-delta_x) + (sy-delta_y)*sample_h]){
                                                coexist_check_result = 0;
                                                goto finish_coexist_check;
                                            }
                                        }
                                    }
                                    finish_coexist_check:;
                                    
                                    if (coexist_check_result){
                                        can_coexist = 1;
                                        break;
                                    }
                                }
                            }
                            
                            if (!can_coexist){
                                if (cell2_is_finished){
                                    failed = 1;
                                    goto finished;
                                }
                                
                                // TODO(allen): check for a newly finished cell / contradictions
                                cell2[t2] = 0;
                                cell_2_changed = 1;
                            }
                            else{
                                ++potential_state_count;
                            }
                        }
                    }
                    
                    if (!cell2_is_finished && potential_state_count == 1){
                        cell_finished[x2 + y2*w] = 1;
                    }
                    
                    if (cell_2_changed){
                        if (change_count == change_max){
                            int32_t n = i;
                            memmove(changes, changes + n, change_count - n);
                            change_count -= n;
                            i = 0;
                        }
                        
                        Wave2D_Change new_change;
                        new_change.x = x2;
                        new_change.y = y2;
                        
                        int32_t added_change = wave2d_add_change(new_change, changes, &change_count, change_max, on_change_queue, w);
                        if (added_change == AddChange_Failed){
                            failed = 1;
                            goto finished;
                        }
                    }
                }
            }
            
            on_change_queue[change.x + change.y*w] = 0;
        }
        
        change_count = 0;
        
#if defined(DEBUGGING)
        for (uint32_t y = 0; y < h; ++y){
            for (uint32_t x = 0; x < w; ++x){
                if (cell_finished[x + y*w]){
                    uint8_t *cell = &coefficients[(x + y*w)*sample_count];
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
    }
    
    finished:;
    if (!failed){
        for (uint32_t x = 0; x < w; ++x){
            for (uint32_t y = 0; y < h; ++y){
                uint8_t *cell = &coefficients[(x + y*w)*sample_count];
                uint32_t t = 0;
                for (; t < sample_count; ++t){
                    if (cell[t]) break;
                }
                
                uint32_t *sample = samples->samples[t].data;
                out[x + y*w] = sample[0];
            }
        }
        result = 1;
    }
    
    return(result);
}

#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv){
    
    // SETUP SAMPLES
    Wave2D_Samples samples = {0};
    
    int32_t sample_memory_size = (1 << 10)*64;
    void *sample_memory = malloc(sample_memory_size);
    
    wave2d_samples_memory(&samples, sample_memory, sample_memory_size);
    
    wave2d_begin_samples(&samples, 2, 2);
    for (int32_t i = 0; i < ArrayCount(test_samples); ++i){
        wave2d_add_sample(&samples, &(test_samples[i][0][0]), (float)occurences[i]);
    }
    wave2d_end_samples(&samples);
    
    // RUN THE WAVE COLLAPSE
    Wave2D_State state = {0};
    uint32_t out[30*30];
    memset(out, 0, sizeof(out));
    
    int32_t state_memory_size = (1 << 10)*64;
    void *state_memory = malloc(state_memory_size);
    wave2d_state_memory(&state, state_memory, state_memory_size);
    
    wave2d_set_output_size(&state, &samples, 30, 30);
    
    int32_t scratch_memory_size = (1 << 10)*64;
    void *scratch_memory = malloc(scratch_memory_size);
    
    Random rng = {0};
    rng.inc = 17;
    
    for (int32_t i = 0; i < 10; ++i){
        rng.state = 12 + i;
        
        int32_t result = wave2d_generate_output(&state, &samples, &rng, out, scratch_memory, scratch_memory_size);
        
        // SHOW OUTPUT
        if (result){
            uint32_t *line = &out[0];
            for (int32_t y = 0; y < 30; ++y){
                uint32_t *src = line;
                for (int32_t x = 0; x < 30; ++x){
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
    
    return(0);
}

