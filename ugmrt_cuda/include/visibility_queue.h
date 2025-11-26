#ifndef VISIBILITY_QUEUE_H
#define VISIBILITY_QUEUE_H

#include <pthread.h>
#include <stdbool.h>
#include "cuda_types.h"  // For CudaVisibility struct

// Thread-safe visibility queue for parallel file processing
typedef struct {
    CudaVisibility* buffer;      // Circular buffer
    size_t capacity;             // Total capacity
    size_t head;                 // Write position
    size_t tail;                 // Read position
    size_t count;                // Current count

    pthread_mutex_t mutex;       // Protect queue state
    pthread_cond_t not_empty;    // Signal for consumers
    pthread_cond_t not_full;     // Signal for producers

    bool done;                   // All producers finished
    int active_producers;        // Count of active producer threads
} VisibilityQueue;

// Initialize queue with given capacity
int vis_queue_init(VisibilityQueue* queue, size_t capacity);

// Destroy queue and free resources
void vis_queue_destroy(VisibilityQueue* queue);

// Push single visibility (blocks if full)
void vis_queue_push(VisibilityQueue* queue, const CudaVisibility* vis);

// Pop batch of visibilities (blocks if empty, returns actual count)
size_t vis_queue_pop_batch(VisibilityQueue* queue, CudaVisibility* batch, size_t max_count);

// Signal that a producer thread has finished
void vis_queue_producer_done(VisibilityQueue* queue);

// Check if queue is finished (no more data will arrive)
bool vis_queue_is_finished(VisibilityQueue* queue);

// Get current queue size (for monitoring)
size_t vis_queue_size(VisibilityQueue* queue);

#endif // VISIBILITY_QUEUE_H
