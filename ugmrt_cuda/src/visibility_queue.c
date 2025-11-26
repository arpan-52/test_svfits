#include "visibility_queue.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int vis_queue_init(VisibilityQueue* queue, size_t capacity) {
    queue->buffer = (CudaVisibility*)malloc(capacity * sizeof(CudaVisibility));
    if (!queue->buffer) {
        fprintf(stderr, "Failed to allocate visibility queue buffer\n");
        return -1;
    }

    queue->capacity = capacity;
    queue->head = 0;
    queue->tail = 0;
    queue->count = 0;
    queue->done = false;
    queue->active_producers = 0;

    pthread_mutex_init(&queue->mutex, NULL);
    pthread_cond_init(&queue->not_empty, NULL);
    pthread_cond_init(&queue->not_full, NULL);

    return 0;
}

void vis_queue_destroy(VisibilityQueue* queue) {
    pthread_mutex_destroy(&queue->mutex);
    pthread_cond_destroy(&queue->not_empty);
    pthread_cond_destroy(&queue->not_full);
    free(queue->buffer);
    queue->buffer = NULL;
}

void vis_queue_push(VisibilityQueue* queue, const CudaVisibility* vis) {
    pthread_mutex_lock(&queue->mutex);

    // Wait while queue is full
    while (queue->count >= queue->capacity) {
        pthread_cond_wait(&queue->not_full, &queue->mutex);
    }

    // Add visibility to queue
    queue->buffer[queue->head] = *vis;
    queue->head = (queue->head + 1) % queue->capacity;
    queue->count++;

    // Signal consumer that data is available
    pthread_cond_signal(&queue->not_empty);

    pthread_mutex_unlock(&queue->mutex);
}

size_t vis_queue_pop_batch(VisibilityQueue* queue, CudaVisibility* batch, size_t max_count) {
    pthread_mutex_lock(&queue->mutex);

    // Wait while queue is empty AND not done
    while (queue->count == 0 && !queue->done) {
        pthread_cond_wait(&queue->not_empty, &queue->mutex);
    }

    // If queue is empty and done, return 0
    if (queue->count == 0 && queue->done) {
        pthread_mutex_unlock(&queue->mutex);
        return 0;
    }

    // Pop as many as possible (up to max_count)
    size_t to_pop = (queue->count < max_count) ? queue->count : max_count;

    for (size_t i = 0; i < to_pop; i++) {
        batch[i] = queue->buffer[queue->tail];
        queue->tail = (queue->tail + 1) % queue->capacity;
    }

    queue->count -= to_pop;

    // Signal producers that space is available
    pthread_cond_broadcast(&queue->not_full);

    pthread_mutex_unlock(&queue->mutex);

    return to_pop;
}

void vis_queue_producer_done(VisibilityQueue* queue) {
    pthread_mutex_lock(&queue->mutex);

    queue->active_producers--;

    // If all producers are done, mark queue as done
    if (queue->active_producers == 0) {
        queue->done = true;
        pthread_cond_broadcast(&queue->not_empty);  // Wake up consumer
    }

    pthread_mutex_unlock(&queue->mutex);
}

bool vis_queue_is_finished(VisibilityQueue* queue) {
    pthread_mutex_lock(&queue->mutex);
    bool finished = (queue->done && queue->count == 0);
    pthread_mutex_unlock(&queue->mutex);
    return finished;
}

size_t vis_queue_size(VisibilityQueue* queue) {
    pthread_mutex_lock(&queue->mutex);
    size_t size = queue->count;
    pthread_mutex_unlock(&queue->mutex);
    return size;
}
