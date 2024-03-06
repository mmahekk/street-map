#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"
#include <float.h>


// Node structure
struct node {
    int id;             // Node ID
    double lat;
    double lon;
    int num_ways;       // Number of ways associated with this node
    int *way_ids;       // Array of way IDs associated with this node
};


// Way structure
struct way {
    int id;             // Way ID
    char *name;         // Name of the road segment
    float maxspeed;     // Maximum speed limit of road segment, in km/hr
    bool oneway;        // True if the road segment is one-way only, false otherwise
    int num_nodes;      // Number of nodes associated with this way object
    int *node_ids;      // Array of node IDs associated with this way object
};


// Main map structure
struct ssmap {
    int nr_nodes;       // Number of nodes
    int nr_ways;        // Number of ways
    struct node *nodes; // Array of nodes
    struct way *ways;   // Array of ways
};


struct ssmap * 
ssmap_create(int nr_nodes, int nr_ways)
{
    // If either parameter is zero, return NULL as there is no map to create
    if (nr_nodes <= 0 || nr_ways <= 0) {
        return NULL;
    }

    struct ssmap *map = malloc(sizeof(struct ssmap));
    if (map == NULL) {  // If memory allocation fails, return NULL
        return NULL;    
    }

    map->nodes = malloc(nr_nodes * sizeof(struct node));
    if (map->nodes == NULL) {   // If memory allocation fails, free the map and return NULL
        free(map);
        return NULL;
    }

    map->ways = malloc(nr_ways * sizeof(struct way));
    if (map->ways == NULL) { // If memory allocation fails, free the nodes and map 
        free(map->nodes);
        free(map);
        return NULL;
    }

    // If both allocations succeeded, initialize other members as needed
    map->nr_nodes = nr_nodes;
    map->nr_ways = nr_ways;
    return map;
}


bool
ssmap_initialize(struct ssmap * m)
{
    return true; // No other items to initialize
}


void
ssmap_destroy(struct ssmap * m)
{
    // Free each dynamically allocated name in ways and node_ids arrays from ssmap_add_way.
    for (int i = 0; i < m->nr_ways; i++) {
        free(m->ways[i].name);
        free(m->ways[i].node_ids);
    }

    // Free the way_ids arrays for each node from ssmap_add_nodes.
    for (int j = 0; j < m->nr_nodes; j++) {
        free(m->nodes[j].way_ids);
    }

    // Free the nodes and ways arrays themselves from ssmap_create.
    free(m->nodes);
    free(m->ways);

    // Free the ssmap structure.
    free(m);
}


struct way * 
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed, bool oneway, 
              int num_nodes, const int node_ids[num_nodes])
{
    struct way *w = &(m->ways[id]);  
    
    // Copy the name with strdup
    w->name = strdup(name);
    if (w->name == NULL) {
        // Handle strdup failure due to memory allocation error.
        return NULL;
    }
  
    // Set the rest of the fields.
    w->id = id;
    w->maxspeed = maxspeed;
    w->oneway = oneway;
    w->num_nodes = num_nodes;

    // Allocate and copy node_ids, handling potential malloc failure.
    w->node_ids = malloc(num_nodes * sizeof(int));
    if (w->node_ids == NULL) {
        // Cleanup due to memory allocation error.
        free(w->name); // Free the strdup allocated memory.
        return NULL;
    }
    memcpy(w->node_ids, node_ids, num_nodes * sizeof(int));

    return w; // Successful creation and initialization.
}


struct node * 
ssmap_add_node(struct ssmap * m, int id, double lat, double lon, 
               int num_ways, const int way_ids[num_ways])
{
    struct node *n = &(m->nodes[id]);    
    
    // Set the node's fields
    n->id = id;
    n->lat = lat;
    n->lon = lon;
    n->num_ways = num_ways;

    // Allocate memory for the way_ids array and copy the way IDs
    n->way_ids = malloc(num_ways * sizeof(int));
    if (n->way_ids == NULL) {
        // Handle malloc failure due to memory allocation error
        return NULL;
    }
    memcpy(n->way_ids, way_ids, num_ways * sizeof(int));

    // Return a pointer to the newly added node
    return n;
}


void
ssmap_print_way(const struct ssmap * m, int id)
{
    // Check if the provided id is within the valid range
    if (id < 0 || id >= m->nr_ways) {
        printf("error: way %d does not exist.\n", id);
        return; // Exit the function if the id is invalid
    }

    struct way *w = &(m->ways[id]);
    printf("Way %d: %s\n", id, w->name);
}


void
ssmap_print_node(const struct ssmap * m, int id)
{
    // Check if the provided id is within the valid range
    if (id < 0 || id >= m->nr_nodes) {
        printf("error: node %d does not exist.\n", id);
        return; // Exit the function if the id is invalid
    }

    struct node *n = &(m->nodes[id]);
    printf("Node %d: (%.7lf, %.7lf)\n", id, n->lat, n->lon);
}


void 
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{
    bool found = false; // Variable to track if any way is found

    for (int i = 0; i < m->nr_ways; i++) {
        struct way *w = &(m->ways[i]);
        // Check if the way's name contains the keyword
        if (strstr(w->name, name) != NULL) {
            // Print the way id followed by a space
            printf("%d ", w->id);
            found = true; // Set found flag to true
        }
    }

    // If at least one way is found, print a newline at the end
    if (found) {
        printf("\n");
    } else {
        // If no way is found that matches the keyword, print newline
        printf("\n");
    }
}


/**
 * Counts the number of ways within a given map structure that match one or both of 
 * the specified names.
 * 
 * This function iterates through all ways in the provided map structure and checks 
 * each way's name for the presence of `name1` and optionally `name2`. A way is 
 * counted if it contains either `name1` or `name2`.
 * If a way is found that contains both `name1` and `name2`, and it is the only 
 * way found, the function returns -1 to indicate this special case. Otherwise, 
 * the function returns the count of matching ways.
 *
 * @param m A pointer to the `ssmap` structure containing the ways to be searched.
 * @param name1 The first name to search for within the way names. 
 * @param name2 The second name to search for within the way names. 
 * This parameter is optional and can be NULL.
 * If `name2` is provided, a way containing both `name1` and `name2` is considered a special case.
 *
 * @return The number of ways that match `name1` or `name2`. If a way is found that 
 * contains both `name1` and `name2`, and it is the only way found, returns -1 to indicate 
 * this special case. Otherwise, returns the count of matching ways. If no ways match, returns 0.
 */
int count_matching_ways(const struct ssmap *m, const char *name1, const char *name2) {
    int count = 0; // Initialize counter for matching ways
    bool found_combined = false; // Flag to indicate if a way with both names was found

    for (int i = 0; i < m->nr_ways; i++) {
        bool found_name1 = strstr(m->ways[i].name, name1) != NULL;
        bool found_name2 = name2 && strstr(m->ways[i].name, name2) != NULL;

        if (found_name1 || found_name2) {
            count++;
            if (found_name1 && found_name2) {
                // This way contains both names
                found_combined = true;
            }
        }
    }

    // If a way is found with both names and it's the only way found, return -1 to 
    // indicate this special case
    if (found_combined && count == 1) {
        return -1;
    } else {
        return count;
    }
}


void ssmap_find_node_by_names(const struct ssmap *m, const char *name1, const char *name2) {
    // Check for the case where both name1 and name2 are part of the same way
    int matching_ways = count_matching_ways(m, name1, name2);
    if (matching_ways == -1) {
        printf("\n");
        return;
    }

    // Iterate through all nodes in the map
    for (int i = 0; i < m->nr_nodes; i++) {
        bool found_name1 = false, found_name2 = false;
        int distinct_ways_count = 0;

        // Check each way associated with the current node for name1 and name2
        for (int j = 0; j < m->nodes[i].num_ways; j++) {
            struct way *way = &m->ways[m->nodes[i].way_ids[j]];

            // If name1 matches, mark found_name1 as true
            if ((strstr(way->name, name1) != NULL) && !found_name1) {
                found_name1 = true;
                distinct_ways_count++; // Count this as a distinct way
            } else if ((strstr(way->name, name1) != NULL) && (name2 == NULL || 
            (name2 != NULL && strcmp(name1, name2) == 0))) {
                // If name1 and name2 are the same, and we already found name1 before,
                // increment distinct_ways_count only if this is a new match in a different way
                distinct_ways_count++;
            }
        
            // If name2 is not NULL and different from name1, and name2 matches, 
            // mark found_name2 as true
            if (name2 != NULL && strcmp(name1, name2) != 0 && strstr(way->name, name2) != NULL) {
                found_name2 = true;
                distinct_ways_count++; // Count this as a distinct way
            }
        }

        // Criteria for printing the node ID:
        // - If name1 and name2 are the same, we need at least two distinct ways with name1
        // - If name1 and name2 are different, the node must be connected to both ways
        if ((name2 == NULL && found_name1) || 
            (name2 != NULL && strcmp(name1, name2) == 0 && distinct_ways_count >= 2) || 
            (name2 != NULL && strcmp(name1, name2) != 0 && found_name1 && found_name2)) {
            printf("%d ", m->nodes[i].id);
        }
    }

    printf("\n");
}


/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)


/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) {
    double R = 6371.;       
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1); 
    double dlon = d2r(lon2-lon1); 
    double a = pow(sin(dlat/2), 2) + cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    return R * c; 
}


/**
 * Checks if there is a direct connection between two nodes in a map.
 * 
 * A direct connection exists if both nodes are found within the same way. The function 
 * returns `true` if such a connection is found, indicating that the nodes are directly 
 * connected through at least one way in the map. If no direct connection is found, the 
 * function prints an error message and returns `false`.
 *
 * @param m A pointer to the `ssmap` structure representing the map to search within.
 * @param node_a The ID of the first node to check for a direct connection.
 * @param node_b The ID of the second node to check for a direct connection.
 * 
 * @return Returns `true` if a direct connection between node_a and node_b exists within any
 *         of the ways in the map. Otherwise, prints an error message and returns `false`.
 */
bool has_direct_connection(const struct ssmap *m, int node_a, int node_b) {
    // Iterate through all ways to check for a direct connection between node_a and node_b
    for (int i = 0; i < m->nr_ways; i++) {
        struct way *way = &m->ways[i];
        bool found_a = false, found_b = false;

        // Check if both nodes are part of the current way
        for (int j = 0; j < way->num_nodes; j++) {
            if (way->node_ids[j] == node_a) found_a = true;
            if (way->node_ids[j] == node_b) found_b = true;
        }

        // If both nodes are found within the same way, a direct connection exists
        if (found_a && found_b) {
            return true;
        }
    }

    // If no direct connection is found after checking all ways, print the error message
    printf("error: there are no roads between node %d and node %d.\n", node_a, node_b);
    return false; 
}


/**
 * Checks if two nodes are adjacent in any way within a map.
 * Iterates through all ways to find if node_a and node_b are directly next to each other or. 
 * Prints an error for non-adjacent nodes and returns `false`; otherwise, returns 
 * `true` if adjacent.
 *
 * @param m Pointer to the map structure (`ssmap`).
 * @param node_a ID of the first node.
 * @param node_b ID of the second node.
 * 
 * @return `true` if nodes are adjacent; otherwise, `false`.
 */
bool are_nodes_adjacent_in_way(const struct ssmap *m, int node_a, int node_b) {
    for (int i = 0; i < m->nr_ways; i++) {
        struct way *way = &m->ways[i];
        bool is_circular_path = (way->node_ids[0] == way->node_ids[way->num_nodes - 1]);

        // Initialize positions as not found
        int posA = -1, posB = -1;

        // Search for the positions of node_a and node_b in the current way
        for (int j = 0; j < way->num_nodes; j++) {
            if (way->node_ids[j] == node_a) posA = j;
            if (way->node_ids[j] == node_b) posB = j;
        }

        // Ensure both nodes are found in the current way
        if ((posA != -1 && posB != -1) && is_circular_path) {
            // For circular paths, check if nodes are adjacent or are start and end nodes
            if (abs(posA - posB) == 1 || abs(posA - posB) == way->num_nodes - 2) {
                return true;
            } else {
                return false;
            }
        } else if ((posA != -1 && posB != -1) && !is_circular_path) {
            // For non-circular paths, simply check if nodes are adjacent
            if (abs(posA - posB) == 1) {
                return true;
            } else {
                printf("error: cannot go directly from node %d to node %d.\n", node_a, node_b);
                return false;
            }
        }
    }

    // If we reach this point, the nodes are not adjacent in any way
    return false;
}


/**
 * Verifies if travel from node_a to node_b follows the correct direction on one-way paths.
 * If node_b follows node_a in a one-way, returns `true`. If node_a follows node_b, 
 * prints an error and returns `false`, indicating reverse travel on a one-way. 
 * If no one-way or no direct connection, returns `true`.
 *
 * @param m Pointer to the map structure (`ssmap`).
 * @param node_a ID of the start node.
 * @param node_b ID of the end node.
 * @return `true` if direction is correct or no one-way restriction; otherwise, `false`.
 */
bool is_travel_direction_correct_on_one_way(const struct ssmap *m, int node_a, int node_b) {
    for (int i = 0; i < m->nr_ways; i++) {
        struct way *way = &m->ways[i];
        // Only proceed if the way is a one-way
        if (!way->oneway) continue;

        for (int j = 0; j < way->num_nodes - 1; j++) {
            if (way->node_ids[j] == node_a && way->node_ids[j + 1] == node_b) {
                return true; // Correct direction on one-way
            }
            if (way->node_ids[j] == node_b && way->node_ids[j + 1] == node_a) {
                printf("error: cannot go in reverse from node %d to node %d.\n", node_a, node_b);
                return false; // Incorrect direction on one-way
            }
        }
    }

    // No one-way restriction found or direction is correct
    return true;
}


/**
 * Calculates the travel time between two connected nodes in a map.
 *
 * This function assumes node_id_a and node_id_b are connected and finds the 
 * travel time based on the distance between them and the max speed of the connecting way. 
 * The result is stored in `travel_time`.
 *
 * @param m Pointer to the map structure (`ssmap`).
 * @param node_id_a ID of the first node.
 * @param node_id_b ID of the second node.
 * @param travel_time Pointer to store the calculated travel time in minutes.
 * @return `true` if a connecting way is found and travel time is calculated, otherwise `false`.
 */
bool calculate_distance(const struct ssmap *m, int node_id_a, int node_id_b, double *travel_time) {
    // Assume that node_id_a and node_id_b are connected by a way and in the correct order.
    const struct node *node_a = &m->nodes[node_id_a];
    const struct node *node_b = &m->nodes[node_id_b];

    // Calculate the distance between the nodes
    double distance = distance_between_nodes(node_a, node_b);

    // Find the way that connects these nodes to determine the max speed.
    for (int i = 0; i < m->nr_ways; i++) {
        struct way *way = &m->ways[i];
        for (int j = 0; j < way->num_nodes - 1; j++) {
            if ((way->node_ids[j] == node_id_a && way->node_ids[j + 1] == node_id_b) ||
                (way->node_ids[j] == node_id_b && way->node_ids[j + 1] == node_id_a)) {
                // We found the way connecting node_id_a and node_id_b.
                // Now calculate travel time using the maximum speed of this way.
                *travel_time = (distance / way->maxspeed) * 60.0; // Convert from hours to minutes
                return true;   
            }
        }
    }

    return false; // If no way is found connecting the nodes, return false.
}


/**
 * Calculates the total travel time for a path through specified nodes in a map.
 *
 * Validates the path by ensuring all node IDs exist, are unique, and each pair of 
 * consecutive nodes is directly connected, adjacent in a way, and follows correct 
 * travel direction on one-way paths. The travel time is computed based on distances 
 * and maximum speeds along the path. Errors for invalid paths (non-existent nodes, 
 * duplicates, lack of direct connection, non-adjacency, or incorrect direction
 * on one-ways) are printed, and the function returns -1.0 in such cases.
 *
 * @param m Pointer to the map structure (`ssmap`).
 * @param size Number of nodes in the path.
 * @param node_ids Array of node IDs representing the path.
 * @return Total travel time for the path in minutes, or -1.0 if the path is invalid.
 */
double ssmap_path_travel_time(const struct ssmap *m, int size, int node_ids[size]) {
    if (size == 1) {
        return 0.0; // Single node path has no travel time
    }
    
    double total_time = 0.0;

    // Validate each node ID and check for duplicates
    for (int i = 0; i < size; i++) {
        if (node_ids[i] < 0 || node_ids[i] >= m->nr_nodes) {
            printf("error: node %d does not exist.\n", node_ids[i]);
            return -1.0;
        }

        for (int j = i + 1; j < size; j++) {
            if (node_ids[i] == node_ids[j]) {
                printf("error: node %d appeared more than once.\n", node_ids[i]);
                return -1.0;
            }
        }
    }

    // Iterate through path to calculate total travel time
    for (int i = 0; i < size - 1; i++) {
        if (!has_direct_connection(m, node_ids[i], node_ids[i + 1])) {
            // If there's no direct connection, has_direct_connection function will print the error.
            return -1.0;
        }
        if (!are_nodes_adjacent_in_way(m, node_ids[i], node_ids[i + 1])) {
            // Error message for non-adjacency is handled within the function.
            return -1.0;
        }
        if (!is_travel_direction_correct_on_one_way(m, node_ids[i], node_ids[i + 1])) {
            // Error message for incorrect direction on one-way is handled within the function.
            return -1.0;
        }
        
       double segment_time = 0.0;
       calculate_distance(m, node_ids[i], node_ids[i + 1], &segment_time);
       total_time += segment_time;
   }

   return total_time;
}


// Node in the priority queue/min heap
struct pq_node {
    int node_id;    // Node ID
    double cost;    // Cost to reach this node from the start node
};


// Min heap structure
struct min_heap {
    struct pq_node *elements; // Array of heap nodes
    int *positions;           // Tracks the position of nodes within the heap
    int size;                 // Current size of the heap
    int capacity;             // Maximum capacity of the heap
};


/**
 * Initializes a minimum heap with a specified capacity.
 *
 * Initializes all positions to -1, indicating that nodes are not in the heap. 
 *
 * @param capacity The maximum number of elements the heap can hold.
 * @return Pointer to the initialized min heap structure, or NULL if memory allocation fails.
 */
struct min_heap* init_min_heap(int capacity) {
    struct min_heap *heap = (struct min_heap *)malloc(sizeof(struct min_heap)); 
    if (!heap) return NULL; // Return NULL if memory allocation fails

    // Allocate memory for the heap's elements array and positions array
    heap->elements = (struct pq_node *)malloc(sizeof(struct pq_node) * capacity); 
    heap->positions = (int *)malloc(sizeof(int) * capacity);
    if (!heap->elements || !heap->positions) {
        free(heap->elements); 
        free(heap->positions);
        free(heap);
        return NULL;
    }

    heap->size = 0;
    heap->capacity = capacity;
    for (int i = 0; i < capacity; i++) {
        heap->positions[i] = -1; // Indicates that the node is not in the heap
    }
    return heap;
}


/**
 * Performs the heapify-up operation on a min heap.
 *
 * This function corrects the heap property starting from a given index up to the root of the heap,
 * ensuring that parent nodes have a lower cost than their children. It swaps elements as necessary
 * and updates their positions within the heap.
 *
 * @param heap Pointer to the min heap structure.
 * @param index The starting index for the heapify-up operation.
 * 
 * This code used the sources: 
 * https://devashish-iitg.medium.com/heap-sort-heapify-up-or-down-5fd35adfff39
 * and https://www.digitalocean.com/community/tutorials/min-heap-binary-tree
 */
void heapify_up(struct min_heap *heap, int index) {
    while (index != 0 && heap->elements[(index - 1) / 2].cost > heap->elements[index].cost) {
        // Swap the current element with its parent if the current element's cost is less
        struct pq_node temp = heap->elements[(index - 1) / 2];
        heap->elements[(index - 1) / 2] = heap->elements[index];
        heap->elements[index] = temp;

        // Update the positions of the current element and its parent in the positions array
        heap->positions[heap->elements[index].node_id] = index;
        heap->positions[temp.node_id] = (index - 1) / 2;

        // Move up to the parent node for the next iteration
        index = (index - 1) / 2;
    }
}


/**
 * Inserts a new node with a specified cost into the min heap.
 *
 * Adds the node to the heap, ensuring the heap property is maintained by 
 * performing a heapify-up operation.
 * If the heap is full, it prints a message and does not insert the node.
 *
 * @param heap Pointer to the min heap structure.
 * @param node_id The ID of the node to insert.
 * @param cost The cost associated with the node.
 * 
 * This code used the source: 
 * https://www.digitalocean.com/community/tutorials/min-heap-binary-tree
 */
void insert_min_heap(struct min_heap *heap, int node_id, double cost) {
    // Check if the heap has reached its capacity
    if (heap->size >= heap->capacity) {
        printf("Heap is full\n");
        return; // Do not insert if the heap is full
    }

    // Add the new node at the end of the heap
    heap->elements[heap->size] = (struct pq_node){node_id, cost};
    heap->positions[node_id] = heap->size; // Update the position
    heap->size++;

    // Adjust the heap to maintain the heap property
    heapify_up(heap, heap->size - 1);
}


/**
 * Performs the heapify-down operation to maintain the min-heap property.
 *
 * Starting from a specified index, this function ensures that the current node 
 * and its descendants satisfy the min-heap property by swapping nodes as necessary. 
 * It recursively adjusts the positions of nodes downwards until the heap is 
 * correctly ordered.
 *
 * @param heap Pointer to the min heap structure.
 * @param index The starting index for the heapify-down operation.
 * 
 * This code used the sources: 
 * https://devashish-iitg.medium.com/heap-sort-heapify-up-or-down-5fd35adfff39
 * and https://www.digitalocean.com/community/tutorials/min-heap-binary-tree
 */
void heapify_down(struct min_heap *heap, int index) {
    int smallest = index;
    int left = 2 * index + 1; // Calculate left child index
    int right = 2 * index + 2; // Calculate right child index

    // If the left child exists and is smaller than the current node, update smallest
    if (left < heap->size && heap->elements[left].cost < heap->elements[smallest].cost) {
        smallest = left;
    }

    // If the right child exists and is smaller than the smallest so far, update smallest
    if (right < heap->size && heap->elements[right].cost < heap->elements[smallest].cost) {
        smallest = right;
    }

    // If the smallest is not the current index, swap and continue heapifying down
    if (smallest != index) {
        // Swap with the smaller child
        struct pq_node temp = heap->elements[index];
        heap->elements[index] = heap->elements[smallest];
        heap->elements[smallest] = temp;

        // Heapify down the subtree
        heapify_down(heap, smallest);
    }
}


/**
 * Extracts the minimum node from the min heap.
 *
 * Removes and returns the root of the heap, which is the minimum element. 
 * Adjusts the heap to maintain the min-heap property by moving the last element to 
 * the root and performing a heapify-down operation. If the heap is empty, 
 * returns a pq_node with -1 for both node_id and cost as an indicator.
 *
 * @param heap Pointer to the min heap structure.
 * @return The pq_node with the minimum cost in the heap, or an indicator node if the heap is empty.
 * 
 * This code used the source: https://www.geeksforgeeks.org/c-program-to-implement-min-heap/
 */
struct pq_node extract_min(struct min_heap *heap) {
    if (heap->size <= 0) {
        return (struct pq_node){-1, -1}; // Empty heap indicator
    }

    // Store the root node to return
    struct pq_node root = heap->elements[0];

    // Move the last element to the root and decrease heap size
    heap->elements[0] = heap->elements[heap->size - 1];
    heap->size--;

    // Reorganize the heap to maintain the min-heap property
    heapify_down(heap, 0);

    return root;
}


/**
 * Frees the memory allocated for a min heap.
 *
 * Releases the memory used by the heap's elements, positions array, and the heap structure itself.
 * If the heap pointer is NULL, the function does nothing.
 *
 * @param heap Pointer to the min heap structure to be freed.
 */
void free_min_heap(struct min_heap *heap) {
    if (heap) { // Ensure the heap pointer is not NULL
        free(heap->elements);
        free(heap->positions);
        free(heap);
    }
}


/**
 * Checks if a min heap is empty.
 *
 * @param heap Pointer to the min heap structure to check.
 * @return `true` if the heap is empty or NULL, otherwise `false`.
 */
bool is_empty_min_heap(struct min_heap *heap) {
    if (heap == NULL) {
        return true; // If the heap is NULL, it is considered empty
    }
    return heap->size == 0;
}


/**
 * Updates the distance for a neighbor node if a shorter path is found.
 *
 * Calculates the alternative distance to a neighbor node considering the current node's distance
 * and the speed of the way connecting them. If this alternative distance is shorter than the
 * currently recorded distance for the neighbor, updates the neighbor's distance and predecessor,
 * and inserts the neighbor into the min heap with the updated distance.
 *
 * @param heap The min heap used for prioritizing node processing based on distance.
 * @param m Pointer to the ssmap structure representing the graph.
 * @param dist Array of distances from the start node to every other node.
 * @param prev Array of predecessors for each node to reconstruct the path.
 * @param current The current node being processed.
 * @param neighbor The neighbor node to potentially update.
 * @param way_speed The speed limit of the way connecting the current node and the neighbor.
 */
void update_distance_if_shorter(struct min_heap *heap, const struct ssmap* m, 
double dist[], int prev[], int current, int neighbor, double way_speed) {
    double distance = distance_between_nodes(&m->nodes[current], &m->nodes[neighbor]) / (way_speed * 60.0);
    double alt = dist[current] + distance;
    if (alt < dist[neighbor]) {
        dist[neighbor] = alt;
        prev[neighbor] = current;
        insert_min_heap(heap, neighbor, alt);
    }
}


/**
 * Processes all neighbors of a given node in a specific way.
 *
 * Iterates over all nodes in a way to find the given node (u) and then checks its neighbors to
 * potentially update their distances. 
 *
 * @param heap The min heap used for prioritizing node processing based on distance.
 * @param m Pointer to the ssmap structure representing the graph.
 * @param dist Array of distances from the start node to every other node.
 * @param prev Array of predecessors for each node to reconstruct the path.
 * @param u The node whose neighbors are being processed.
 * @param way Pointer to the way structure in which neighbors are being checked.
 */
void process_node_neighbors(struct min_heap *heap, const struct ssmap* m, 
double dist[], int prev[], int u, struct way* way) {
    for (int pos = 0; pos < way->num_nodes; pos++) {
        if (way->node_ids[pos] == u) {
            // If not the first node, or the way is not one-way, consider the previous node
            if (pos > 0 && !way->oneway) {
                int neighbor = way->node_ids[pos - 1];
                update_distance_if_shorter(heap, m, dist, prev, u, neighbor, way->maxspeed);
            }
            // If not the last node, consider the next node
            if (pos < way->num_nodes - 1) {
                int neighbor = way->node_ids[pos + 1];
                update_distance_if_shorter(heap, m, dist, prev, u, neighbor, way->maxspeed);
            }
        }
    }
}


void ssmap_path_create(const struct ssmap* m, int start_id, int end_id) {
    if (start_id == end_id) { // If start and end nodes are the same, print the start node and return
        printf("%d\n", start_id);
        return;
    }
   
    // Check if the provided node IDs are within the valid range
    if (start_id < 0 || start_id >= m->nr_nodes || end_id < 0 || end_id >= m->nr_nodes) {
        printf("error: node %d does not exist.\n", start_id < 0 || start_id >= m->nr_nodes ? start_id : end_id);
        return;
    }

    double* dist = (double*)malloc(m->nr_nodes * sizeof(double)); // Array to hold the distance from start
    int* prev = (int*)malloc(m->nr_nodes * sizeof(int));  // Array to hold the previous node in the optimal path
    struct min_heap* q = init_min_heap(m->nr_nodes); // Priority queue as min heap

    // Initialization
    for (int i = 0; i < m->nr_nodes; i++) {
        dist[i] = INFINITY;
        prev[i] = -1; // UNDEFINED
    }

    dist[start_id] = 0; // Distance from start to itself is 0
    insert_min_heap(q, start_id, dist[start_id]);

    while (!is_empty_min_heap(q)) {
        int u = extract_min(q).node_id;
        for (int i = 0; i < m->nodes[u].num_ways; i++) {
            struct way* way = &m->ways[m->nodes[u].way_ids[i]];
            process_node_neighbors(q, m, dist, prev, u, way);
        }
    }

    // Reconstruct the path from end_id to start_id
    if (prev[end_id] == -1) {
        printf("error: could not find a path from node %d to node %d.\n", start_id, end_id);
    } else {
        // Reconstruct and print the path
        int* path = (int*)malloc(m->nr_nodes * sizeof(int));
        int path_idx = 0, at;
        for (at = end_id; at != -1; at = prev[at]) {
            path[path_idx++] = at;
        }
        // Print the path
        for (int i = path_idx - 1; i >= 0; i--) {
            printf("%d ", path[i]);
        }
        printf("\n");
        free(path); // Free the path array
    }

    // Clean up
    free(dist);
    free(prev);
    free_min_heap(q);
}

