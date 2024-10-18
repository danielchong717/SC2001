#include <vector>
#include <cstdlib>
#include <chrono>
#include <limits.h>
#include <random>
using namespace std;
using namespace std::chrono;

// Function to generate adjacency matrix
vector<vector<int>> generateAdjacencyMatrix(int vertices, int edges) {
    vector<vector<int>> matrix(vertices, vector<int>(vertices, 0));
    int u = 0;
    int edgeCount = 0;
    while (edgeCount < edges) {
        int v = 0;
        while (v < vertices) {
            if (u != v && matrix[u][v] == 0 && matrix[v][u] == 0) {
                int weight = rand() % 100 + 1; // Random weight between 1 and 100
                matrix[u][v] = weight;
                matrix[v][u] = weight;
                edgeCount++;
            }
            v++;
        }
        u++;
    }
    return matrix;
}

#include <vectpr>
#include <limits>
#include <utility>
#include <chrono>
#include <queue>   
#include <cstdlib> 
#include <random>
using namespace std;
using namespace std::chrono;

// Function to generate adjacency list
vector<vector<pair<int, int>>> generateAdjacencyList(int vertices, int edges) {
    vector<vector<pair<int, int>>> adjacencyList(vertices);
    int edgeCount = 0;
    while (edgeCount < edges) {
        int u = rand() % vertices; 
        int v = rand() % vertices; 
        if (u != v) {
            bool exists = false;
            for (const auto& edge : adjacencyList[u]) {
                /**
                 * Here, edge.first is the index of the connected vertex, 
                 * where edge.second refers to the weight of the edge.
                */
                if (edge.first == v) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                int weight = rand() % 100 + 1; // Random weight between 1 and 100
                adjacencyList[u].emplace_back(v, weight);
                adjacencyList[v].emplace_back(u, weight);
                edgeCount++;
            }
        }
    }
    return adjacencyList;
}

// Function to implement Dijkstra's algorithm using an array-based priority queue
pair<int, double> dijkstraArray(const vector<vector<int>> graph, int src) {
    int vertices = graph.size();
    vector<int> dist(vertices, INT_MAX);
    vector<bool> visited(vertices, false);
    dist[src] = 0;
    int comparisons = 0;
    auto start_time = high_resolution_clock::now();

    for (int i = 0; i < vertices - 1; i++) {
        // Find the vertex with the minimum distance value from the set of vertices not yet visited
        int u = -1;
        for (int v = 0; v < vertices; v++) {
            if (!visited[v] && (u == -1 || dist[v] < dist[u])) {
                u = v;
                comparisons++;
            }
        }
        visited[u] = true;
        // Matches back the corresponding vertex selected with the minimum cost from the other vertices.
        for (int v = 0; v < vertices; v++) {
            if (!visited[v] && graph[u][v] && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
                comparisons++;
            }
        }
    }
    auto end_time = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(end_time - start_time).count() / 1000.0;
    return {comparisons, duration};
}
#include <vector>
#include <utility> // for std::pair
#include <queue>   // for std::priority_queue
#include <limits>  // for INT_MAX
#include <chrono>

using namespace std;
using namespace std::chrono;

#include <vector>
#include <utility> // for std::pair
#include <queue>   // for std::priority_queue
#include <limits>  // for INT_MAX
#include <chrono>

using namespace std;
using namespace std::chrono;

// Function to implement Dijkstra's algorithm using a min-heap (priority queue)
pair<int, double> dijkstraMinHeap(const vector<vector<pair<int, int>>>& graph, int src) {
    int vertices = graph.size();
    vector<int> dist(vertices, INT_MAX);
    vector<bool> visited(vertices, false);
    dist[src] = 0;
    int comparisons = 0;
    auto start_time = high_resolution_clock::now();

    // Min-heap (priority queue) to store (distance, vertex)
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> minHeap;
    minHeap.push({0, src}); 
    while (!minHeap.empty()) {
        // Extract the vertex with the minimum distance
        int u = minHeap.top().second;
        minHeap.pop();
        if (visited[u]) continue;
        visited[u] = true;
        // Checks for a more efficient path from its neighbors than the current
        for (const auto& edge : graph[u]) {
            int v = edge.first;
            int weight = edge.second;
            if (!visited[v] && dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                minHeap.push({dist[v], v}); // Add the updated distance to the min-heap
                comparisons++;
            }
        }
    }

    auto end_time = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(end_time - start_time).count() / 1000.0;
    return {comparisons, duration};
}

// Main function to generate multiple sparse graphs and dense graphs as adjacency matrices, all saved to two files.
int main() {
    ofstream adjacencyMatrixData("adjacencyMatrixData.csv");
    ofstream adjacencyListData("adjacencyListData.csv");
    if (!adjacencyMatrixData.is_open() && !adjacencyListData.is_open()) {
        cerr << "Error opening file!" << endl;
        return 1;
    }
    else {
        adjacencyMatrixData << "Vertices,Comparisons(sparse),Duration(sparse),Comparisons(dense),Duration(dense)\n";
        for (int vertices = 2; vertices <= 800; vertices++) {
            //Sparse graph adjacency matrix generation.
            int edge = vertices - 1;
            vector<vector<int>> sparseMatrix = generateAdjacencyMatrix(vertices, edge);
            //Dense graph adjacency matrix generation.
            edge = (vertices - 1)*(vertices) / 2;
            vector<vector<int>> denseMatrix = generateAdjacencyMatrix(vertices, edge);
            auto sparseRow = dijkstraMatrix(sparseMatrix, 0);
            auto denseRow = dijkstraMatrix(denseMatrix, 0);
            adjacencyMatrixData << vertices << "," << sparseRow.first << "," << sparseRow.second << "," << denseRow.first << "," << denseRow.second;
            if (vertices < 800)
                adjacencyMatrixData << "\n";
        }
    }
    adjacencyMatrixData.close();
    cout << "Finished Dijkstra's Algorithm on 2D Array Adjacency Matrix." << endl;
    return 0;
}