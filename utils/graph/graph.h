/* =============================================================================

 Copyright (C) 2009-2025 Valerii Sukhorukov <vsukhorukov@yahoo.com>

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

============================================================================== */

/**
 * \file graph.h
 * \brief Abstract graphs.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_GRAPH_GRAPH_H
#define UTILS_GRAPH_GRAPH_H

#include <deque>
#include <list>
#include <ostream>
#include <set>
#include <utility>      // std::pair
#include <vector>

#include "../common/misc.h"
#include "vertex_edge.h"

/// Abstract graph circuitry.
namespace utils::graph {

/**
 * \brief Stores routines essential for an abstract graph.
 * \details Graphs are collections of nodes (vertexes) interconnected by edges.
 * \tparam ET Edge type.
 */
template<typename ET>
class Graph {

    using vec2int = vec2<int>;

public:

    // typedefs
    using EdgeT = ET;               ///< Type of the graph edges.
    using vertex_t = ET::vertex_t;  ///< Type of the graph vertexes.
    using weight_t = ET::weight_t;  ///< Weight type of the edges.

    /// Type alias for path over consecutively connected vertexes.
    using pathT = std::vector<vertex_t>;

    /// Type alias for the graph adjacency list.
    using adjLT = vec2<ET>;

    /**
     * \brief Breadth first search on the graph
     * \details Determines if vertex \p tar belongs to the graph.
     * \param[in] ajl Adjacency list of the graph.
     * \param[inout] q Auxiliary deque.
     * \param[inout] visited Auxiliary vector of flags deniting visited status
     * of graph vertexes.
     * \param[in] tar Searched vertex.
     * \return 1/0 if \p tar is found/not found respectively.
     */
    auto bfs(
        const adjLT& ajl,
        std::deque<vertex_t>& q,        // by reference
        std::vector<bool>& visited,     // by reference
        const vertex_t tar
    ) const;

    /**
     * \brief Computes paths connecting a vertex in the graph
     * \details Computes paths starting at vertex \p source to other vertexes
     * in the connected component of the graph specified
     * by the djacency list \p ajl .
     * Implements Dijkstra's algorithm
     * \ref https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm.
     * \param[in] source Vertex from which the paths are computed.
     * \param[in] ajl Adjacency list of the graph.
     * \param[inout] min_distance Minimal distances.
     * \param[inout] previous Path vertexes.
     */
    void compute_paths(
        vertex_t source,
        const adjLT& ajl,
        std::vector<weight_t>& min_distance,   // by reference
        std::vector<vertex_t>& previous        // by reference
    ) const;

    /**
     * \brief The shortest path between two graph vertexes.
     * \details Calculates the shortest path between graph vertexes \p v1
     * and \p v2 using \ref compute_paths(vertex_t, const adjLT&) .
     * \param[in] v1 First vertex.
     * \param[in] v2 Second vertex.
     * \param[in] ajl Adjacency list of the graph.
     * \return The shortest path between vertexes \p v1 and \p v2 .
     */
    auto shortest_path(
        vertex_t v1,
        vertex_t v2,
        const adjLT& ajl
    ) -> pathT;

    /**
     * \brief Resets adjacency list of the graph.
     * \details Resets components of adjacency list \p al
     * to weight values \p w .
     * \param[out] al Adjacency list of the graph.
     * \param[in]  w Desired edge weights.
     */
    void reset_al(
        adjLT& al,        // by reference
        weight_t w
    );

    /**
     * \brief Converts list to a std::vector.
     * \param[in] l Input list.
     */
    static auto list2vector(
        const std::list<vertex_t>& l
    ) -> std::vector<vertex_t>;

    /**
     * \brief Creates graph adjacency matrix from its adjacency list.
     * \param[in] ajl Adjacency list of the graph.
     * \param[out] agm Adjacency matrix of the graph.
     */
    void adjacency_matrix(
        const adjLT& ajl,
        vec2int& agm        // by reference
    ) const;

    /**
     * \brief Creates graph laplacian matrix from its adjacency matrix.
     * \param[in] agm Adjacency matrix of the graph.
     * \param[out] lm Laplacian matrix of the graph.
     */
    void laplacian_matrix(
        const vec2int& agm,
        vec2int& lm
    );

    /**
     * \brief Prints adjacency list of the graph.
     * \details Prints full adjacency list of the graph \p ajl
     * to output stream \p os .
     * \note Please use cautiosly with large graphs.
     * \param ajl Adjacency list of the graph.
     * \param os Output stream.
     */
    static void print_adjacency_list(
        const adjLT& ajl,
        std::ostream& os    // by reference
    );

    /**
     * \brief Prints adjacency list of the graph.
     * \details Prints neighbours of the graph vertex \p v
     * to output stream \p os .
     * \param[in] v Vertex, which neighbours are printed.
     * \param[in] ajl Adjacency list of the graph.
     * \param[in] os Output stream.
     */
    static void print_adjacency_list_line(
        vertex_t v,
        const adjLT& ajl,
        std::ostream& os    // by reference
    );
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename ET>
auto Graph<ET>::
bfs(
    const adjLT& ajl,
    std::deque<vertex_t>& q,     // by reference
    std::vector<bool>& visited,  // by reference
    const vertex_t tar
) const
{
    if (q.empty()) return 0;    // reached end of the search wt finding tar

    auto s {std::move(q.front())};
    q.pop_front();
    if (s == tar) return 1;

    for (const auto& o : ajl[s]) {
        const auto trg {o.target};
        if (!visited[trg]) {
            visited[trg] = true;
            q.push_back(trg);
        }
    }

    return bfs(ajl, q, visited, tar);
}


template<typename ET>
void Graph<ET>::
compute_paths(
    const vertex_t source,
    const adjLT& ajl,
    std::vector<weight_t>& min_distance,  // by reference
    std::vector<vertex_t>& previous       // by reference
    ) const
{
    const auto n = ajl.size();

    min_distance.resize(n);
    std::fill(min_distance.begin(), min_distance.end(), ET::max_weight);
    min_distance[source] = 0;

    previous.resize(n);
    std::fill(previous.begin(), previous.end(), huge<vertex_t>);

    std::set<std::pair<weight_t, vertex_t>> vertex_queue;
    vertex_queue.insert(std::make_pair(min_distance[source], source));

    while (!vertex_queue.empty())    {
        const auto dist = vertex_queue.begin()->first;
        const auto u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());

        // Visit each edge exiting u
        const auto& neighbours = ajl[u];
        for (auto iter = neighbours.begin();
                  iter != neighbours.end();
                  iter++) {
            const auto v = iter->target;
            const auto w = iter->weight;
            const auto distance_through_u = dist + w;

            if (distance_through_u < min_distance[v]) {
                vertex_queue.erase(std::make_pair(min_distance[v], v));

                min_distance[v] = distance_through_u;
                previous[v] = u;
                vertex_queue.insert(std::make_pair(min_distance[v], v));
            }
        }
    }
}


template<typename ET>
auto Graph<ET>::
shortest_path(
    const vertex_t v1,
    const vertex_t v2,
    const adjLT& ajl
) -> pathT
{
    pathT previous;
    std::vector<weight_t> min_distance;
    compute_paths(v1, ajl, min_distance, previous);

    // Shortest path to a specific vertex u:
    auto get_shortest_path_to = [&](vertex_t u) {
        std::list<vertex_t> path;
        for (; u != huge<vertex_t>; u = previous[u])
            path.push_front(u);

        return path;
    };

    if (min_distance[v2] < ET::max_weight) {
        // The shortest path edge sequence from ind1 to v2:
        const auto path_l = get_shortest_path_to(v2);
        return list2vector(path_l);
    }
    return pathT();
}


template<typename ET>
void Graph<ET>::
reset_al(
    adjLT &al,
    const weight_t w
)
{
    for (auto& o : al)
        for (auto& oo : o)
            oo = {huge<vertex_t>, w};
}


template<typename ET>
auto Graph<ET>::
list2vector(
    const std::list<vertex_t>& l
) -> std::vector<vertex_t>
{
    std::vector<vertex_t> v(l.size());
    typename std::list<vertex_t>::const_iterator itl;
    typename std::vector<vertex_t>::iterator itv;
    for (itl=l.begin(), itv=v.begin(); itl!=l.end(); itl++, itv++)
        *itv = *itl;

    return v;
}


template<typename ET>
void Graph<ET>::
adjacency_matrix(
    const adjLT& ajl,
    vec2int &agm
) const
{
    agm.resize(ajl.size());
    for (auto& o : agm)
        o.resize(ajl.size(), 0);

    for (szt i=0; i<ajl.size(); i++)
        for (const auto& o : ajl[i])
            agm[i][o.target] = 1;
}


template<typename ET>
void Graph<ET>::
laplacian_matrix(
    const vec2int& agm,
    vec2int& lm
)
{
    lm = common::Vec2::array_like<int,int>(agm);

    const auto s = lm.size();
    for (szt j=0; j<s; j++) {
        int d {};
        for (szt i=0; i<s; i++) d += agm[j][i];       // j-th node degree
        for (szt i=0; i<s; i++) lm[j][i] = -agm[j][i];
        lm[j][j] = d;
    }
}


template<typename ET>
void Graph<ET>::
print_adjacency_list(
    const adjLT& ajl,
    std::ostream& os
)
{
    for (szt i1=0; i1<ajl.size(); i1++)
        print_adjacency_list_line(i1, ajl, os);
}


template<typename ET>
void Graph<ET>::
print_adjacency_list_line(
    const vertex_t v,
    const adjLT& ajl,
    std::ostream& os
)
{
    os << v;
    for (const auto& o : ajl[v])
        os << " [ " << o.target << " " << o.weight << " ] ";
    os << std::endl;
}

}  // namespace utils::graph

#endif  // UTILS_GRAPH_GRAPH_H
