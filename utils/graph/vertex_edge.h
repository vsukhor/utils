/* ==============================================================================

 Copyright (C) 2009-2019, Valerii Sukhorukov <vsukhorukov@yahoo.com>

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
 * \file vertex_edge.h
 * \brief Description of vertexes and edges for use in graphs.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_GRAPH_VERTEX_EDGE_H
#define UTILS_GRAPH_VERTEX_EDGE_H

#include "../common/misc.h"

/// Library-wide.
namespace Utils {
/// Abstract graph circuitry.
namespace Graph {

using namespace Utils::Common;

 /// Weight classification of graph edges.
enum class EdgeMode {
	weighted,
	unweighted
};


/**
  * \brief Implements a graph edge connected to a vertex.
  * \tparam fT Floating point type.
  * \tparam iT Integer type.
  * \tparam LinkMode Weitghted/unweighted mode.
 */
template <typename fT, typename iT, auto LinkMode>
struct EdgeType {};


/**
  * \brief Specification of \a EdgeType for weighted edges.
  * \details Implements a graph edge connected to a vertex in weighted mode.
  * \tparam fT Floating point type.
  * \tparam iT Integer type.
 */
template <typename fT, typename iT>
struct EdgeType<fT,iT,EdgeMode::weighted> {

	using vertex_t = iT;	///< Type alias for vertex index.
	using weight_t = fT;	///< Type alias for edge weight.

	static constexpr weight_t max_weight {huge<weight_t>};	///< maximal edge weight allowed

	vertex_t target;	///< Target neighbour index.
	weight_t weight;	///< Weight of the connecting edge.

	/**
	 * \brief Constructor.
	 * \param target Vertex to which this edge is to bind the parent vertex.
	 * \param weight Weight of the edge (default is one).
	 */
	explicit EdgeType(
			vertex_t target,
			weight_t weight=one<weight_t>
			)
		: target {target}
		, weight {weight}
	{}
};


/**
  * \brief Specification of \a EdgeType for unweighted edges.
  * \details Implements a graph edge connected to a vertex in unweighted mode.
  * \tparam fT Floating point type.
  * \tparam iT Integer type.
 */
template <typename fT, typename iT>
struct EdgeType<fT,iT,EdgeMode::unweighted> {

	using vertex_t = iT;	///< Type alias for vertex index.
	using weight_t = fT;	///< Type alias for edge weight.

	static constexpr weight_t max_weight {huge<weight_t>};	///< Maximal edge weight allowed.
	static constexpr weight_t weight {one<weight_t>};		///< Weight of the connecting edge.

	vertex_t target;	///< Target neighbour index.

	/**
	 * \brief Constructor.
	 * \param target Target vertex.
	 */
	explicit EdgeType(vertex_t target)
		: target {target}
	{}
};


}	// namespace Graph
}	// namespace Utils

#endif  // UTILS_GRAPH_VERTEX_EDGE_H