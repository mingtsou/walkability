import osmnx as ox
import networkx as nx
import numpy as np
from heapq import heappush, heappop
from itertools import count

# geometry
import geopandas as gpd
from shapely.geometry import *
from shapely.ops import *

INF = float('inf')

"""
UTILITY FUNCTIONS
"""
def _create_weight_func(G, weight):
     def weight_func(u, v):
         edge = G[u][v]

         if G.is_multigraph():
             return min(att.get(weight, 1) for att in edge.values())

         return edge.get(weight, 1)

     return weight_func


def dijkstra(G, source, target, weight, blacklist=set(), max_dist=INF):
    """Bastardized version of nx.shortest_path (single source and target).
    Added ability to avoid blacklisted nodes.
    """

    G_succ = G.succ if G.is_directed() else G.adj

    push = heappush
    pop = heappop

    # dictionary of optimal path to solved nodes
    paths = { source: [ source ] }
    # dictionary of final distances for each solved node
    dist = {}
    # track which nodes have been seen and the shortest path distance
    seen = {}

    c = count()
    fringe = []

    seen[source] = 0.0
    push(fringe, (0.0, next(c), source))

    weight = _create_weight_func(G, 'length')

    while fringe:
        # this is a priority queue, so will always pop minimum distance
        d, _, v = pop(fringe)
        
        # ignore if node has been visited already or is blacklisted
        if v in dist or v in blacklist:
            continue

        dist[v] = d
        #
        # path found! break
        if v == target:
            break

        # loop through neighbors of v
        for u, e in G_succ[v].items():
            vu_cost = dist[v] + weight(u, v)

            if vu_cost <= max_dist and (u not in seen or vu_cost < seen[u]):
                seen[u] = vu_cost
                paths[u] = paths[v] + [u]
                push(fringe, (vu_cost, next(c), u))

    if target in paths and target in dist:
        return paths[target], dist[target]
    
    return None
        
           
def redundant_paths(G, source, target, weight, coeff, max_dist):
    """Code inspired by implementation in Urban Network Analysis 
    toolkit (cityform.mit.edu/projects/urban-network-analysis.html)
    """

    search_result = dijkstra(G, source, target, weight, max_dist=max_dist)

    if not search_result:
        return []

    path, dist = search_result

    if dist > max_dist:
        return []

    # if max_dist is less than shortest path distance * coeff, use that
    available_dist = min(dist * coeff, max_dist)

    path = [source]
    edges = []
    
    paths = _get_paths(G, path, edges, target, weight, available_dist) 

    return paths

def _get_paths(G, nodes, edges, target, weight, available_dist, cache={}):
    """Recursively get redundant paths."""

    G_succ = G.succ if G.is_directed() else G.adj

    output = []
    visited = set(nodes)
    path_info = []
    source = nodes[-1]
    get_weight = _create_weight_func(G, weight)

    for neighbor, neighbor_edge in G_succ[source].items():

        if neighbor in visited:
            continue

        new_available_dist = available_dist - get_weight(source, neighbor)

        if new_available_dist < 0:
            continue

        # shortest path key for cache
        sp_key = (source, neighbor, target)

        if sp_key not in cache:
            sp = dijkstra(G, neighbor, target, 'length', set([source]), 
                            new_available_dist)

            cache[sp_key] = INF if not sp else sp[1]

        if new_available_dist < cache[sp_key]:
            continue

        path_info.append((neighbor, neighbor_edge, new_available_dist))


    for (neighbor, neighbor_edge, new_available_dist) in path_info:
        new_path = nodes + [neighbor]
        new_edges = edges + [neighbor_edge]

        if neighbor == target:
            output.append((new_path, new_edges))
        else:
            output.extend(_get_paths(G, new_path, new_edges, target, 
                                    weight, new_available_dist, cache))

    return output


def find_optimal_path(paths, path_label='label', criteria=['length'],
                        weights=None, minimize=True):

    gdf = Path.paths_to_gdf(paths)

    if type(minimize) in (list, np.array):
        minimize = np.array(minimize)
    else:
        minimize = np.ones(len(criteria)).astype(bool)

    if type(criteria) is not list:
        criteria = [criteria]

    # if no weights are specified, use equal weighting
    if not weights:
        weights = np.ones(len(criteria)) / (1. * len(criteria))

    scores = gdf.groupby(path_label)[criteria].sum().values
    # min/max normalization
    scores = (scores - scores.min()) / (scores.max() - scores.min())
    # reverse if necessary
    scores[:, ~minimize] = 1. - scores[:, ~minimize]
    #apply weight to scores                    
    scores = (scores * weights).sum(axis=1)

    func = 'argmin'

    return paths[getattr(scores, func)()]

