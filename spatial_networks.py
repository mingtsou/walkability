import networkx as nx
import osmnx as ox
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from algorithms import redundant_paths, find_optimal_score

from shapely.geometry import LineString
from shapely.ops import linemerge

METERS_PER_MILE = 1609.34
INF = float('INF')

"""
UTILITY FUNCTIONS
"""
def _set_if_missing(d, **kwargs):
    kwargs.update(d)

    return kwargs

def _ensure_geometry(G, u, v):
    """This function checks whether there is a "geometry" field, and if it 
    doesn't exist it will create "geometry" that is simply a straight line
    between the two nodes.
    """

    edge = G[u][v][0] # get edge with min length (need to change!!!)

    if 'geometry' in edge:
        return edge['geometry']

    u_node = G.node[u]
    v_node = G.node[v]

    x1, y1 = u_node['x'], u_node['y']
    x2, y2 = v_node['x'], v_node['y']

    return LineString(((x1, y1), (x2, y2)))


class Path(object):

    def __init__(self, G, node_order):
        self.G = G
        self.subG = G.subgraph(node_order)
        self.nodes = node_order
        self.edge_pairs = zip(self.nodes[:-1], self.nodes[1:])

        geoms = [_ensure_geometry(self.subG, u, v) for u, v in self.edge_pairs]
        self.geometry = linemerge(geoms)

    def get_bbox(self):
        xs, yx = self.geometry.xy

        return (max(xs), min(xs), max(ys), min(ys))

    def to_gdf(self, label=None):
        rows = []

        for u, v in self.edge_pairs:
            row = self.G[u][v][0]
            row['name'] = row.get('name', '')
            if type(row['name']) is list:
                row['name'] = '/'.join(row['name'])
            row['from'] = u
            row['to'] = v

            if label is not None:
                row['label'] = label

            if 'geometry' not in row:
                x1, y1 = self.G.node[u]['x'], self.G.node[u]['y']
                x2, y2 = self.G.node[v]['x'], self.G.node[v]['y']

                row['geometry'] = LineString(((x1, y1), (x2, y2)))

            rows.append(row)
            
        return gpd.GeoDataFrame(rows, crs={'init': 'epsg:4326', 'no_defs': True})

    def plot(self, ax=None, **kwargs):
        if not ax:
            close = kwargs.pop('close', False)
            save = kwargs.pop('save', False)
            show = kwargs.pop('show', False)
            node_size = kwargs.pop('node_size', 0)

            fig, ax = ox.plot(self.G, node_size=node_size, close=close,
                            save=save, show=show, **kwargs)


        xs, ys = self.geometry.xy
        coords = zip(xs, ys)
        lines = zip(coords[:-1], coords[1:])
        endpoints_xs = (xs[0], xs[-1])
        endpoints_ys = (ys[0], ys[-1])

        path_nodesize = kwargs.pop('path_nodesize', 30)
        path_nodealpha = kwargs.pop('path_nodealpha', 1.)
        path_linealpha = kwargs.pop('path_linealpha', 1.)
        path_linewidth = kwargs.pop('path_linewidth', 4)
        
        ax.scatter(endpoints_xs, endpoints_ys, s=path_nodesize,
                    alpha=path_nodealpha, edgecolor='#444444',
                    color=['#444444', '#ffffff'], zorder=5)

        ax.plot(xs, ys, '#ff2222', alpha=path_linealpha, 
                linewidth=path_linewidth, solid_capstyle='round', zorder=3)

        return ax

    @classmethod 
    def find_optimal_path(cls, paths, path_label='label', criteria=['length'],
                            weights=None, minimize=True):
        
        gdf = Path.paths_to_gdf(paths)
        scores = gdf.groupby(path_label)[criteria].sum().values

        best = find_optimal_score(scores, path_label=path_label, 
                         criteria=criteria, weights=weights, minimize=minimize)

        return paths[best]

    @classmethod
    def paths_to_gdf(cls, paths):
        gdfs = [p.to_gdf(i) for i, p in enumerate(paths)]

        return gpd.GeoDataFrame(pd.concat(gdfs))


class SpatialNetwork(object):

    def __init__(self, G):
        self.net = G

    def shortest_path(self, source, target):
        nodes = nx.shortest_path(self.net, source, target, 'length')

        return Path(self.net, nodes)

    def redundant_paths(self, source, target, coeff, radius=METERS_PER_MILE):
        paths = redundant_paths(self.net, source, target, 
                                'length', coeff, radius)

        return [Path(self.net, nodes) for nodes, edges in paths]
    
    def add_edge_criteria(self, criteria_name, criteria):
        for point in criteria:
            '''
            # Not implemented!
            edge_id = #get_nearest_edge(self.net, point)
            edge = self.net[u][v]
            '''
            node = self.net.node[n_id]
            node['pois'] = node.get('pois', [])
            node['pois'].append(poi)

    def add_pois(self, pois):
        for poi in pois:
            point = poi['coordinates']
            n_id = ox.get_nearest_node(self.net, point)
            node = self.net.node[n_id]
            node['pois'] = node.get('pois', [])
            node['pois'].append(poi)

    def to_catchment(self, center, radius):
        return Catchment(self.net, center, radius)

    def plot(self, **kwargs):
        return ox.plot_graph(self.net, **kwargs)

    def plot_paths(self, paths, **kwargs):
        graph_kwargs = dict((k, v) for k, v in kwargs.items()
                            if not k.startswith('path'))
        graph_kwargs = _set_if_missing(graph_kwargs, show=False, 
                                        save=False, close=False)

        fig, ax = self.plot(**graph_kwargs)

        path_kwargs = dict((k, v) for k, v in kwargs.items() 
                            if k.startswith('path'))

        for p in paths:
            p.plot(ax=ax, **path_kwargs)

        return fig, ax

    @classmethod
    def from_osm_bbox(cls, bbox):
        """Grab OpenStreetMap road data within bounding box.
        NOTE: Put this in a separate "data" module?
        """
        G = ox.graph_from_bbox(*bbox, network_type='walk')

        return cls(G)


# inherit from SpatialNetwork?
class Catchment(object):
    """A catchment is essentially an ego network (a network centered around
    a specific node) limited by a specified network radius.
    """

    def __init__(self, G, center, radius):
        self.radius = radius
        self.center = center
        self.net = nx.ego_graph(G, center, radius=radius, distance='length')

    def shortest_path(self, target):
        nodes = nx.shortest_path(self.net, self.center, target, 'length')

        return Path(self.net, nodes)

    def redundant_paths(self, target, coeff, radius=METERS_PER_MILE):
        paths = redundant_paths(self.net, self.center, target, 
                                'length', coeff, radius)

        return [Path(self.net, nodes) for nodes, edges in paths]

    def plot(self, **kwargs):
        return ox.plot_graph(self.net, **kwargs)

    def plot_paths(self, paths, **kwargs):
        graph_kwargs = dict((k, v) for k, v in kwargs.items()
                            if not k.startswith('path'))
        graph_kwargs = _set_if_missing(graph_kwargs, show=False,
                                        save=False, close=False)

        fig, ax = self.plot(**graph_kwargs)

        path_kwargs = dict((k, v) for k, v in kwargs.items()
                            if k.startswith('path'))

        for p in paths:
            p.plot(ax=ax, **path_kwargs)

        return fig, ax






        












