import xarray as xr
import networkx as nx
import numpy as np
import hydromt
from pathlib import Path
import logging 
import matplotlib.pyplot as plt
import json
import itertools
import os 
import argparse as AP

'''
Auth: Joost Buitink
Refactor: Mike O'Hanrahan

This has been refactored to be more modular, including logging and snakemake workflow compatibility.
'''


def setup_logger(name):
    # Create a custom logger
    logger = logging.getLogger(name)
    
    os.makedirs('data/0-log', exist_ok=True)

    # Create handlers
    c_handler = logging.StreamHandler()
    f_handler = logging.FileHandler('data/0-log/file.log', mode='w')
    c_handler.setLevel(logging.WARNING)
    f_handler.setLevel(logging.ERROR)

    # Create formatters and add it to handlers
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

    return logger



    
# Description of how the LDD values translate the indices [x-offset, y-offset]
# assumes the grid has increasing x and y coordinates
ldd_direction = {
    1: (-1, -1),  # 1
    2: (0, -1),  # 2
    3: (1, -1),  # 3
    4: (-1, 0),  # 4
    5: (0, 0),  # 5
    6: (1, 0),  # 6
    7: (-1, 1),  # 7
    8: (0, 1),  # 8
    9: (1, 1),  # 9
}

def find_downstream_id(basin_id, 
                       sub_map, 
                       gauge_map, 
                       ldd_map, 
                       x_name, 
                       y_name):
    # Find the location of the outlet
    outlet_location = gauge_map.where(gauge_map == basin_id).argmax(...)

    # Find the ldd at that location
    ldd_value = int(
        ldd_map.isel(
            {x_name: outlet_location[x_name], y_name: outlet_location[y_name]}
        ).values
    )

    # Only look for a downstream location if the ldd is not a pit
    if ldd_value != 5:
        x_down = outlet_location[x_name] + ldd_direction[ldd_value][0]
        y_down = outlet_location[y_name] + ldd_direction[ldd_value][1]
        # Find the downstream basin id
        downstream_basin = float(sub_map.isel({x_name: x_down, y_name: y_down}).values)
    else:
        # Set to nan if it is a pit
        downstream_basin = np.nan
    return downstream_basin


def generate_graph_levels(
    subcatchment_map: xr.DataArray,
    gauge_map: xr.DataArray,
    ldd_map: xr.DataArray,
):
    """
    Parameters
    ----------
    subcatchment_map:
        Wflow model layer with the subcatchment IDs
    gauge_map:
        Wflow model layer with the gauges, must be in line with the subcatchment_map
    ldd_map:
        Wflow model layer with the LDD values

    Returns
    -------
    DG:
        Directed graph that describes the subcatchments and their relations
    levels_dict:
        Dictionary with the grouped catchment ids
    """
    
    # Get dimension names
    x_dim = subcatchment_map.raster.x_dim
    y_dim = subcatchment_map.raster.y_dim

    # Sort maps by x-y dimensions, to ensure correct mapping of LDD values
    subcatchment_map = subcatchment_map.sortby([x_dim, y_dim])
    ldd_map = ldd_map.sortby([x_dim, y_dim])
    gauge_map = gauge_map.sortby([x_dim, y_dim])

    # Prepare a directed graph
    DG = nx.DiGraph()

    # Build a graph with nodes for all subbasins (without connections)
    for id in np.unique(subcatchment_map):
        if not np.isnan(id):
            DG.add_node(id)

    # Add connections for the graph, based on the subcatchment map
    for id in np.unique(subcatchment_map):
        if not np.isnan(id):
            # Find the downstream basin
            downstream = find_downstream_id(
                basin_id=id,
                sub_map=subcatchment_map,
                gauge_map=gauge_map,
                ldd_map=ldd_map,
                x_name=x_dim,
                y_name=y_dim,
            )

            # Add edge that connects the current basin to the downstream basin
            if not np.isnan(downstream):
                DG.add_edge(id, downstream)

    # Generate lists for each "level" in the connected graph
    levels = [sorted(generation) for generation in nx.topological_generations(DG)]
    # Convert this data to a dictionary
    levels_dict = {}
    for idx, values in enumerate(levels):
        levels_dict[f"level{idx}"] = values

    return DG, levels_dict

def remove_nodes_by_id(graph, node_ids):
    """
    Remove nodes from the graph by their IDs.

    Parameters:
    - graph: The graph object from which nodes will be removed.
    - node_ids: A list of node IDs to remove from the graph.
    """
    # Remove the nodes
    graph.remove_nodes_from(node_ids)

    return graph

if __name__ == "__main__":
    try:
        try:
            gridfile = snakemake.input.gridfile
            gaugeset = snakemake.params.gaugeset
            testmode = False

            # Set up logger
            logger = setup_logger()

        except:
            ap = AP.ArgumentParser()
            ap.add_argument("--gridfile", type=str, required=True)
            ap.add_argument("--gaugeset", type=str, required=True)
            ap.add_argument("--testmode", type=bool, default=True)
            args = ap.parse_args()

            gridfile = args.gridfile
            gaugeset = args.gaugeset
            testmode = args.testmode

            print(f'cwd: {os.getcwd()}')
            logger = setup_logger(os.getcwd(), 'create_dependency_graph.py')
            print(f'logger: {logger}')
            print(f'gridfile: {gridfile}')
            print(f'gaugeset: {gaugeset}')
            print(f'testmode: {testmode}')
        ds = xr.open_dataset(gridfile)

        sub = ds[f"wflow_subcatch_{gaugeset}"]
        gauge = ds[f"wflow_gauges_{gaugeset}"]

        ldd = ds["wflow_ldd"]

        graph, levels_dict = generate_graph_levels(
            subcatchment_map=sub,
            gauge_map=gauge,
            ldd_map=ldd,
        )
        
        nodes_to_remove = []
        
        graph = remove_nodes_by_id(graph, nodes_to_remove)
        
        print(graph)
        

        logger.info(f"Graph for {gaugeset} has {len(graph.nodes)} nodes and {len(graph.edges)} edges")
        logger.info(f"Levels for {gaugeset} are: {levels_dict}")
        

        fig = plt.figure("graph", figsize=(18,8), clear=True, tight_layout=True)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        nx.draw(graph, with_labels=True, ax=ax1)

        for layer, nodes in enumerate(nx.topological_generations(graph)):
            # `multipartite_layout` expects the layer as a node attribute, so add the
            # numeric layer value as a node attribute
            for node in nodes:
                graph.nodes[node]["layer"] = layer

        # Compute the multipartite_layout using the "layer" node attribute
        pos = nx.multipartite_layout(graph, subset_key="layer")

        nx.draw_networkx(graph, pos=pos, ax=ax2)
        ax2.set_title(f"{gaugeset} DAG layout in topological order")
        
        # Create a dictionary where each key is a node ID and each value is a dictionary of node attributes
        nodes = {node: data for node, data in graph.nodes(data=True)}

        # Add dependencies to each node's attributes
        for node, data in nodes.items():
            # Perform a BFS from the node and add all visited nodes to the _deps list
            data['_deps'] = [n for n in nx.bfs_tree(graph, node, reverse=True) if n != node]

        # Create a dictionary to store layers and their dependencies
        levels_graph = {}

        # Group nodes by their layer and collect all dependencies for each layer
        for layer, nodes_in_layer in itertools.groupby(sorted(nodes.items(), key=lambda x: x[1]['layer']), lambda x: x[1]['layer']):
            nodes_in_layer = list(nodes_in_layer)
            # Flatten the list of dependencies for all nodes in the layer
            deps = set(itertools.chain.from_iterable(nodes[node]['_deps'] for node, _ in nodes_in_layer))
            # Store the layer and its dependencies in the levels_graph
            levels_graph[f'level{layer}'] = {'deps': list(deps), 'elements': [node for node, _ in nodes_in_layer]}

        if testmode:
            print(json.dumps(levels_graph, indent=4))
            plt.show()
            
        else:
            with open(f'data/2-interim/{gaugeset}_nodes_graph.json', "w") as f:
                json.dump(nodes, f, indent=4)
                logger.info(f"Nodes for {gaugeset} have been saved to data/2-interim/{gaugeset}_nodes_graph.json")
                
            with open(f'data/2-interim/{gaugeset}_levels_graph.json', "w") as f:
                json.dump(levels_graph, f, indent=4)
                logger.info(f"Levels for {gaugeset} have been saved to data/2-interim/{gaugeset}_levels_graph.json")
                
            plt.savefig(f'data/5-visualization/{gaugeset}_dependency_graph.png')
            plt.close()

    except:
        logger.error("An error occurred", exc_info=True)
        raise
