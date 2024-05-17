import xarray as xr
import networkx as nx
import numpy as np
import hydromt
import snakemake
from pathlib import Path

'''
Auth: Joost Buitink
'''

try:
    gridfile = snakemake.input.gridfile

except:
    
    

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

def find_downstream_id(basin_id, sub_map, gauge_map, ldd_map, x_name, y_name):
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


if __name__ == "__main__":
    ds = xr.open_dataset(gridfile)

    layer = "Sall"

    sub = ds[f"wflow_subcatch_{layer}"]
    gauge = ds[f"wflow_gauges_{layer}"]

    # ds = xr.open_dataset(
    #     R"p:\11209265-grade2023\climate_scenarios\scenarios_23\3_Wflow\Rhine\wflow_202307\staticmaps\staticmaps.nc"
    # )
    # layer = "grade_climate"
    # sub = ds[f"wflow_subcatch_{layer}"]
    # gauge = ds[f"wflow_gauges_{layer}"]

    ldd = ds["wflow_ldd"]

    graph, levels_dict = generate_graph_levels(
        subcatchment_map=sub,
        gauge_map=gauge,
        ldd_map=ldd,
    )

    print(levels_dict)

    # Plot results
    import matplotlib.pyplot as plt

    fig = plt.figure("graph", clear=True, tight_layout=True)
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
    ax2.set_title("DAG layout in topological order")
