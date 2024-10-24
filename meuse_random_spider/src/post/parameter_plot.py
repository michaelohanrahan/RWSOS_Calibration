from pathlib import Path
import xarray as xr
import logging
import os
import json
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import geopandas as gpd
import matplotlib.colors as colors
import seaborn as sns

# Set up basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MapLoader:
    def __init__(self, work_dir, filename, name):
        self.work_dir = work_dir
        self.filename = filename
        self.name = name
        self.data = None

    def load_map(self, map_layers=None):
        if self.data is None:
            file_path = Path(self.work_dir, self.filename)
            try:
                if not file_path.exists():
                    raise FileNotFoundError(f"File not found: {file_path}")
                
                if map_layers:
                    self.data = xr.open_dataset(file_path)[map_layers]
                    logging.info(f"Loaded map: {self.name} from {file_path} (layers: {', '.join(map_layers)})")
                else:
                    self.data = xr.open_dataset(file_path)
                    logging.info(f"Loaded all layers of map: {self.name} from {file_path}")
            except FileNotFoundError as e:
                logging.error(f"Error loading {self.name}: {str(e)}")
            except KeyError as e:
                logging.error(f"KeyError when loading {self.name} from {file_path}: {str(e)}")
                logging.info(f"File exists: {file_path.exists()}, Size: {file_path.stat().st_size if file_path.exists() else 'N/A'}")
            except Exception as e:
                logging.error(f"Unexpected error loading {self.name} from {file_path}: {str(e)}")
        return self.data
    
    def get_name(self):
        return self.name

class MapCollection:
    def __init__(self):
        self.maps = {}

    def add_map(self, key, map_loader):
        self.maps[key] = map_loader
        logging.info(f"Added map to collection: {key}")

    def load_all_maps(self, map_layers=None):
        for map_loader in self.maps.values():
            map_loader.load_map(map_layers)
        logging.info(f"Loaded all {len(self.maps)} maps" + (f" (layers: {', '.join(map_layers)})" if map_layers else ""))

    def get_map(self, key):
        return self.maps[key]

class RecipeLoader:
    def __init__(self, recipe_path):
        self.recipe_path = Path(recipe_path)
        self.recipe = None

    def load_recipe(self):
        try:
            with self.recipe_path.open('r') as f:
                self.recipe = json.load(f)
            logging.info(f"Loaded recipe from {self.recipe_path}")
            self._process_recipe()
        except Exception as e:
            logging.error(f"Error loading recipe from {self.recipe_path}: {str(e)}")

    def _process_recipe(self):
        for key, value in self.recipe.items():
            if isinstance(value, str) and ',' in value:
                self.recipe[key] = [item.strip() for item in value.split(',')]
        logging.info("Processed recipe, splitting comma-separated values")

    def get_recipe(self):
        return self.recipe

class MapPlotter:
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load geojson files
        self.subcatch = gpd.read_file(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random_spider\data\3-input\staticgeoms\subcatch_Hall.geojson')

        # Create a custom colormap
        self.diff_cmap = self.create_custom_cmap()

        # Add a new non-diverging colormap for min-max difference
        self.minmax_cmap = plt.cm.viridis

    def create_custom_cmap(self):
        # Create a custom colormap: red (negative) to white (zero) to green (positive)
        colors_list = ['#FF4136', '#FFFFFF', '#2ECC40']
        n_bins = 100  # Number of discrete color levels
        return colors.LinearSegmentedColormap.from_list('custom_diverging', colors_list, N=n_bins)

    def _create_map(self, data, map_name, layer_name, cmap='viridis', vmin=None, vmax=None):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        
        if vmax is None:
            vmax = np.nanpercentile(data.values, 80)  # Set vmax to 80th percentile
        if vmin is None:
            vmin = np.nanmin(data.values)
        
        im = ax.pcolormesh(data.longitude, data.latitude, data, cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, label=layer_name, pad=0.08)
        
        # Add internal borders (subcatchments)
        self.subcatch.boundary.plot(ax=ax, color='white', linewidth=0.5, linestyle=':', alpha=0.2, zorder=2, transform=ccrs.PlateCarree())
        
        # Add country borders
        ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.5)
        
        # Set extent to match the data
        ax.set_extent([data.longitude.min(), data.longitude.max(), data.latitude.min(), data.latitude.max()], crs=ccrs.PlateCarree())
        
        # Add gridlines and labels
        gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        
        # Set title and adjust layout
        ax.set_title(f"{map_name} - {layer_name}")
        
        # Adjust layout to remove whitespace
        plt.tight_layout()
        
        return fig

    def plot_single_maps(self, map_collection, map_layers, output_dir):
        method = 'single'
        for map_name, map_loader in map_collection.maps.items():
            map_data = map_loader.data
            for layer in map_layers:
                if layer in map_data:
                    data = map_data[layer]
                    fig = self._create_map(data, map_name, layer)
                    
                    # Save the plot
                    output_now = Path(output_dir, map_name, method)
                    output_now.mkdir(parents=True, exist_ok=True)
                    fig.savefig(Path(output_now, f"{map_name}_{layer}.png"), dpi=400, bbox_inches='tight')
                    plt.close(fig)
                    
                    logging.info(f"Saved plot for: {map_name}_{layer}")
                else:
                    logging.warning(f"Layer '{layer}' not found in parameter map '{map_name}'")

    def plot_diff_maps(self, map_collection, map_layers, output_dir):
        method = 'diff'
        base_model = map_collection.get_map('base_model').data
        for map_name, map_loader in map_collection.maps.items():
            if map_name == 'base_model':
                continue
            map_data = map_loader.data
            for layer in map_layers:
                if layer in map_data and layer in base_model:
                    diff_data = map_data[layer] - base_model[layer]
                    vmax = np.max([np.abs(np.nanpercentile(diff_data.values, 2)), np.abs(np.nanpercentile(diff_data.values, 98))])
                    vmin = -vmax
                    
                    fig = self._create_map(diff_data, f"{map_name}_diff", layer, cmap=self.diff_cmap, vmin=vmin, vmax=vmax)
                                
                    # Save the plot
                    output_now = Path(output_dir, map_name, method)
                    output_now.mkdir(parents=True, exist_ok=True)
                    fig.savefig(Path(output_now, f"{map_name}_diff_{layer}.png"), dpi=400, bbox_inches='tight')
                    plt.close(fig)
                    
                    logging.info(f"Saved diff plot for: {map_name}_diff_{layer}")
                else:
                    logging.warning(f"Layer '{layer}' not found in parameter map '{map_name}' or base model")

    def plot_minmax_diff_maps(self, map_collection, map_layers, output_dir):
        method = 'minmax_diff'
        top_n_maps = [map_loader.data for name, map_loader in map_collection.maps.items() if name.startswith('Top_')]
        
        for layer in map_layers:
            if all(layer in map_data for map_data in top_n_maps):
                layer_data = [map_data[layer] for map_data in top_n_maps]
                min_data = xr.concat(layer_data, dim='ensemble').min(dim='ensemble')
                max_data = xr.concat(layer_data, dim='ensemble').max(dim='ensemble')
                diff_data = max_data - min_data
                
                vmax = np.nanpercentile(diff_data.values, 98)
                vmin = np.nanpercentile(diff_data.values, 2)
                
                fig = self._create_map(diff_data, "Min-Max Difference", layer, cmap=self.minmax_cmap, vmin=vmin, vmax=vmax)
                
                # Save the plot
                output_now = Path(output_dir, method)
                output_now.mkdir(parents=True, exist_ok=True)
                fig.savefig(Path(output_now, f"minmax_diff_{layer}.png"), dpi=400, bbox_inches='tight')
                plt.close(fig)
                
                logging.info(f"Saved min-max diff plot for: minmax_diff_{layer}")
            else:
                logging.warning(f"Layer '{layer}' not found in all Top_n parameter maps")

    def plot_all_maps(self, map_collection, map_layers, output_dir):
        # self.plot_single_maps(map_collection, map_layers, output_dir)
        # self.plot_diff_maps(map_collection, map_layers, output_dir)
        # self.plot_minmax_diff_maps(map_collection, map_layers, output_dir)
        # self.plot_spatial_correlation_matrix(map_collection, map_layers, output_dir)

    def plot_spatial_correlation_matrix(self, map_collection, map_layers, output_dir):
        base_model = map_collection.get_map('base_model').data
        corr_matrix = np.zeros((len(map_layers), len(map_layers)))
        
        for i, layer1 in enumerate(map_layers):
            for j, layer2 in enumerate(map_layers):
                if layer1 in base_model and layer2 in base_model:
                    corr = np.corrcoef(base_model[layer1].values.flatten(), 
                                       base_model[layer2].values.flatten())[0, 1]
                    corr_matrix[i, j] = corr
        
        fig, ax = plt.subplots(figsize=(12, 10))
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, 
                    xticklabels=map_layers, yticklabels=map_layers, ax=ax)
        ax.set_title('Spatial Correlation Between Parameters')
        
        output_path = Path(output_dir, 'spatial_correlation_matrix.png')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved spatial correlation matrix to {output_path}")

if __name__ == "__main__":
    work_dir = Path(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random_spider')
    map_collection = MapCollection()

    keys = [f"Top_{n}" for n in range(1, 10+1)]
    keys.append('base_model')
    # Add Top_n maps
    for key in keys:
        map_loader = MapLoader(work_dir, f'data/3-input/input_{key}/staticmaps.nc', key)
        map_collection.add_map(key, map_loader)
    
    recipe_path = work_dir / 'config' / 'LHS_calib_recipe.json'
    recipe_loader = RecipeLoader(recipe_path)
    recipe_loader.load_recipe()
    recipe = recipe_loader.get_recipe()
    map_layers = [key for key in recipe.keys()]
    map_layers = [item.strip() for sublist in [key.split(',') for key in map_layers] for item in sublist]
    map_layers = map_layers[:]
    
    # Load all maps
    map_collection.load_all_maps(map_layers)

    # Create MapPlotter instance
    output_dir = work_dir / 'data' / '5-visualization' / 'param'
    map_plotter = MapPlotter(output_dir)

    # Plot all maps
    map_plotter.plot_all_maps(map_collection, map_layers, output_dir)
