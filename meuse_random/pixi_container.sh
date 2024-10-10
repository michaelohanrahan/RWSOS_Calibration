#!/bin/bash

# Define variables
TEMP_DIR=~/temp
SIF_IMAGE="$TEMP_DIR/pixi_image.sif"
DEFINITION_FILE="/mnt/c/git/pixi_container_snake/apptainer_dep_dual_image/pixi_env.def"  # Adjust the path as necessary
EXPORT_PATH="/mnt/c/git/pixi_container_snake/apptainer_dep_dual_image/pixi_image.sif"
PIXILOCK_FILE="/mnt/c/git/pixi_container_snake/apptainer_dep_dual_image/pixi.lock"  # Adjust the path as necessary
PIXITOML_FILE="/mnt/c/git/pixi_container_snake/apptainer_dep_dual_image/pixi.toml"  # Adjust the path as necessary

# Create temporary directory
mkdir -p $TEMP_DIR

# Copy the definition file, lock file, and TOML file to the temporary directory
cp $DEFINITION_FILE $TEMP_DIR
cp $PIXILOCK_FILE $TEMP_DIR
cp $PIXITOML_FILE $TEMP_DIR

# Build the Singularity image with fakeroot
singularity build --fakeroot $SIF_IMAGE "$TEMP_DIR/pixi_env.def"

# Check if the build was successful
if [ $? -eq 0 ]; then
    echo "Singularity image created successfully."
else
    echo "Failed to create Singularity image."
    exit 1
fi

# Move the image to the desired location
cp $SIF_IMAGE $EXPORT_PATH

# Check if the copy was successful
if [ $? -eq 0 ]; then
    echo "Image copied to $EXPORT_PATH successfully."
else
    echo "Failed to copy the image."
    exit 1
fi

# Clean up temporary files
rm -rf $TEMP_DIR  # Remove the temporary directory

echo "Process completed successfully."