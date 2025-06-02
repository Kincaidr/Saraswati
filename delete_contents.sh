#!/bin/bash

# Check if cluster_name is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <cluster_name>"
  exit 1
fi

# Assign the first argument to cluster_name
cluster_name=$1

# Define directories based on cluster_name
sim_catalogs_dir="${cluster_name}/simulation/sim_catalogs"
sim_images_dir="${cluster_name}/simulation/sim_images"
sim_plots_dir="${cluster_name}/simulation/sim_plots"

# Function to delete contents of a directory without removing the directory itself
delete_contents() {
  local dir=$1
  if [ -d "$dir" ]; then
    # Delete all contents within the directory
    rm -rf "${dir:?}/"*
    echo "Deleted contents of $dir"
  else
    echo "Directory $dir does not exist, skipping."
  fi
}

# Delete contents in each target directory
delete_contents "$sim_catalogs_dir"
delete_contents "$sim_images_dir"
delete_contents "$sim_plots_dir"

