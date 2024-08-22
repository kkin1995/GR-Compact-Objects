import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import glob
import argparse

M_sun = 1.989e33 # g
R_sun = 69.634e9 # cm
DATA_BASE_PATH = "/Users/karankinariwala/Dropbox/KARAN/2 Areas/Education/PhD/gr_compact_objects/research_project/data"
PLOT_BASE_PATH = "/Users/karankinariwala/Dropbox/KARAN/2 Areas/Education/PhD/gr_compact_objects/research_project/plots"

def plot_exists(filename):
    return os.path.isfile(filename)

def load_data(file_path):
    data = pd.read_csv(file_path)

    if 'r' in data.columns and 'm' in data.columns and 'P' in data.columns:
        data = data.astype({'r': float, 'm': float, 'P': float})
    elif 'Central_Density' in data.columns and 'Mass' in data.columns and 'Radius' in data.columns:
        data = data.astype({'Central_Density': float, 'Mass': float, 'Radius': float})
    else:
        print(f"Unexpected column names in {file_path}. Columns found: {data.columns.tolist()}")

    data = data.dropna()

    return data

def parse_profile_filename(filename):
    pattern = r'(\w+)_profile_central_density_(\d+\.\d+e[+-]?\d+)\.csv'
    match = re.search(pattern, filename)
    if match:
        gas_type = match.group(1)
        central_density = float(match.group(2))
        return gas_type, central_density
    return None, None

def parse_summary_filename(filename):
    pattern = r'(\w+)_summary_(\d+\.\d+e[+-]?\d+)_to_(\d+\.\d+e[+-]?\d+)\.csv'
    match = re.search(pattern, filename)
    if match:
        gas_type = match.group(1)
        start_density = float(match.group(2))
        end_density = float(match.group(3))
        return gas_type, start_density, end_density
    return None, None, None

def plot_stellar_profile(data, gas_type, central_density):
    filename = f"{PLOT_BASE_PATH}/{gas_type}_profile_{central_density:.2e}.png"
    if plot_exists(filename):
        return  

    fig, ax1 = plt.subplots(figsize = (10, 6))

    ax1.scatter(data['r'] * 1e-5, data['m'] / M_sun, s=10, c='red', label='Mass')
    ax1.set_xlabel('Radius (r) [km]')
    ax1.set_ylabel('M / M_sun', color = 'red')
    ax1.tick_params(axis='y', labelcolor='red')

    ax2 = ax1.twinx()
    ax2.scatter(data['r'] * 1e-5, data['P'], s=10, c='blue', label='Pressure')
    ax2.set_ylabel('Pressure (P) [dynn/cm^2]', color='blue')
    ax2.tick_params(axis='y', labelcolor='blue')

    plt.title(f'{gas_type.replace("_", " ").title()}\nCentral Density: {central_density:.2e} g/cm^3')
    fig.tight_layout()
    plt.grid(True)
    plt.savefig(filename)
    plt.close()

def plot_mass_radius_relation(data, gas_type, start_density, end_density):
    filename = f"{PLOT_BASE_PATH}/{gas_type}_mass_radius_relation_{start_density:.2e}_to_{end_density:.2e}.png"
    if plot_exists(filename):
        return
    
    plt.figure(figsize=(10,6))
    plt.scatter(data['Radius'] * 1e-5, data['Mass'] / M_sun, s=15, c='red')
    plt.xlabel('Radius (R) [km]')
    plt.ylabel('Mass (M / M_sun)')
    plt.title(f"Mass-Radius Relation for {gas_type.replace("_", " ").title()}\nCentral Density: {start_density:.2e} to {end_density:.2e} g/cm^3")
    plt.grid(True)
    plt.savefig(filename)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Plot Stellar Structure Data")
    parser.add_argument("--replot", action="store_true", help="Force Replotting of All Data")
    args = parser.parse_args()

    global plot_exists
    if args.replot:
        plot_exists = lambda filename: False

    # Plot Individual Stellar Profiles
    csv_files = glob.glob(f"{DATA_BASE_PATH}/*.csv", recursive=True)
    print(f"Found {len(csv_files)} CSV Files")
    profile_files = []
    summary_files = []

    # Categorize files based on their names
    for file in csv_files:
        filename = os.path.basename(file)
        if '_profile_central_density_' in filename:
            profile_files.append(file)
        elif '_summary_' in filename:
            summary_files.append(file)

    print(f"Categorized {len(profile_files)} profile files and {len(summary_files)} summary files")
    
    for idx, filename in enumerate(profile_files, 1):
        print(f"Processing Stellar Profile File: {idx} / {len(profile_files)}")
        data = load_data(filename)
        gas_type, central_density = parse_profile_filename(os.path.basename(filename))
        if gas_type and central_density:
            plot_stellar_profile(data, gas_type, central_density)
        else:
            print(f"Couldn't Parse Filename: {filename}")

    for idx, filename in enumerate(summary_files, 1):
        print(f"Processing Summary File: {idx} / {len(profile_files)}")
        data = load_data(filename)
        gas_type, start_density, end_density = parse_summary_filename(os.path.basename(filename))
        if gas_type and start_density and end_density:
            plot_mass_radius_relation(data, gas_type, start_density, end_density)
        else:
            print(f"Couldn't Parse Summary Filename: {filename}")

if __name__ == "__main__":
    main()