#This script plots the rarefaction curve of the samples, by the desired taxonomic level

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from collections import Counter

def load_and_filter_data(file_path, taxonomic_level):
    """
    Load Kraken report and filter data by taxonomic level.
    """
    # Kraken2 report format: percentage, clade_reads, direct_reads, level, tax_id, name
    kraken_data = pd.read_csv(file_path, sep='\t', header=None,
                              names=['percentage', 'clade_reads', 'direct_reads', 'level', 'tax_id', 'name'])
    
    # Filter by taxonomic level and get direct reads (not cumulative)
    filtered_data = kraken_data[kraken_data['level'] == taxonomic_level.upper()]
    
    # Use direct_reads instead of clade_reads for proper counts
    read_counts = filtered_data['direct_reads'].values
    
    # Remove zero counts
    read_counts = read_counts[read_counts > 0]
    
    print(f"Found {len(read_counts)} taxa at level {taxonomic_level.upper()}")
    print(f"Total reads: {sum(read_counts)}")
    
    return read_counts

def rarefaction_curve_improved(read_counts, max_depth, step=None, iterations=100):
    """
    Generate a rarefaction curve using multinomial sampling.
    More efficient and accurate than the expansion method.
    """
    if len(read_counts) == 0:
        return [], []
    
    total_reads = sum(read_counts)
    
    # Start from a small depth and use smaller steps for better curve resolution
    min_depth = max(1, total_reads // 1000)  # Start much earlier
    if step is None:
        step = max(1, total_reads // 100)  # Smaller steps for better resolution
    
    # Create sampling depths starting from min_depth
    depths = list(range(min_depth, min(max_depth, total_reads) + 1, step))
    
    # Add some early points for better curve shape
    early_points = [1, 10, 50, 100, 500, 1000]
    for point in early_points:
        if point < min_depth and point <= total_reads:
            depths.insert(0, point)
    
    depths = sorted(list(set(depths)))  # Remove duplicates and sort
    
    if len(depths) == 0:
        return [], []
    
    # Convert to probabilities
    probabilities = read_counts / total_reads
    
    species_counts = []
    
    print(f"Sampling at {len(depths)} depth points from {min(depths)} to {max(depths)}")
    
    for depth in depths:
        observed_species = []
        
        for _ in range(iterations):
            # Multinomial sampling
            sampled_counts = np.random.multinomial(depth, probabilities)
            # Count how many species have at least one read
            observed = np.sum(sampled_counts > 0)
            observed_species.append(observed)
        
        species_counts.append(np.mean(observed_species))
    
    return depths, species_counts

def extract_sample_name(file_path):
    """
    Extract the sample name from the file path.
    """
    return file_path.split('/')[-1].split('_kraken2_report')[0]

def main(input_file, taxonomic_level, output_file):
    plt.figure(figsize=(12, 8))

    sample_name = extract_sample_name(input_file)
    read_counts = load_and_filter_data(input_file, taxonomic_level)

    if len(read_counts) == 0:
        print(f"⚠️ No reads at level '{taxonomic_level.upper()}' in {input_file}")
        # Create empty plot with warning
        plt.text(0.5, 0.5, f'No data for taxonomic level {taxonomic_level.upper()}', 
                ha='center', va='center', transform=plt.gca().transAxes, fontsize=14)
        plt.title(f'Rarefaction Curve - {sample_name}')
        plt.xlabel('Sequencing Depth')
        plt.ylabel('Unique Taxa')
        plt.savefig(output_file, format='jpg', dpi=300, bbox_inches='tight')
        return

    # Use 80% of total reads as max depth to see the curve better
    total_reads = sum(read_counts)
    max_depth = int(total_reads * 0.8)  # Use 80% of total reads
    
    print(f"Total reads: {total_reads}, Max depth for rarefaction: {max_depth}")
    print(f"Number of taxa: {len(read_counts)}")
    
    depths, mean_species_counts = rarefaction_curve_improved(read_counts, max_depth)

    if len(depths) == 0 or len(mean_species_counts) == 0:
        print(f"⚠️ No valid data for rarefaction curve in {sample_name}")
        return

    # Plot the curve
    plt.plot(depths, mean_species_counts, 'b-', linewidth=2, label=sample_name)
    plt.scatter(depths[::max(1, len(depths)//20)], 
               mean_species_counts[::max(1, len(depths)//20)], 
               color='red', s=30, alpha=0.7)

    plt.title(f'Rarefaction Curve - {sample_name}', fontsize=14)
    plt.xlabel('Sequencing Depth', fontsize=12)
    plt.ylabel(f'Unique Taxa ({taxonomic_level.upper()} level)', fontsize=12)
    plt.legend(title="Sample", loc='lower right')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, format='jpg', dpi=300, bbox_inches='tight')
    print(f"✅ Rarefaction plot saved to {output_file}")
    print(f"📊 Final diversity estimate: {mean_species_counts[-1]:.1f} taxa")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate and save rarefaction curves from multiple Kraken reports.')
    parser.add_argument('--input', type=str, required=True, help='File containing paths to Kraken report files')
    parser.add_argument('--taxlevel', type=str, required=True, help='Taxonomic level to analyze (e.g., S for species, F for family)')
    parser.add_argument('--output', type=str, required=True, help='Output filename for the plot (JPEG format)')

    args = parser.parse_args()
    main(args.input, args.taxlevel, args.output)

#Run individually
#python3 rarefaction_curve.py --input '/home/.../kraken_report_paths.txt' --taxlevel 'C' --output "plot.jpg"

