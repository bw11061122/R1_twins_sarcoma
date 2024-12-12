# Generate plots from PURPLE / AMBER / COBALT output
# Credits to Angus (ah39)

# Imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip
import os
import shutil

# Define specific colors for each chromosome
chromosome_colors = {
    'chr1': '#2b2d42',  # Dark slate gray
    'chr2': '#d35400',  # Darker orange
    'chr3': '#4b0082',  # Dark indigo
    'chr4': '#800000',  # Oxblood red
    'chr5': '#9467bd',  # Muted purple
    'chr6': '#1b4f72',  # Deep blue
    'chr7': '#0072b2',  # Sanger blue
    'chr8': '#e377c2',  # Soft pink
    'chr9': '#bcbd22',  # Olive green
    'chr10': '#7f7f7f', # Neutral gray
    'chr11': '#2ca02c', # Green
    'chr12': '#ff9f1c', # Rich amber
    'chr13': '#d90429', # Bright red
    'chr14': '#8d99ae', # Cool gray
    'chr15': '#8c564b', # Brown
    'chr16': '#264653', # Charcoal blue
    'chr17': '#e9c46a', # Sandy yellow
    'chr18': '#e76f51', # Coral
    'chr19': '#6c757d', # Slate gray
    'chr20': '#f4a261', # Light tan
    'chr21': '#1d3557', # Dark blue
    'chr22': '#ff6347', # Tomato red
    'chrX': '#3a0ca3',  # Royal blue
    'chrY': '#a8dadc'   # Light blue
}

# Define necessary functions 

# Function to load gzipped TSV file
def load_gzipped_tsv(file_path):
    with gzip.open(file_path, 'rt') as f:
        data = pd.read_csv(f, sep='\t')
    return data

# Function to plot AMBER BAF
def plot_amber_baf(patient_id, ax=None):
    base_dir = '/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/amber_out'
    file_path = f"{base_dir}/{patient_id}/{patient_id}.amber.baf.tsv.gz"
    data = load_gzipped_tsv(file_path)

    # Filter out invalid data
    filtered_df = data[data['tumorBAF'] != -1]

    # Convert chromosome to a categorical type with a specific order
    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    filtered_df['chromosome'] = pd.Categorical(filtered_df['chromosome'], categories=chromosomes, ordered=True)

    # Sort by chromosome and position
    filtered_df = filtered_df.sort_values(['chromosome', 'position'])

    # Create a new column for the combined x-axis positions
    filtered_df['combined_position'] = np.nan

    # Calculate the combined position for each chromosome
    start = 0
    chrom_midpoints = {}
    for chrom in chromosomes:
        chrom_data = filtered_df[filtered_df['chromosome'] == chrom]
        end = start + len(chrom_data)
        filtered_df.loc[filtered_df['chromosome'] == chrom, 'combined_position'] = range(start, end)
        chrom_midpoints[chrom] = (start + end) // 2
        start = end

    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 7))

    # Plotting tumor BAF with alpha for points
    for chrom in chromosomes:
        chrom_data = filtered_df[filtered_df['chromosome'] == chrom]
        ax.scatter(chrom_data['combined_position'], chrom_data['tumorBAF'], s=0.05, alpha=0.3, color=chromosome_colors[chrom])

    # Plotting legend without alpha
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=chromosome_colors[chrom], markersize=10, label=chrom) for chrom in chromosomes]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small', ncol=1, title="Chromosomes", frameon=False)

    ax.set_title(f'AMBER {patient_id} Tumor BAF', fontsize=24)  # Increased the font size for the title
    ax.set_ylabel('B-Allele Frequency')
    ax.set_xticks(list(chrom_midpoints.values()))
    ax.set_xticklabels(list(chrom_midpoints.keys()), rotation='vertical')

    if ax is None:
        plt.show()


def plot_cobalt_data(patient_id, ax=None):
    base_dir = '/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/cobalt_output'
    file_path = f"{base_dir}/{patient_id}/{patient_id}.cobalt.ratio.tsv"
    gz_file_path = f"{file_path}.gz"
    
    # Unzipping the file if necessary
    if os.path.exists(gz_file_path):
        with gzip.open(gz_file_path, 'rb') as f_in:
            with open(file_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    
    # Reading the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Filter out invalid data
    filtered_df = df[df['tumorReadCount'] != -1]
    
    # Convert chromosome to a categorical type with a specific order
    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    filtered_df['chromosome'] = pd.Categorical(filtered_df['chromosome'], categories=chromosomes, ordered=True)
    
    # Sort by chromosome and position
    filtered_df = filtered_df.sort_values(['chromosome', 'position'])
    
    # Create a new column for the combined x-axis positions
    filtered_df['combined_position'] = np.nan
    
    # Calculate the combined position for each chromosome
    start = 0
    chrom_midpoints = {}
    for chrom in chromosomes:
        chrom_data = filtered_df[filtered_df['chromosome'] == chrom]
        end = start + len(chrom_data)
        filtered_df.loc[filtered_df['chromosome'] == chrom, 'combined_position'] = range(start, end)
        chrom_midpoints[chrom] = (start + end) // 2
        start = end
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 7))

    # Plotting tumor GC ratios
    for chrom in chromosomes:
        chrom_data = filtered_df[filtered_df['chromosome'] == chrom]
        filtered_df.loc[filtered_df['chromosome'] == chrom, 'tumorGCRatio'] = np.clip(chrom_data['tumorGCRatio'], 0, 5)  # Clip values to a maximum of 5
        ax.scatter(chrom_data['combined_position'], filtered_df.loc[filtered_df['chromosome'] == chrom, 'tumorGCRatio'], s=1, alpha=0.6, color=chromosome_colors[chrom], label=chrom if chrom_data.size > 0 else "")
    ax.axhline(y=1, color='red', linestyle='--', linewidth=1)
    ax.set_title(f'COBALT {patient_id} GC Normalised Tumor Reads', fontsize=24)  # Increased the font size for the title
    ax.set_xlabel('')
    ax.set_ylabel('GC Ratio')
    ax.set_ylim(0, 5)  # Set y-axis limit for GC ratios
    ax.set_xticks(list(chrom_midpoints.values()))
    ax.set_xticklabels(list(chrom_midpoints.keys()), rotation='vertical')
    
    # Custom legend similar to AMBER
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=chromosome_colors[chrom], markersize=10, label=chrom) for chrom in chromosomes]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small', ncol=1, title="Chromosomes", frameon=False)

    if ax is None:
        plt.show()

def plot_allele_copy_numbers_segment(patient_id, ax=None):
    base_dir = '/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output'
    file_path_segment = f"{base_dir}/{patient_id}/without_GRIDDS/{patient_id}.purple.segment.tsv"
    file_path_purity = f"{base_dir}/{patient_id}/without_GRIDDS/{patient_id}.purple.purity.tsv"
    
    # Read the segment file
    df_segment = pd.read_csv(file_path_segment, sep='\t')
    
    # Ensure chromosomes are ordered
    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    df_segment['chromosome'] = pd.Categorical(df_segment['chromosome'], categories=chromosomes, ordered=True)
    df_segment = df_segment.sort_values(['chromosome', 'start'])
    
    # Create a new column for the combined x-axis positions
    df_segment['combined_position'] = np.nan
    start = 0
    chrom_midpoints = []
    for chrom in chromosomes:
        chrom_data = df_segment[df_segment['chromosome'] == chrom]
        end = start + len(chrom_data)
        df_segment.loc[df_segment['chromosome'] == chrom, 'combined_position'] = range(start, end)
        chrom_midpoints.append((start + end) / 2)
        start = end
    
    # Read the purity file
    df_purity = pd.read_csv(file_path_purity, sep='\t')
    purity = df_purity['purity'].iloc[0] * 100
    min_purity = df_purity['minPurity'].iloc[0] * 100
    max_purity = df_purity['maxPurity'].iloc[0] * 100
    ploidy = df_purity['ploidy'].iloc[0]

    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 7))

    # Plot Minor Allele Copy Number
    ax.scatter(df_segment['combined_position'], df_segment['minorAlleleCopyNumber'], color='#003366', s=1, label='Minor Allele Copy Number')
    
    # Plot Major Allele Copy Number
    ax.scatter(df_segment['combined_position'], df_segment['tumorCopyNumber'], color='#800020', s=1, label='Total Copy Number')
    
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Copy Number')
    title = f'Copy Number Call of {patient_id}: Purity - {purity:.2f}% ({min_purity:.2f}-{max_purity:.2f}% range); Ploidy - {ploidy:.2f}'
    ax.set_title(title, fontsize=24)  # Increased the font size for the title
    
    # Set y-axis limit
    ax.set_ylim(0, 6)
    
    # Set x-ticks to chromosome names in the middle and oriented vertically
    ax.set_xticks(chrom_midpoints)
    ax.set_xticklabels(chromosomes, rotation='vertical')
    
    ax.legend()
    if ax is None:
        plt.show()

def plot_allele_copy_numbers_smooth(patient_id, ax=None):
    base_dir = '/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/purple/purple_output'
    file_path_segment = f"{base_dir}/{patient_id}/without_GRIDDS/{patient_id}.purple.cnv.somatic.tsv"
    file_path_purity = f"{base_dir}/{patient_id}/without_GRIDDS/{patient_id}.purple.purity.tsv"
    
    # Read the segment file
    df_segment = pd.read_csv(file_path_segment, sep='\t')
    
    # Ensure chromosomes are ordered
    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    df_segment['chromosome'] = pd.Categorical(df_segment['chromosome'], categories=chromosomes, ordered=True)
    df_segment = df_segment.sort_values(['chromosome', 'start'])
    
    # Round the copy number values
    df_segment['minorAlleleCopyNumber'] = df_segment['minorAlleleCopyNumber'].round()
    df_segment['copyNumber'] = df_segment['copyNumber'].round()  # Use copyNumber instead of majorAlleleCopyNumber
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 7))
    
    # Plot Minor and Total Copy Number for each chromosome
    cumulative_offset = 0
    chromosome_boundaries = []
    buffer_between_segments = 50000  # Small buffer between segments
    buffer_between_chromosomes = 1000000  # Larger buffer between chromosomes
    
    for chrom in chromosomes:
        chrom_data = df_segment[df_segment['chromosome'] == chrom]
        if not chrom_data.empty:
            chromosome_boundaries.append((chrom, cumulative_offset))
            for index, row in chrom_data.iterrows():
                start_pos = cumulative_offset + row['start']
                end_pos = cumulative_offset + row['end']
                ax.plot([start_pos, end_pos], [row['minorAlleleCopyNumber'], row['minorAlleleCopyNumber']], color='#003366', lw=4, label='Minor Allele Copy Number' if index == 0 and chrom == chromosomes[0] else "")
                ax.plot([start_pos, end_pos], [row['copyNumber'], row['copyNumber']], color='#800020', lw=4, label='Total Copy Number' if index == 0 and chrom == chromosomes[0] else "")
            cumulative_offset += chrom_data['end'].max() + buffer_between_chromosomes
    
    # Add vertical lines to mark chromosome boundaries
    for chrom, boundary in chromosome_boundaries:
        ax.axvline(x=boundary, color='grey', linestyle='--')
    
    # Read the purity file
    df_purity = pd.read_csv(file_path_purity, sep='\t')
    purity = df_purity['purity'].iloc[0] * 100
    min_purity = df_purity['minPurity'].iloc[0] * 100
    max_purity = df_purity['maxPurity'].iloc[0] * 100
    ploidy = df_purity['ploidy'].iloc[0]

    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Copy Number')
    title = f'Copy Number Call of {patient_id}: Purity - {purity:.2f}% ({min_purity:.2f}-{max_purity:.2f}% range); Ploidy - {ploidy:.2f}'
    ax.set_title(title, fontsize=24)
    
    # Set y-axis limit
    ax.set_ylim(0, 6)
    
    # Set x-ticks to chromosome names in the middle and oriented vertically
    chrom_midpoints = [(boundary + (df_segment[df_segment['chromosome'] == chrom]['end'].max() / 2)) for chrom, boundary in chromosome_boundaries]
    ax.set_xticks(chrom_midpoints)
    ax.set_xticklabels([chrom for chrom, _ in chromosome_boundaries], rotation='vertical')
    
    # Add legend
    ax.legend()
    
    if ax is None:
        plt.show()

# generate plots for all samples
samples = ["PD62341aa", "PD62341ad", "PD62341h", "PD62341n", "PD62341q", "PD62341v", "PD63383t", "PD63383u" "PD63383w", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341ap", "PD62341am", "PD62341b", "PD62341u", "PD63383ap", "PD63383aq"]

for sample in samples:
    plot_amber_baf(sample)
    plot_cobalt_data(sample)
    plot_allele_copy_numbers_segment(sample)
    plot_allele_copy_numbers_smooth(sample)
