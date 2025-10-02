#!/usr/bin/env python3
"""
Data optimization script to combine CSV, TSV, and Excel files into optimized format.
Combines gene annotations, TradDIS data, and creates optimized file for fast loading.
"""

import pandas as pd
import numpy as np
import pickle
import os
from pathlib import Path

def load_and_process_data():
    """Load all data sources and combine into optimized structure."""
    data_dir = Path("/workspaces/blank-app/data")
    
    print("🔄 Loading data files...")
    
    # Load gene annotations (TSV)
    print("📊 Loading gene annotations...")
    gene_data = pd.read_csv(data_dir / "ncbi_dataset.tsv", sep='\t')
    print(f"   Loaded {len(gene_data):,} gene annotations")
    
    # Load TradDIS data (Excel)
    print("🧪 Loading TradDIS data...")
    tradis_data = pd.read_excel(data_dir / "Tradis_TAsites_coverage_in_all_36samples_june2013_NoRmDup intergenic.xlsx")
    print(f"   Loaded {len(tradis_data):,} TradDIS entries")
    
    # Load supplementary data (CSV)
    print("📋 Loading supplementary data...")
    try:
        supp_data = pd.read_csv(data_dir / "Big_data_sup2.csv")
        print(f"   Loaded {len(supp_data):,} supplementary entries")
    except Exception as e:
        print(f"   Warning: Could not load supplementary data: {e}")
        supp_data = None
    
    return gene_data, tradis_data, supp_data

def process_tradis_data(tradis_df):
    """Process TradDIS data to extract statistics and p-values."""
    print("🔬 Processing TradDIS statistics...")
    
    # Get position column (assuming first column)
    position_col = tradis_df.columns[0]
    
    # Get viability columns I-N (indices 8-13)
    if len(tradis_df.columns) >= 14:
        viability_cols = tradis_df.columns[8:14]
        print(f"   Using viability columns: {list(viability_cols)}")
        
        processed_sites = []
        
        for idx, row in tradis_df.iterrows():
            try:
                if pd.notna(row[position_col]):
                    position = int(float(row[position_col]))
                    
                    # Extract viability values
                    values = []
                    for col in viability_cols:
                        if pd.notna(row[col]):
                            try:
                                val = float(row[col])
                                if val >= 0:
                                    values.append(val)
                            except (ValueError, TypeError):
                                continue
                    
                    if len(values) >= 3:  # Minimum for statistics
                        mean_val = np.mean(values)
                        std_val = np.std(values)
                        
                        # Calculate p-value using t-test against null hypothesis (mean = 0)
                        if std_val > 0:
                            from scipy import stats
                            t_stat, p_val = stats.ttest_1samp(values, 0)
                            p_val = abs(p_val)  # Use absolute p-value
                        else:
                            p_val = 1.0  # No variance = no significance
                        
                        # Look for ortholog information
                        ortholog_info = {}
                        for col in tradis_df.columns:
                            col_name = str(col).lower()
                            if any(keyword in col_name for keyword in ['mmar', 'ortholog', 'homolog', 'rv']):
                                if pd.notna(row[col]) and str(row[col]).strip():
                                    ortholog_info[col] = str(row[col]).strip()
                        
                        processed_sites.append({
                            'position': position,
                            'mean': mean_val,
                            'std': std_val,
                            'p_value': p_val,
                            'values': values,
                            'count': len(values),
                            'ortholog_info': ortholog_info
                        })
                        
            except (ValueError, TypeError) as e:
                continue
        
        print(f"   Processed {len(processed_sites):,} TA sites with statistics")
        return processed_sites
    
    return []

def combine_gene_and_ortholog_data(gene_data, processed_sites):
    """Combine gene data with ortholog information from TradDIS."""
    print("🧬 Combining gene and ortholog data...")
    
    # Create ortholog lookup by position
    ortholog_lookup = {}
    for site in processed_sites:
        if site['ortholog_info']:
            ortholog_lookup[site['position']] = site['ortholog_info']
    
    # Add ortholog information to genes
    enhanced_genes = []
    for idx, gene in gene_data.iterrows():
        gene_dict = gene.to_dict()
        
        # Check if gene position has ortholog data
        gene_start = gene['Begin']
        gene_end = gene['End']
        
        # Look for ortholog data in gene region
        gene_orthologs = {}
        for pos in range(gene_start, gene_end + 1, 1000):  # Sample every 1kb
            if pos in ortholog_lookup:
                gene_orthologs.update(ortholog_lookup[pos])
        
        gene_dict['ortholog_data'] = gene_orthologs
        enhanced_genes.append(gene_dict)
    
    print(f"   Enhanced {len([g for g in enhanced_genes if g['ortholog_data']]):,} genes with ortholog data")
    return enhanced_genes

def save_optimized_data(enhanced_genes, processed_sites, output_dir):
    """Save data in optimized formats."""
    print("💾 Saving optimized data...")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Save as pickle for fastest loading
    optimized_data = {
        'genes': enhanced_genes,
        'ta_sites': processed_sites,
        'metadata': {
            'gene_count': len(enhanced_genes),
            'ta_site_count': len(processed_sites),
            'created': pd.Timestamp.now().isoformat()
        }
    }
    
    pickle_file = output_dir / "optimized_genome_data.pkl"
    with open(pickle_file, 'wb') as f:
        pickle.dump(optimized_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"   Saved pickle file: {pickle_file}")
    
    # Also save as parquet for cross-platform compatibility
    genes_df = pd.DataFrame(enhanced_genes)
    ta_sites_df = pd.DataFrame(processed_sites)
    
    genes_parquet = output_dir / "genes_enhanced.parquet"
    ta_sites_parquet = output_dir / "ta_sites_processed.parquet"
    
    genes_df.to_parquet(genes_parquet, compression='snappy')
    ta_sites_df.to_parquet(ta_sites_parquet, compression='snappy')
    
    print(f"   Saved genes parquet: {genes_parquet}")
    print(f"   Saved TA sites parquet: {ta_sites_parquet}")
    
    # Save summary statistics
    stats = {
        'total_genes': len(enhanced_genes),
        'genes_with_orthologs': len([g for g in enhanced_genes if g['ortholog_data']]),
        'total_ta_sites': len(processed_sites),
        'significant_sites': len([s for s in processed_sites if s['p_value'] < 0.05]),
        'mean_coverage': np.mean([s['mean'] for s in processed_sites]),
        'file_sizes': {
            'pickle_mb': pickle_file.stat().st_size / (1024*1024),
            'genes_parquet_mb': genes_parquet.stat().st_size / (1024*1024),
            'ta_parquet_mb': ta_sites_parquet.stat().st_size / (1024*1024)
        }
    }
    
    stats_file = output_dir / "optimization_stats.json"
    import json
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"   Optimization complete!")
    print(f"   📊 Statistics saved to: {stats_file}")
    
    return stats

if __name__ == "__main__":
    # Install required packages
    try:
        from scipy import stats
    except ImportError:
        print("Installing required packages...")
        import subprocess
        subprocess.check_call(["pip", "install", "scipy", "pyarrow"])
        from scipy import stats
    
    # Process data
    gene_data, tradis_data, supp_data = load_and_process_data()
    processed_sites = process_tradis_data(tradis_data)
    enhanced_genes = combine_gene_and_ortholog_data(gene_data, processed_sites)
    
    # Save optimized data
    stats = save_optimized_data(enhanced_genes, processed_sites, "/workspaces/blank-app/data")
    
    print("\n🎉 Data optimization complete!")
    print(f"   Total genes: {stats['total_genes']:,}")
    print(f"   Genes with orthologs: {stats['genes_with_orthologs']:,}")
    print(f"   Total TA sites: {stats['total_ta_sites']:,}")
    print(f"   Significant sites (p<0.05): {stats['significant_sites']:,}")
    print(f"   Pickle file size: {stats['file_sizes']['pickle_mb']:.1f} MB")