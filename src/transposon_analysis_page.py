"""
Transposon Analysis Page
======================
Analysis of transposon distributions and effects in M. marinum genome.
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import numpy as np

# Data paths
DATA_DIR = Path(__file__).parent.parent / "data"

def load_transposon_data():
    """Load transposon annotation and analysis data."""
    transposon_files = {
        "Transposon Annotations": DATA_DIR / "transposon_annotations.csv",
        "Insertion Sites": DATA_DIR / "transposon_insertions.csv",
        "Gene Disruptions": DATA_DIR / "gene_disruptions.csv"
    }
    
    data = {}
    for source, filepath in transposon_files.items():
        if filepath.exists():
            try:
                df = pd.read_csv(filepath)
                data[source] = df
                st.success(f"Loaded {len(df)} records from {source}")
            except Exception as e:
                st.warning(f"Could not load {source}: {str(e)}")
        else:
            st.info(f"{source} file not found: {filepath}")
    
    return data

def analyze_transposon_distribution(data):
    """Analyze transposon distribution patterns."""
    if not data:
        st.warning("No transposon data available for analysis.")
        return
    
    # Analysis tabs
    tab1, tab2, tab3, tab4 = st.tabs(["Distribution", "Types", "Gene Impact", "Hotspots"])
    
    with tab1:
        st.subheader("Genomic Distribution")
        
        if "Transposon Annotations" in data:
            df = data["Transposon Annotations"]
            
            # Position distribution
            if 'position' in df.columns:
                fig = px.histogram(df, x='position', bins=100,
                                 title="Distribution of Transposons Across Genome",
                                 labels={'x': 'Genome Position', 'y': 'Count'})
                st.plotly_chart(fig, use_container_width=True)
            
            # Length distribution
            if 'length' in df.columns:
                fig = px.histogram(df, x='length', bins=50,
                                 title="Transposon Length Distribution",
                                 labels={'x': 'Length (bp)', 'y': 'Count'})
                st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        st.subheader("Transposon Types")
        
        if "Transposon Annotations" in data:
            df = data["Transposon Annotations"]
            
            # Type distribution
            if 'type' in df.columns:
                type_counts = df['type'].value_counts()
                fig = px.pie(values=type_counts.values, names=type_counts.index,
                            title="Distribution of Transposon Types")
                st.plotly_chart(fig, use_container_width=True)
            
            # Family distribution
            if 'family' in df.columns:
                family_counts = df['family'].value_counts().head(10)
                fig = px.bar(x=family_counts.index, y=family_counts.values,
                            title="Top 10 Transposon Families",
                            labels={'x': 'Family', 'y': 'Count'})
                st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.subheader("Gene Impact Analysis")
        
        if "Gene Disruptions" in data:
            df = data["Gene Disruptions"]
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.metric("Total Gene Disruptions", len(df))
            
            with col2:
                if 'gene_function' in df.columns:
                    essential_genes = df[df['gene_function'].str.contains('essential', case=False, na=False)]
                    st.metric("Essential Gene Disruptions", len(essential_genes))
            
            # Disruption types
            if 'disruption_type' in df.columns:
                disruption_counts = df['disruption_type'].value_counts()
                fig = px.bar(x=disruption_counts.index, y=disruption_counts.values,
                            title="Types of Gene Disruptions",
                            labels={'x': 'Disruption Type', 'y': 'Count'})
                st.plotly_chart(fig, use_container_width=True)
    
    with tab4:
        st.subheader("Insertion Hotspots")
        
        if "Insertion Sites" in data:
            df = data["Insertion Sites"]
            
            # Hotspot analysis
            if 'position' in df.columns:
                # Create bins for hotspot analysis
                bins = np.linspace(df['position'].min(), df['position'].max(), 50)
                df['bin'] = pd.cut(df['position'], bins=bins)
                hotspots = df['bin'].value_counts().head(10)
                
                fig = px.bar(x=range(len(hotspots)), y=hotspots.values,
                            title="Top 10 Insertion Hotspots",
                            labels={'x': 'Genomic Region', 'y': 'Insertion Count'})
                st.plotly_chart(fig, use_container_width=True)

def show_transposon_page():
    """Main function to display the transposon analysis page."""
    st.title("🧬 Transposon Analysis")
    st.markdown("""
    Analyze transposon distributions, types, and their impact on gene function 
    in the *M. marinum* genome.
    """)
    
    # Load data
    with st.spinner("Loading transposon data..."):
        transposon_data = load_transposon_data()
    
    if not transposon_data:
        st.error("No transposon data could be loaded. Please check data files.")
        st.markdown("""
        Expected data files in the `data/` directory:
        - `transposon_annotations.csv`
        - `transposon_insertions.csv`
        - `gene_disruptions.csv`
        
        Expected columns:
        - `position`: Genomic position
        - `type`: Transposon type
        - `family`: Transposon family
        - `length`: Transposon length
        - `gene_function`: Function of disrupted genes
        - `disruption_type`: Type of gene disruption
        """)
        return
    
    # Analysis
    analyze_transposon_distribution(transposon_data)
    
    # Summary statistics
    st.subheader("Summary Statistics")
    
    if transposon_data:
        total_elements = sum(len(df) for df in transposon_data.values())
        st.metric("Total Transposable Elements", total_elements)
        
        # Coverage analysis
        if "Transposon Annotations" in transposon_data:
            df = transposon_data["Transposon Annotations"]
            if 'length' in df.columns:
                total_length = df['length'].sum()
                st.metric("Total Transposon Coverage", f"{total_length:,} bp")

if __name__ == "__main__":
    show_transposon_page()
