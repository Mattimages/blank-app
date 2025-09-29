"""
Phosphorylation Analysis Page
==========================
Database-driven analysis of protein phosphorylation patterns in M. marinum.
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

# Data paths
DATA_DIR = Path(__file__).parent.parent / "data"

def load_phosphorylation_data():
    """Load phosphorylation data from various sources."""
    phospho_files = {
        "Database 1": DATA_DIR / "phosphorylation_db1.csv",
        "Database 2": DATA_DIR / "phosphorylation_db2.csv",
        "Experimental": DATA_DIR / "experimental_phospho.csv"
    }
    
    data = {}
    for source, filepath in phospho_files.items():
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

def analyze_phosphorylation_patterns(data):
    """Analyze phosphorylation patterns across datasets."""
    if not data:
        st.warning("No phosphorylation data available for analysis.")
        return
    
    # Combine all datasets
    all_data = []
    for source, df in data.items():
        df_copy = df.copy()
        df_copy['source'] = source
        all_data.append(df_copy)
    
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Analysis tabs
    tab1, tab2, tab3, tab4 = st.tabs(["Overview", "Site Analysis", "Protein Distribution", "Functional Analysis"])
    
    with tab1:
        st.subheader("Dataset Overview")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Total Phosphosites", len(combined_df))
        with col2:
            st.metric("Unique Proteins", combined_df['protein_id'].nunique() if 'protein_id' in combined_df.columns else "N/A")
        with col3:
            st.metric("Data Sources", len(data))
        
        # Distribution by source
        if len(data) > 1:
            source_counts = combined_df['source'].value_counts()
            fig = px.pie(values=source_counts.values, names=source_counts.index, 
                        title="Phosphorylation Sites by Data Source")
            st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        st.subheader("Phosphorylation Site Analysis")
        
        # Amino acid distribution
        if 'amino_acid' in combined_df.columns:
            aa_counts = combined_df['amino_acid'].value_counts()
            fig = px.bar(x=aa_counts.index, y=aa_counts.values,
                        title="Distribution of Phosphorylated Amino Acids",
                        labels={'x': 'Amino Acid', 'y': 'Count'})
            st.plotly_chart(fig, use_container_width=True)
        
        # Position analysis
        if 'position' in combined_df.columns:
            st.subheader("Position Distribution")
            fig = px.histogram(combined_df, x='position', bins=50,
                             title="Distribution of Phosphorylation Positions")
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        st.subheader("Protein Distribution Analysis")
        
        if 'protein_id' in combined_df.columns:
            # Proteins with most phosphorylation sites
            protein_counts = combined_df['protein_id'].value_counts().head(20)
            fig = px.bar(x=protein_counts.values, y=protein_counts.index,
                        orientation='h',
                        title="Top 20 Proteins by Phosphorylation Sites",
                        labels={'x': 'Number of Sites', 'y': 'Protein ID'})
            st.plotly_chart(fig, use_container_width=True)
    
    with tab4:
        st.subheader("Functional Analysis")
        
        # GO term analysis if available
        if 'go_term' in combined_df.columns:
            go_counts = combined_df['go_term'].value_counts().head(15)
            fig = px.bar(x=go_counts.values, y=go_counts.index,
                        orientation='h',
                        title="Top GO Terms for Phosphorylated Proteins",
                        labels={'x': 'Count', 'y': 'GO Term'})
            st.plotly_chart(fig, use_container_width=True)
        
        # Pathway analysis
        if 'pathway' in combined_df.columns:
            pathway_counts = combined_df['pathway'].value_counts().head(10)
            fig = px.pie(values=pathway_counts.values, names=pathway_counts.index,
                        title="Pathway Distribution of Phosphorylated Proteins")
            st.plotly_chart(fig, use_container_width=True)

def show_phosphorylation_page():
    """Main function to display the phosphorylation analysis page."""
    st.title("🧬 Phosphorylation Analysis")
    st.markdown("""
    Analyze protein phosphorylation patterns in *M. marinum* using integrated databases
    and experimental data.
    """)
    
    # Load data
    with st.spinner("Loading phosphorylation databases..."):
        phospho_data = load_phosphorylation_data()
    
    if not phospho_data:
        st.error("No phosphorylation data could be loaded. Please check data files.")
        st.markdown("""
        Expected data files in the `data/` directory:
        - `phosphorylation_db1.csv`
        - `phosphorylation_db2.csv` 
        - `experimental_phospho.csv`
        
        Expected columns in CSV files:
        - `protein_id`: Protein identifier
        - `amino_acid`: Phosphorylated amino acid (S/T/Y)
        - `position`: Position in protein sequence
        - `go_term`: Associated GO terms (optional)
        - `pathway`: Associated pathways (optional)
        """)
        return
    
    # Analysis
    analyze_phosphorylation_patterns(phospho_data)
    
    # Data export
    st.subheader("Data Export")
    if st.button("Download Combined Dataset"):
        # Combine all data for export
        all_data = []
        for source, df in phospho_data.items():
            df_copy = df.copy()
            df_copy['source'] = source
            all_data.append(df_copy)
        
        combined_df = pd.concat(all_data, ignore_index=True)
        csv = combined_df.to_csv(index=False)
        
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name="combined_phosphorylation_data.csv",
            mime="text/csv"
        )

if __name__ == "__main__":
    show_phosphorylation_page()
