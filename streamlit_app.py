"""
Main Application for Genomic Analysis
----------------------------------
Integrates phosphorylation analysis, transposon analysis, and genome browser functionality.
"""

import streamlit as st
from pathlib import Path
from src.phosphorylation_page import show_phosphorylation_page
from src.transposon_analysis_page import show as show_transposon_analysis
from src.genome_browser import create_genome_browser_page

# Configure the main page
st.set_page_config(page_title="Genomic Analysis", page_icon="🧬", layout="wide")

# Sidebar for navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Phosphorylation Analysis", "Transposon Analysis", "Genome Browser"])

if page == "Home":
    st.title("🧬 Genomic Analysis Platform")
    st.markdown("""
    Welcome to the Genomic Analysis Platform. Select a tool from the sidebar to begin:
    
    - **Phosphorylation Analysis**: Analyze protein phosphorylation patterns
    - **Transposon Analysis**: Study transposon distributions and effects
    - **Genome Browser**: Visualize genomic features and sequences
    """)

elif page == "Phosphorylation Analysis":
    show_phosphorylation_page()

elif page == "Transposon Analysis":
    show_transposon_analysis()

elif page == "Genome Browser":
    create_genome_browser_page()
