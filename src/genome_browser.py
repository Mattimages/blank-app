#!/usr/bin/env python3
"""Interactive Genome Browser with scroll-based navigation."""

import logging
from pathlib import Path
import streamlit as st
import pandas as pd
from Bio import SeqIO
import plotly.graph_objects as go
from Bio.Seq import Seq

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Data paths
DATA_DIR = Path(__file__).parent.parent / "data"
GENOME_FILE = DATA_DIR / "Marinum_e11_genome.fasta"
ANNOTATION_FILE = DATA_DIR / "ncbi_dataset.tsv"

def load_genome():
    """Load the M. marinum genome sequence."""
    try:
        return next(SeqIO.parse(GENOME_FILE, "fasta"))
    except Exception as e:
        logger.error(f"Error loading genome sequence: {e}")
        st.error(f"Fout bij het laden van het genoom: {str(e)}")
        return None

def load_annotations():
    """Load gene annotations from TSV file."""
    try:
        df = pd.read_csv(ANNOTATION_FILE, sep="\t")
        required_columns = ['Begin', 'End', 'Orientation', 'Name']
        if not all(col in df.columns for col in required_columns):
            st.error(f"Vereiste kolommen ontbreken. Gevonden: {list(df.columns)}")
            return None
        
        # Convert Orientation to standard strand notation
        df['strand'] = df['Orientation'].map({'plus': '+', 'minus': '-'})
        # Use Name as product description
        df['product'] = df['Name'].fillna(df.get('Symbol', 'Unknown'))
        
        return df
    except Exception as e:
        logger.error(f"Error loading annotations: {e}")
        st.error(f"Fout bij het laden van annotaties: {str(e)}")
        return None

def create_arrow_shape(start_x, end_x, y, strand, height=0.3):
    """Create coordinates for an arrow shape."""
    arrow_head_size = min((end_x - start_x) * 0.2, 100)
    
    if strand == "+":
        x = [start_x, end_x - arrow_head_size, end_x - arrow_head_size, end_x, 
             end_x - arrow_head_size, end_x - arrow_head_size, start_x, start_x]
        y_coords = [y - height/2, y - height/2, y - height, y, y + height, 
                   y + height/2, y + height/2, y - height/2]
    else:
        x = [start_x, start_x + arrow_head_size, start_x + arrow_head_size, start_x, 
             start_x + arrow_head_size, start_x + arrow_head_size, end_x, end_x]
        y_coords = [y, y + height, y + height/2, y, y - height/2, y - height, 
                   y - height/2, y + height/2]
    
    return x, y_coords

def get_codon_colors():
    """Define colors for different codon types."""
    return {
        'start': '#00FF00',  # Green for start codons
        'stop': '#FF0000',   # Red for stop codons
        'A': '#FFB6C1',      # Light pink for A
        'T': '#87CEEB',      # Sky blue for T
        'G': '#98FB98',      # Pale green for G
        'C': '#DDA0DD'       # Plum for C
    }



def translate_coding_sequence(sequence, gene_start, gene_end, strand):
    """Translate a gene's coding sequence."""
    try:
        gene_seq = sequence[gene_start-1:gene_end]  # Convert to 0-based indexing
        
        if strand == '-':
            gene_seq = gene_seq.reverse_complement()
        
        # Find the best reading frame (look for start codon)
        best_translation = None
        best_frame = 0
        
        for frame in range(3):
            frame_seq = gene_seq[frame:]
            if len(frame_seq) >= 3:
                protein = frame_seq.translate()
                # Prefer translations that start with M
                if str(protein).startswith('M'):
                    best_translation = protein
                    best_frame = frame
                    break
                elif best_translation is None:
                    best_translation = protein
                    best_frame = frame
        
        return best_translation, best_frame
    except Exception as e:
        logger.error(f"Translation error: {e}")
        return None, 0

def create_interactive_genome_browser(sequence, annotations, start_pos=0, window_size=50000, show_codons=False, show_translations=False, reading_frame=0):
    """Create a fully interactive genome browser with scroll-based navigation."""
    sequence_length = len(sequence)
    end_pos = min(start_pos + window_size, sequence_length)
    
    fig = go.Figure()
    
    # Add genome baseline
    fig.add_trace(go.Scatter(
        x=[start_pos, end_pos],
        y=[0, 0],
        mode='lines',
        line=dict(color='black', width=2),
        name='Genome',
        showlegend=False,
        hoverinfo='skip'
    ))
    
    # Add DNA sequence visualization if enabled and window is small enough
    if show_codons and window_size <= 300:
        visible_seq = sequence[start_pos:end_pos]
        
        if reading_frame > 0 and reading_frame < len(visible_seq):
            visible_seq = visible_seq[reading_frame:]
            codon_start_pos = start_pos + reading_frame
        else:
            codon_start_pos = start_pos
        
        colors = get_codon_colors()
        
        # Show individual bases
        for i, base in enumerate(str(visible_seq)):
            x_pos = codon_start_pos + i
            
            # Color by base type
            base_color = colors.get(base.upper(), '#CCCCCC')
            
            # Add base as text on center line
            fig.add_trace(go.Scatter(
                x=[x_pos],
                y=[0.1],
                mode='text',
                text=[base],
                textfont=dict(size=12, color='black'),
                showlegend=False,
                hovertemplate=f"Base: {base}<br>Position: {x_pos}<extra></extra>",
                name=f"Base_{i}"
            ))
            
            # Add colored background for base
            fig.add_trace(go.Scatter(
                x=[x_pos-0.4, x_pos+0.4, x_pos+0.4, x_pos-0.4, x_pos-0.4],
                y=[0.05, 0.05, 0.15, 0.15, 0.05],
                fill='toself',
                fillcolor=base_color,
                line=dict(width=0),
                showlegend=False,
                hoverinfo='skip',
                opacity=0.7
            ))
        
        # Show codons if sequence is long enough
        if len(visible_seq) >= 3:
            codons = [str(visible_seq[i:i+3]) for i in range(0, len(visible_seq)-2, 3)]
            codon_positions = [codon_start_pos + i*3 for i in range(len(codons))]
            
            for codon, pos in zip(codons, codon_positions):
                if len(codon) == 3:
                    # Determine codon type
                    if codon.upper() in ['ATG']:
                        codon_color = colors['start']
                        codon_type = 'Start'
                    elif codon.upper() in ['TAA', 'TAG', 'TGA']:
                        codon_color = colors['stop']
                        codon_type = 'Stop'
                    else:
                        codon_color = '#E6E6FA'  # Light lavender for regular codons
                        codon_type = 'Coding'
                    
                    # Add codon bracket above sequence
                    fig.add_trace(go.Scatter(
                        x=[pos, pos+3],
                        y=[0.25, 0.25],
                        mode='lines',
                        line=dict(color=codon_color, width=3),
                        showlegend=False,
                        hovertemplate=f"Codon: {codon}<br>Type: {codon_type}<br>Position: {pos}<extra></extra>",
                        name=f"Codon_{pos}"
                    ))
    
    # Add gene annotations
    if annotations is not None:
        # Filter annotations in visible range
        mask = (annotations["Begin"] < end_pos) & (annotations["End"] > start_pos)
        visible_genes = annotations[mask]
        
        for _, gene in visible_genes.iterrows():
            gene_start = max(gene["Begin"], start_pos)
            gene_end = min(gene["End"], end_pos)
            strand = gene["strand"]
            product = gene["product"]
            
            if strand == "+":
                y_pos = 0.5
                color = "blue"
            else:
                y_pos = -0.5
                color = "red"
            
            arrow_x, arrow_y = create_arrow_shape(gene_start, gene_end, y_pos, strand)
            
            fig.add_trace(go.Scatter(
                x=arrow_x,
                y=arrow_y,
                fill="toself",
                fillcolor=color,
                line=dict(color=color, width=0.5),
                mode="lines",
                name=product,
                showlegend=False,
                hovertemplate=(
                    f"<b>{product}</b><br>" +
                    f"Position: {gene['Begin']:,} - {gene['End']:,}<br>" +
                    f"Length: {gene['End'] - gene['Begin']:,} bp<br>" +
                    f"Strand: {strand}<br>" +
                    "<extra></extra>"
                )
            ))
            
            # Add translation if enabled and window is small
            if show_translations and window_size <= 1000:
                protein, best_frame = translate_coding_sequence(
                    sequence, gene["Begin"], gene["End"], gene["strand"]
                )
                
                if protein and len(str(protein)) > 0:
                    # Show first few amino acids
                    protein_str = str(protein)[:20] + ("..." if len(str(protein)) > 20 else "")
                    
                    # Position translation text
                    text_y = y_pos + (0.2 if strand == "+" else -0.2)
                    text_x = (gene_start + gene_end) / 2
                    
                    fig.add_trace(go.Scatter(
                        x=[text_x],
                        y=[text_y],
                        mode='text',
                        text=[protein_str],
                        textfont=dict(size=8, color='darkgreen'),
                        showlegend=False,
                        hovertemplate=f"Protein: {str(protein)}<br>Length: {len(protein)} aa<extra></extra>",
                        name=f"Protein_{gene['Begin']}"
                    ))
    
    # Configure layout with proper zoom controls
    fig.update_layout(
        title=f"🧬 M. marinum Genome Browser (Position {start_pos:,} - {end_pos:,})",
        xaxis=dict(
            title="Genome Position (bp)",
            range=[start_pos, end_pos],
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            tickformat=',.0f'
        ),
        yaxis=dict(
            title="",
            range=[-1.2, 1.2],
            showticklabels=False,
            showgrid=False,
            zeroline=True,
            zerolinecolor='black',
            zerolinewidth=2
        ),
        height=500,
        plot_bgcolor="white",
        dragmode='pan',
        showlegend=False,
        hovermode='closest',
        # Enable proper zooming
        xaxis_rangeslider_visible=False,
        uirevision='constant'  # Fixed: was layout_uirevision
    )
    
    # Update config for better scroll interaction
    fig.update_layout(
        xaxis=dict(
            autorange=False,
            fixedrange=False
        ),
        yaxis=dict(
            fixedrange=True  # Lock Y-axis, only allow X-axis zoom/pan
        )
    )
    
    return fig

def create_genome_browser_page():
    """Create the main genome browser interface."""
    st.set_page_config(
        page_title="M. marinum Genome Browser",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 M. marinum Genome Browser")
    
    # Initialize session state
    if 'show_codons' not in st.session_state:
        st.session_state.show_codons = False
    if 'show_translations' not in st.session_state:
        st.session_state.show_translations = False
    if 'selected_gene' not in st.session_state:
        st.session_state.selected_gene = None
    if 'current_start' not in st.session_state:
        st.session_state.current_start = 0
    if 'current_window' not in st.session_state:
        st.session_state.current_window = 50000
    if 'reading_frame' not in st.session_state:
        st.session_state.reading_frame = 0
    
    sequence_record = load_genome()
    if not sequence_record:
        return
    
    sequence = sequence_record.seq
    sequence_length = len(sequence)
    
    annotations = load_annotations()
    if annotations is None:
        st.warning("Geen annotaties geladen - alleen genoomsequentie beschikbaar")
    
    # Navigation and View Controls
    st.subheader("Navigation & Display Controls")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.session_state.current_start = st.number_input(
            "Start Position", 
            min_value=0, 
            max_value=sequence_length-1000,
            value=st.session_state.current_start,
            step=1000
        )
    
    with col2:
        st.session_state.current_window = st.selectbox(
            "Window Size",
            [1000, 5000, 10000, 50000, 100000],
            index=3
        )
    
    with col3:
        if st.button("🧬 Toggle Bases"):
            st.session_state.show_codons = not st.session_state.show_codons
            st.rerun()
    
    with col4:
        if st.button("🔤 Toggle Translations"):
            st.session_state.show_translations = not st.session_state.show_translations
            st.rerun()
    
    # Additional controls when base view is enabled
    if st.session_state.show_codons:
        col1, col2, col3 = st.columns(3)
        with col1:
            st.session_state.reading_frame = st.selectbox("Reading Frame", [0, 1, 2], index=st.session_state.reading_frame)
        with col2:
            st.info(f"Base View: {'ON' if st.session_state.show_codons else 'OFF'}")
        with col3:
            st.info(f"Translations: {'ON' if st.session_state.show_translations else 'OFF'}")
    
    # Navigation buttons
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        if st.button("⬅️⬅️ Jump Back"):
            st.session_state.current_start = max(0, st.session_state.current_start - st.session_state.current_window)
            st.rerun()
    
    with col2:
        if st.button("⬅️ Step Back"):
            st.session_state.current_start = max(0, st.session_state.current_start - st.session_state.current_window//4)
            st.rerun()
    
    with col3:
        if st.button("Step Forward ➡️"):
            max_start = sequence_length - st.session_state.current_window
            st.session_state.current_start = min(max_start, st.session_state.current_start + st.session_state.current_window//4)
            st.rerun()
    
    with col4:
        if st.button("Jump Forward ➡️➡️"):
            max_start = sequence_length - st.session_state.current_window
            st.session_state.current_start = min(max_start, st.session_state.current_start + st.session_state.current_window)
            st.rerun()
    
    # Main genome browser with integrated views
    fig = create_interactive_genome_browser(
        sequence, 
        annotations, 
        start_pos=st.session_state.current_start,
        window_size=st.session_state.current_window,
        show_codons=st.session_state.show_codons,
        show_translations=st.session_state.show_translations,
        reading_frame=st.session_state.reading_frame
    )
    
    # Configure plot for proper interaction
    config = {
        'scrollZoom': True,
        'displayModeBar': True,
        'displaylogo': False,
        'modeBarButtonsToRemove': ['select2d', 'lasso2d', 'autoScale2d'],
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'genome_browser',
            'height': 500,
            'width': 1200,
            'scale': 1
        }
    }
    
    st.plotly_chart(fig, use_container_width=True, config=config)
    
    # Gene translations panel
    if st.session_state.show_translations and annotations is not None:
        st.subheader("Gene Translations")
        
        # Gene selection
        gene_names = annotations['product'].dropna().unique()[:20]  # Limit to first 20 genes
        selected_gene = st.selectbox("Select Gene for Translation", 
                                   options=[''] + list(gene_names))
        
        if selected_gene:
            gene_info = annotations[annotations['product'] == selected_gene].iloc[0]
            
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**Gene:** {gene_info['product']}")
                st.write(f"**Position:** {gene_info['Begin']:,} - {gene_info['End']:,}")
                st.write(f"**Strand:** {gene_info['strand']}")
                st.write(f"**Length:** {gene_info['End'] - gene_info['Begin']:,} bp")
            
            with col2:
                if 'Symbol' in gene_info:
                    st.write(f"**Symbol:** {gene_info.get('Symbol', 'N/A')}")
                if 'Gene Type' in gene_info:
                    st.write(f"**Type:** {gene_info.get('Gene Type', 'N/A')}")
                if 'Locus tag' in gene_info:
                    st.write(f"**Locus:** {gene_info.get('Locus tag', 'N/A')}")
            
            # Translate the gene
            protein, best_frame = translate_coding_sequence(
                sequence, gene_info['Begin'], gene_info['End'], gene_info['strand']
            )
            
            if protein:
                st.write(f"**Translation (Frame {best_frame+1}):**")
                st.code(str(protein), language=None)
                
                # Protein stats
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Protein Length", f"{len(protein)} aa")
                with col2:
                    st.metric("Molecular Weight", f"{len(protein) * 110:.1f} Da")
                with col3:
                    stop_codons = str(protein).count('*')
                    st.metric("Stop Codons", stop_codons)
    
    with st.sidebar:
        st.write("### Genoom Informatie")
        st.write(f"**Totale Lengte:** {sequence_length:,} bp")
        
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = (gc_count / sequence_length) * 100
        st.write(f"**GC Gehalte:** {gc_content:.1f}%")
        
        if annotations is not None:
            st.write(f"**Totaal Genen:** {len(annotations):,}")
            
            # Gene type distribution
            if 'Gene Type' in annotations.columns:
                gene_types = annotations['Gene Type'].value_counts()
                st.write("**Gen Types:**")
                for gene_type, count in gene_types.head(5).items():
                    st.write(f"- {gene_type}: {count}")
        
        st.write("### Navigatie")
        st.write("🔍 **Zoomen:** Muiswiel op plot")
        st.write("⬅️➡️ **Pannen:** Sleep met muis")
        st.write("🎯 **Reset:** Dubbelklik op plot")
        st.write("⏭️ **Navigeer:** Gebruik knoppen boven plot")
        
        # Current view info
        st.write("### Huidige Weergave")
        st.write(f"**Positie:** {st.session_state.current_start:,} - {st.session_state.current_start + st.session_state.current_window:,}")
        st.write(f"**Venster:** {st.session_state.current_window:,} bp")
        
        if st.session_state.show_codons:
            st.write(f"**Reading Frame:** {st.session_state.reading_frame + 1}")
        
        st.write("### Weergave Opties")
        st.write(f"🧬 **Bases:** {'ON' if st.session_state.show_codons else 'OFF'}")
        st.write(f"🔤 **Translations:** {'ON' if st.session_state.show_translations else 'OFF'}")
        
        # Legend for base colors
        if st.session_state.show_codons:
            st.write("### Base/Codon Kleuren")
            st.write("🟢 Start codon (ATG)")
            st.write("🔴 Stop codon (TAA/TAG/TGA)")
            st.write("🌸 Adenine (A)")
            st.write("🔷 Thymine (T)")  
            st.write("🍃 Guanine (G)")
            st.write("🟣 Cytosine (C)")
            
        # Quick navigation
        st.write("### Snelle Navigatie")
        gene_position = st.number_input(
            "Ga naar gen positie:",
            min_value=1,
            max_value=sequence_length,
            value=st.session_state.current_start + 1,
            step=1
        )
        if st.button("Ga Naar Positie"):
            st.session_state.current_start = max(0, gene_position - 1)
            st.rerun()

if __name__ == "__main__":
    create_genome_browser_page()
