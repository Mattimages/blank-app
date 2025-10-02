#!/usr/bin/env python3
"""Interactive Genome Browser with scroll-based navigation."""

import logging
from pathlib import Path
import streamlit as st

# Zorg dat navigatie helper beschikbaar is voordat UI opgebouwd wordt
def _jump_to(pos: int):
    if 'sequence_length' not in st.session_state:
        return
    genome_len = st.session_state.sequence_length
    pos = max(0, min(genome_len - 100, int(pos)))
    if st.session_state.get('current_start') != pos:
        st.session_state.current_start = pos
        st.session_state.cached_plot = None

def _handle_jump_query(query: str):
    """Parseer een gen- of locus query en update start/window.
    Ondersteunt:
      - "123456" -> center rond positie
      - "123000-127000" -> venster exact die range (afgekapt tot max 30k)
      - genenaam prefix (case-insensitive) zoekt in preloaded all_genes op 'Name' of 'Gene ID'
    """
    if not query:
        return
    seq_len = st.session_state.get('sequence_length') or 0
    genes = st.session_state.get('all_genes') or []
    q = query.strip()
    # Range pattern
    import re
    m = re.match(r"^(\d+)-(\d+)$", q)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        if a > b:
            a, b = b, a
        a = max(0, min(seq_len-100, a))
        b = max(0, min(seq_len, b))
        span = b - a
        if span <= 0:
            return
        # Beperk extreem grote ranges voor performance (nu 1M)
        if span > 1_000_000:
            b = a + 1_000_000
            span = 1_000_000
        st.session_state.current_start = a
        st.session_state.current_window = span
        st.session_state.cached_plot = None
        return
    # Single positie
    if q.isdigit():
        pos = int(q)
        pos = max(0, min(seq_len-100, pos))
        win = st.session_state.get('current_window', 5_000)
        half = win // 2
        start = max(0, min(seq_len - win, pos - half))
        st.session_state.current_start = start
        st.session_state.cached_plot = None
        return
    # Genenaam/prefix zoeken
    qlow = q.lower()
    match = None
    for g in genes:
        name = str(g.get('Name') or '').lower()
        gid = str(g.get('Gene ID') or '').lower()
        if name.startswith(qlow) or (gid and gid.startswith(qlow)):
            match = g
            break
    if match:
        gstart = match.get('Begin', match.get('start_position',0))
        gend = match.get('End', match.get('end_position', gstart+100))
        center = (gstart + gend)//2
        win = st.session_state.get('current_window', 5_000)
        half = win // 2
        start = max(0, min(seq_len - win, center - half))
        st.session_state.current_start = start
        st.session_state.cached_plot = None
        st.toast(f"Sprong naar gen: {match.get('Name') or match.get('Gene ID')}")
    else:
        st.warning("Geen match gevonden voor query")
import streamlit.components.v1 as components
from Bio import SeqIO
import plotly.graph_objects as go
from Bio.Seq import Seq
import pickle
import os
import numpy as np
from functools import lru_cache
import math

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Data paths
DATA_DIR = Path(__file__).parent.parent / "data"
CACHE_DIR = DATA_DIR / "cache"
CACHE_DIR.mkdir(exist_ok=True)

@lru_cache(maxsize=32)
def get_cached_gene_data(start, end):
    """Cache gene data for performance."""
    cache_file = CACHE_DIR / f"genes_{start}_{end}.pkl"
    if cache_file.exists():
        try:
            with open(cache_file, 'rb') as f:
                return pickle.load(f)
        except:
            pass

    return None

def save_gene_data_cache(start, end, data):
    """Save gene data to cache."""
    try:
        cache_file = CACHE_DIR / f"genes_{start}_{end}.pkl"
        with open(cache_file, 'wb') as f:
            pickle.dump(data, f)
    except:
        pass
        
# Only FASTA file needed - all other data comes from database
GENOME_FILE = DATA_DIR / "Marinum_e11_genome.fasta"

def load_genome():
    """Load the M. marinum genome sequence."""
    try:
        return next(SeqIO.parse(GENOME_FILE, "fasta"))
    except Exception as e:
        logger.error(f"Error loading genome sequence: {e}")
        st.error(f"Fout bij het laden van het genoom: {str(e)}")
        return None

# load_annotations function removed - use load_optimized_ta_data() directly

def get_genes_in_region(database, start_pos, end_pos):
    """Get genes within genomic region using SQL database.
    Normalise: always populate both 'symbol' (lower) and 'Symbol' (Title) keys for downstream code.
    Label precedence elsewhere relies on symbol -> locus -> product.
    """
    if not database:
        return []
    try:
        genes = database.get_genes_in_region(start_pos, end_pos)
        if genes:
            out = []
            for gene in genes:
                sym = gene.get('symbol') or gene.get('Symbol') or ''
                std = {
                    'Name': gene.get('name', ''),
                    'Begin': gene.get('start_position', 0),
                    'End': gene.get('end_position', 0),
                    'Strand': gene.get('strand', '+'),
                    'Product': gene.get('product', ''),
                    'Type': gene.get('gene_type', 'CDS'),
                    'Organism': gene.get('organism', ''),
                    'Locus tag': gene.get('locus_tag', ''),
                    'Symbol': sym,  # canonical display key
                    'symbol': sym,  # lower-case convenience
                    'Gene Type': gene.get('gene_type', 'CDS'),
                    'Orientation': 'plus' if gene.get('strand', '+') == '+' else 'minus',
                    'strand': gene.get('strand', '+'),
                    'product': gene.get('product', ''),
                    'gene_name': gene.get('name', ''),
                    'locus_tag': gene.get('locus_tag', ''),
                    'ortholog_data': gene.get('ortholog_data', {}),
                    'Protein accession': gene.get('protein_accession', ''),
                    'Protein length': gene.get('protein_length', '')
                }
                out.append(std)
            return out
        return []
    except Exception as e:
        logging.error(f"Error querying genes: {e}")
        return []

def load_optimized_ta_data():
    """Initialize high-performance database and return database instance."""
    try:
        try:
            from .database_v2 import HighPerformanceGenomicDB
        except ImportError:
            from database_v2 import HighPerformanceGenomicDB
        
        # Initialize high-performance database if needed (cached after first call)
        if not hasattr(st.session_state, 'hp_database'):
            with st.spinner('⚡ Initializing high-performance genomic database...'):
                # Create new multi-database architecture instance
                db = HighPerformanceGenomicDB()
                stats = db.get_database_stats()
                
                # Cache the database instance for reuse
                st.session_state.hp_database = db
                st.session_state.database_initialized = True
                
                logging.info(f"🚀 High-Performance Database ready: {stats['ta_sites_count']:,} TA sites, {stats['genes_count']:,} genes, {stats.get('orthologs_count', 0):,} orthologs")
                st.success(f"⚡ Database loaded: {stats['ta_sites_count']:,} TA sites, {stats['genes_count']:,} genes")
                
                return db
        else:
            return st.session_state.hp_database
            
    except Exception as e:
        logging.error(f"Error loading high-performance database: {e}")
        st.error(f"❌ Database error: {e}")
        return None

    # ---- TA Gradient Aggregatie & Cache ----

        def _aggregate_ta_segments(ta_sites, start_pos, end_pos, segments):
            """Reduce raw TA sites to fixed number of segments with averaged scaled_y & p.
            Returns list of dicts: {x0,x1,y0,y1,p_mid}."""
            if not ta_sites or segments < 2:
                return []
            # Filter region
            region = [s for s in ta_sites if start_pos <= s['position'] <= end_pos]
            if not region:
                region = ta_sites
            region.sort(key=lambda s: s['position'])
            first = region[0]['position']; last = region[-1]['position']
            span = max(1, last - first)
            bins = [[] for _ in range(segments)]
            for s in region:
                idx = int((s['position'] - first) / span * (segments-1))
                idx = max(0, min(segments-1, idx))
                bins[idx].append(s)
            # Compute averages and build contiguous segment list
            points = []
            for i, bucket in enumerate(bins):
                if not bucket:
                    continue
                mean_y = float(np.mean([b.get('_scaled_y',0) for b in bucket]))
                mean_p = float(np.mean([b.get('_p',0.5) for b in bucket]))
                mid_pos = float(np.mean([b['position'] for b in bucket]))
                points.append((mid_pos, mean_y, mean_p))
            segments_out = []
            for i in range(len(points)-1):
                x0,y0,p0 = points[i]
                x1,y1,p1 = points[i+1]
                p_mid = (p0+p1)/2
                segments_out.append({'x0':x0,'x1':x1,'y0':y0,'y1':y1,'p_mid':p_mid})
            return segments_out
def _p_to_color(p):
    """Continuous gradient (adjusted):
    0.00 -> dark blue (0,70,180)
    0.05 -> yellow (255,255,0)
    0.05 - 1.00 -> yellow -> red (255,0,0)
    Values >1.0 are clamped. This centers the visual significance threshold at p=0.05 (yellow)."""
    if p is None:
        p = 0.5
    p = max(0.0, min(1.5, p))
    threshold = 0.05
    if p <= threshold:
        t = p / threshold if threshold > 0 else 1.0
        r = int(0 + t * 255)
        g = int(70 + t * (255-70))
        b = int(180 + t * (0-180))
    else:
        span = 1.0 - threshold
        t = 0 if span <= 0 else min(1.0, (p - threshold) / span)
        r = 255
        g = int(255 - t * 255)
        b = 0
    return f'rgb({r},{g},{b})'

def _p_vec_to_colors(p_arr: np.ndarray):
    """Vectorized mapping p -> rgb consistent met _p_to_color (yellow at p=0.05)."""
    if p_arr.size == 0:
        return []
    threshold = 0.05
    p = np.clip(p_arr, 0.0, 1.5)
    colors = np.empty(p.shape, dtype=object)
    low = p <= threshold
    high = ~low
    if np.any(low):
        t = np.divide(p[low], threshold, out=np.zeros_like(p[low]), where=threshold>0)
        r = (0 + t * 255).astype(int)
        g = (70 + t * (255-70)).astype(int)
        b = (180 + t * (0-180)).astype(int)
        low_idx = np.flatnonzero(low)
        for i,(ri,gi,bi) in enumerate(zip(r,g,b)):
            colors[low_idx[i]] = f'rgb({ri},{gi},{bi})'
    if np.any(high):
        span = (1.0 - threshold) if (1.0 - threshold) > 0 else 1.0
        t = np.clip((p[high]-threshold)/span, 0, 1)
        r = np.full(t.shape, 255, dtype=int)
        g = (255 - t * 255).astype(int)
        b = np.zeros(t.shape, dtype=int)
        high_idx = np.flatnonzero(high)
        for i,(ri,gi,bi) in enumerate(zip(r,g,b)):
            colors[high_idx[i]] = f'rgb({ri},{gi},{bi})'
    return colors.tolist()

    def get_cached_ta_segments(region_start, region_end, ta_sites, resolution):
        """Cache TA aggregated segments in session_state keyed by region & resolution."""
        key = f"ta_seg_{region_start}_{region_end}_{resolution}"
        cache = st.session_state.setdefault('ta_segment_cache', {})
        if key in cache:
            return cache[key]
        segs = _aggregate_ta_segments(ta_sites, region_start, region_end, resolution)
        cache[key] = segs
        # Limit cache size
        if len(cache) > 100:
            # rudimentary LRU removal (pop first inserted)
            first_key = next(iter(cache.keys()))
            if first_key != key:
                cache.pop(first_key, None)
        return segs

def get_ta_sites_in_region(database, start_pos, end_pos):
    """Get TA sites within genomic region using SQL database."""
    if not database:
        return []
    
    try:
        return database.get_ta_sites_in_region(start_pos, end_pos)
    except Exception as e:
        logger.error(f"Error querying TA sites: {e}")
        return []

def create_smooth_tradis_line(ta_sites, start_pos, end_pos, smoothing_window=2000):
    """Create smooth continuous TradDIS line for genome browser."""
    if not ta_sites:
        return [], [], []
    
    # Create position array for interpolation every 500bp for performance
    positions = np.arange(start_pos, end_pos + 1, step=500)
    heights = np.zeros(len(positions))
    p_values = np.ones(len(positions))  # Default to 1 (non-significant)
    
    # For each position, find nearby TA sites and interpolate
    for i, pos in enumerate(positions):
        nearby_sites = [
            site for site in ta_sites 
            if abs(site['position'] - pos) <= smoothing_window
        ]
        
        if nearby_sites:
            # Weight by distance (inverse distance weighting)
            weights = []
            values = []
            p_vals = []
            
            for site in nearby_sites:
                distance = abs(site['position'] - pos)
                weight = 1 / (distance + 1)  # +1 to avoid division by zero
                weights.append(weight)
                values.append(site.get('viability_mean', 0))
                p_vals.append(site.get('p_value', 0))
            
            # Weighted average
            total_weight = sum(weights)
            if total_weight > 0:
                heights[i] = sum(v * w for v, w in zip(values, weights)) / total_weight
                p_values[i] = sum(p * w for p, w in zip(p_vals, weights)) / total_weight
    
    return positions.tolist(), heights.tolist(), p_values.tolist()

def get_color_from_pvalue(p_value):
    """Convert p-value to RGB color string for smooth transitions."""
    if p_value >= 0.05:
        return 'rgb(0, 100, 255)'  # Blue (not significant)
    elif p_value >= 0.01:
        # Interpolate between blue and orange
        t = (0.05 - p_value) / (0.05 - 0.01)
        r = int(0 + t * 255)
        g = int(100 + t * 65)
        b = int(255 - t * 255)
        return f'rgb({r}, {g}, {b})'
    elif p_value >= 0.001:
        # Interpolate between orange and red
        t = (0.01 - p_value) / (0.01 - 0.001)
        r = 255
        g = int(165 - t * 165)
        b = 0
        return f'rgb({r}, {g}, {b})'
    else:
        # Red to dark red
        t = min(1.0, (0.001 - p_value) / 0.0005)
        r = int(255 - t * 116)
        g = 0
        b = 0
        return f'rgb({r}, {g}, {b})'

def create_arrow_shape(start_x, end_x, y, strand, height=0.3):
    """Create coordinates for arrows with ABSOLUTELY FIXED sizing - ZERO dynamic scaling."""
    gene_length = end_x - start_x
    
    # COMPLETELY FIXED arrow head size - NEVER CHANGES REGARDLESS OF ZOOM
    arrow_head_size = 25  # Fixed 25 pixels - DONE!
    
    # FIXED height - never changes
    body_height = height * 0.5
    point_height = body_height * 1.25
    
    if strand == "+":
        # Forward arrow (right-pointing): prominent arrowhead
        x = [start_x, end_x - arrow_head_size, end_x - arrow_head_size, end_x, 
             end_x - arrow_head_size, end_x - arrow_head_size, start_x, start_x]
        y_coords = [y - body_height, y - body_height, y - point_height, y, y + point_height, 
                   y + body_height, y + body_height, y - body_height]
    else:
        # Reverse arrow (left-pointing): prominent arrowhead  
        x = [start_x, start_x + arrow_head_size, start_x + arrow_head_size, 
             end_x, end_x, start_x + arrow_head_size, start_x + arrow_head_size, start_x]
        y_coords = [y, y + point_height, y + body_height, y + body_height, 
                   y - body_height, y - body_height, y - point_height, y]
    
    return x, y_coords

def calculate_gene_levels(genes_list, ta_analysis_active=False):
    """Calculate y-levels for genes to avoid overlap with guaranteed spacing.
    When TA analysis is active, push + strand genes below centerline."""
    if not genes_list:
        return {}
    
    # Sort genes by start position
    sorted_genes = sorted(genes_list, key=lambda g: g.get('Begin', g.get('start_position', 0)))
    gene_levels = {}
    
    # Separate forward and reverse strand genes
    forward_genes = [g for g in sorted_genes if g.get('strand', g.get('Strand', '+')).strip() == '+']
    reverse_genes = [g for g in sorted_genes if g.get('strand', g.get('Strand', '+')).strip() == '-']
    
    # Minimum spacing between genes (in bases)
    min_spacing = 100
    
    # Assign levels for forward genes (BELOW centerline, negative y values)
    forward_levels = []
    for gene in forward_genes:
        level = 0
        gene_start = gene.get('Begin', gene.get('start_position', 0))
        gene_end = gene.get('End', gene.get('end_position', 0))
        
        # Find a level where this gene doesn't overlap (with spacing)
        while level < len(forward_levels):
            overlap = False
            for existing_start, existing_end in forward_levels[level]:
                # Check for overlap with spacing buffer
                if not (gene_end + min_spacing < existing_start or gene_start > existing_end + min_spacing):
                    overlap = True
                    break
            if not overlap:
                break
            level += 1
        
        # Create new level if needed
        while level >= len(forward_levels):
            forward_levels.append([])
        
        forward_levels[level].append((gene_start, gene_end))
        gene_name = gene.get('Name', gene.get('gene_name', ''))
        
        # FIXED: + strand genes are ALWAYS BELOW centerline (y < 0)
        # TA analysis goes ABOVE centerline (y > 0) 
        if ta_analysis_active:
            # When TA analysis active, + genes go further below to make more room
            gene_levels[gene_name] = -2.0 - (level * 1.2)  # Much bigger spacing!
        else:
            # Default: + strand genes BELOW centerline, never above (y < 0)
            gene_levels[gene_name] = -1.5 - (level * 1.2)  # Much bigger spacing!
    
    # Assign levels for reverse genes (further below zero)
    reverse_levels = []
    for gene in reverse_genes:
        level = 0
        gene_start = gene.get('Begin', gene.get('start_position', 0))
        gene_end = gene.get('End', gene.get('end_position', 0))
        
        # Find a level where this gene doesn't overlap (with spacing)
        while level < len(reverse_levels):
            overlap = False
            for existing_start, existing_end in reverse_levels[level]:
                # Check for overlap with spacing buffer
                if not (gene_end + min_spacing < existing_start or gene_start > existing_end + min_spacing):
                    overlap = True
                    break
            if not overlap:
                break
            level += 1
        
        # Create new level if needed
        while level >= len(reverse_levels):
            reverse_levels.append([])
        
        reverse_levels[level].append((gene_start, gene_end))
        gene_name = gene.get('Name', gene.get('gene_name', ''))
        
        # FIXED: - strand genes are FURTHER BELOW + strand genes
        if ta_analysis_active:
            # When TA analysis active, - genes go even further below
            gene_levels[gene_name] = -4.0 - (level * 1.2)  # Much further down + bigger spacing!
        else:
            # Default: - strand genes FURTHER below + strand genes
            gene_levels[gene_name] = -3.0 - (level * 1.2)  # Much further down + bigger spacing!
    
    return gene_levels
def get_codon_colors():
    return {}

def translate_coding_sequence(sequence, gene_start, gene_end, strand):
    return None, 0

    # Configure layout with pan functionality
    fig.update_layout(
        title=f"🧬 M. marinum Genome Browser (Position {visible_start:,} - {visible_end:,})",
        xaxis=dict(
            title="Genome Position (bp)",
            range=[visible_start, visible_end],
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            tickformat=',.0f'
        ),
        yaxis=dict(
            title="",
            range=[-15.0, 28.0],  # MUCH bigger area below X-axis for genes (was -8.0) + TA sites above
            showticklabels=False,
            showgrid=False,
            zeroline=True,
            zerolinecolor='black',
            zerolinewidth=2
        ),
        height=600,
        plot_bgcolor="white",
        dragmode='pan',  # Enable pan by default
        showlegend=True,  # Enable legend for p-value colors
        hovermode='closest',
        xaxis_rangeslider_visible=False,
        uirevision='constant'
    )
    
    # Configure axes for scroll-wheel panning (clean)
    fig.update_xaxes(autorange=False)
    fig.update_yaxes(fixedrange=True, scaleanchor=None)
    return fig

# ---------------- VEREENVOUDIGDE VOLLEDIGE GENOOM WEERGAVE -----------------

# Tile/adaptieve code volledig verwijderd per laatste vereiste.
if 'fast_metrics' not in st.session_state:
    st.session_state.fast_metrics = {}

def _fast_assign_levels(genes):
    plus = [g for g in genes if g.get('strand', '+') == '+']
    minus = [g for g in genes if g.get('strand', '+') == '-']
    plus.sort(key=lambda g: g.get('Begin', g.get('start_position', 0)))
    minus.sort(key=lambda g: g.get('Begin', g.get('start_position', 0)))
    def pack(arr):
        levels = []  # list of end positions
        assigned = {}
        for g in arr:
            s = g.get('Begin', g.get('start_position', 0))
            e = g.get('End', g.get('end_position', 0))
            placed = False
            for li, lane_end in enumerate(levels):
                if s > lane_end + 50:  # 50bp buffer
                    levels[li] = e
                    assigned[id(g)] = li
                    placed = True
                    break
            if not placed:
                levels.append(e)
                assigned[id(g)] = len(levels)-1
        return assigned, len(levels)
    plus_assigned, plus_h = pack(plus)
    minus_assigned, minus_h = pack(minus)
    level_map = {}
    for g in plus:
        level_map[id(g)] = -(plus_assigned[id(g)] + 1) * 1.0
    for g in minus:
        level_map[id(g)] = -(plus_h + minus_assigned[id(g)] + 2) * 1.0
    return level_map
def _arrow_path(gene, y, window_bp, min_head=18, head_fraction=0.12):
    start_x = gene['Begin']
    end_x = gene['End']
    length = max(1, end_x - start_x)
    strand = gene.get('strand', '+')
    head_size = max(min_head, int(length * head_fraction))
    head_size = min(head_size, int(length * 0.45))  # cap head to 45% of gene
    # Increased arrow height for better visibility
    body_height = 0.585  # 30% dikker (voorheen 0.45)
    point_height = body_height * 1.8
    if strand == '+':
        shaft_end = end_x - head_size
        x = [start_x, shaft_end, shaft_end, end_x, shaft_end, shaft_end, start_x, start_x]
        yv = [y-body_height, y-body_height, y-point_height, y, y+point_height, y+body_height, y+body_height, y-body_height]
    else:
        shaft_start = start_x + head_size
        x = [start_x, shaft_start, shaft_start, end_x, end_x, shaft_start, shaft_start, start_x]
        yv = [y, y+point_height, y+body_height, y+body_height, y-body_height, y-body_height, y-point_height, y]
    return x, yv

def _choose_label(gene):
    locus = gene.get('Locus tag') or gene.get('locus_tag')
    product = gene.get('Product') or gene.get('product')
    if locus and str(locus).strip():
        lc = str(locus).strip()
        if lc.startswith('MMARE11_'):
            lc = lc.replace('MMARE11_', '')
        return lc
    if product and str(product).strip():
        return str(product).strip()[:18]
    return ''

def _truncate_label(label, max_chars):
    if len(label) <= max_chars:
        return label
    if max_chars <= 3:
        return label[:max_chars]
    return label[:max_chars-1] + '…'

GLOBAL_Y_MIN = -22  # vaste ondergrens voor gen pijlen
GLOBAL_Y_MAX = 28   # bovengrens voor geschaalde TA pieken

def create_fast_genome_browser(sequence, db, center_pos, window_size):
    # Volledige chunk (750k) renderen; zichtbare viewport limiteren via x-axis range
    region_start = st.session_state.load_region_start
    region_end = st.session_state.load_region_end
    full_len = len(sequence)
    # We tonen alle genen/TA sites in de chunk zodat panning client-side gebeurt
    visible_start = region_start
    visible_end = region_end
    # Load region genes & TA
    genes = get_genes_in_region(db, region_start, region_end)
    ta_sites = []
    if st.session_state.show_tradis_analysis:
        ta_sites = get_ta_sites_in_region(db, region_start, region_end)
        for site in ta_sites:
            mean = site.get('viability_mean', 0) or 0
            site['_coverage'] = mean
            site['_scaled_y'] = mean / 10.0  # eenvoudige schaal (cov/10)
            std_val = site.get('viability_std', 1.0)
            mean_val = mean
            if std_val == 0 and mean_val == 0:
                site['_p'] = 0.001
            elif std_val == 0 and mean_val > 0:
                site['_p'] = 0.01
            else:
                site['_p'] = site.get('p_value', 0.5)
    # Enforce primary genes only (MMARE11_RS0 prefix) if flag set
    if st.session_state.get('filter_primary_genes', True):
        filtered = []
        for g in genes:
            lt = g.get('Locus tag') or g.get('locus_tag') or ''
            if str(lt).startswith('MMARE11_RS0'):
                filtered.append(g)
        genes = filtered
    st.session_state.fast_metrics['genes_loaded'] = len(genes)
    st.session_state.fast_metrics['ta_loaded'] = len(ta_sites) if ta_sites else 0
    level_map = _fast_assign_levels(genes)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, full_len], y=[0,0], mode='lines', line=dict(color='black', width=2), hoverinfo='skip', showlegend=False))
    # Gene arrows rendering (re-added)
    if genes:
        window_bp = region_end - region_start
        for g in genes:
            g_begin = g.get('Begin', g.get('start_position', 0))
            g_end = g.get('End', g.get('end_position', 0))
            lane_y = level_map.get(id(g), -1.0)
            strand = g.get('strand', '+')
            x_coords, y_coords = _arrow_path(g, lane_y, window_bp)
            fill_color = 'steelblue' if strand == '+' else 'tomato'
            fig.add_trace(go.Scatter(
                x=x_coords, y=y_coords, mode='lines', fill='toself',
                fillcolor=fill_color, line=dict(color=fill_color, width=1),
                hoverinfo='skip', showlegend=False, name='gene'))
    if st.session_state.show_tradis_analysis and ta_sites:
        # Volledige ruwe punten
        raw_pos = np.array([s['position'] for s in ta_sites], dtype=float)
        raw_y = np.array([s.get('_scaled_y', 0) for s in ta_sites], dtype=float)
        raw_p = np.array([s.get('_p', 0.5) for s in ta_sites], dtype=float)
        order = np.argsort(raw_pos)
        raw_pos, raw_y, raw_p = raw_pos[order], raw_y[order], raw_p[order]

        # Geen filtering of downsampling: toon ALLE TA sites (kan performance beïnvloeden)
        raw_pos_vis = raw_pos
        raw_y_vis = raw_y
        raw_p_vis = raw_p
        st.session_state.fast_metrics['ta_downsample_step'] = '-'

        raw_mode = st.session_state.get('ta_raw_line_mode', False)
        # Snelle modus optie in session (kan elders via UI worden gezet, anders standaard False)
        if raw_mode:
            # Eén polyline: gebruik kleur o.b.v. mediane p zodat Plotly slechts 1 trace hoeft te tekenen
            med_p = float(np.median(raw_p)) if raw_p.size else 0.5
            fig.add_trace(go.Scattergl(
                x=raw_pos, y=raw_y, mode='lines+markers',
                line=dict(width=1.2, color=_p_to_color(med_p)),
                marker=dict(size=4, color=_p_to_color(med_p), opacity=0.7),
                hoverinfo='skip', showlegend=False, name='TA raw line'))
            # Voeg aparte marker trace toe voor hover per TA site
            marker_colors = _p_vec_to_colors(raw_p)
            fig.add_trace(go.Scattergl(
                x=raw_pos, y=raw_y, mode='markers',
                marker=dict(size=5, color=marker_colors, opacity=0.9, line=dict(width=0)),
                customdata=np.stack([raw_p, raw_y], axis=1),
                hovertemplate='Pos: %{x:,}<br>Scaled(c/10): %{y:.2f}<br>p-value: %{customdata[0]:.4f}<extra></extra>',
                showlegend=False, name='TA sites'))
        else:
            # Aggregatie per resolutie-bin: gemiddelde positie, y en p binnen elke bin
            grid_step = float(st.session_state.get('ta_line_resolution', 25))
            full_grid_start = max(region_start, raw_pos[0])
            full_grid_end = min(region_end, raw_pos[-1])
            bin_edges = np.arange(full_grid_start, full_grid_end + grid_step, grid_step)
            if bin_edges[-1] < full_grid_end:
                bin_edges = np.append(bin_edges, full_grid_end)
            # Digitize raw_pos
            idx = np.digitize(raw_pos, bin_edges) - 1  # bin index
            # Beperk idx binnen bereik
            idx = np.clip(idx, 0, len(bin_edges)-2)
            agg_x = []
            agg_y = []
            agg_p = []  # median p per bin
            # We itereren uniek over gebruikte bins voor snelheid
            for b in np.unique(idx):
                mask = idx == b
                # Neem gemiddelde (pos & y) en median (p) voor representatieve knoop
                bx = float(raw_pos[mask].mean())
                by = float(raw_y[mask].mean())
                bp = float(np.median(raw_p[mask]))
                agg_x.append(bx)
                agg_y.append(by)
                agg_p.append(bp)
            agg_x = np.array(agg_x)
            agg_y = np.array(agg_y)
            agg_p = np.array(agg_p)
            # Sorteren op positie (digitize + unique kan volgorde wisselen)
            order_bins = np.argsort(agg_x)
            agg_x = agg_x[order_bins]
            agg_y = agg_y[order_bins]
            agg_p = agg_p[order_bins]
            st.session_state.fast_metrics['ta_bins'] = agg_x.size
            # Opbouwen segmenten: kleur bepaald door median p van start-bin
            if agg_x.size > 1:
                start_colors = _p_vec_to_colors(agg_p[:-1])  # één kleur per segment op basis van start bin
                col_arr = np.array(start_colors)
                change = np.where(col_arr[1:] != col_arr[:-1])[0] + 1
                run_starts = np.concatenate(([0], change))
                run_ends = np.concatenate((change, [len(start_colors)]))
                batched = []
                for rs, re in zip(run_starts, run_ends):
                    # Segment run omvat knopen rs..re, gebruik kleur van rs
                    xs = agg_x[rs:re+1].tolist()
                    ys = agg_y[rs:re+1].tolist()
                    batched.append((xs, ys, start_colors[rs]))
                st.session_state.fast_metrics['ta_segments'] = len(start_colors)
            else:
                batched = []

        # Geen kunstmatige cap meer; gekozen resolutie wordt direct gebruikt zonder bridging

        if not raw_mode:
            # Voeg batched lijnen toe
            for xs, ys_line, col in batched:
                fig.add_trace(go.Scattergl(
                    x=xs, y=ys_line, mode='lines',
                    line=dict(width=1.4, color=col),
                    hoverinfo='skip', showlegend=False, name='TA line'))

            # Markers groter voor zichtbaarheid (alle ruwe punten)
            marker_colors = _p_vec_to_colors(raw_p)
            fig.add_trace(go.Scattergl(
                x=raw_pos, y=raw_y, mode='markers',
                marker=dict(size=5, color=marker_colors, opacity=0.9, line=dict(width=0)),
                customdata=np.stack([raw_p, raw_y], axis=1),
                hovertemplate='Pos: %{x:,}<br>Scaled(c/10): %{y:.2f}<br>p-value: %{customdata[0]:.4f}<extra></extra>',
                showlegend=False, name='TA sites'))
    # Gen pijlen + labels (unified hover)
    # Normaliseer genes lijst voor label rendering en gebruik symbol prioriteit
    text_x=[]; text_y=[]; text_meta=[]; hover_x=[]; hover_y=[]; hover_text=[]
    for g in genes:
        g_begin = g.get('Begin', g.get('start_position',0))
        g_end = g.get('End', g.get('end_position',0))
        lane_y = level_map.get(id(g), -1.0)
        strand = g.get('strand','+')
        product = (g.get('Product') or g.get('product') or '').strip()
        locus = (g.get('Locus tag') or g.get('locus_tag') or '').strip()
        # Label op pijl: locus (prefix weg) of ingekorte product fallback
        primary_label = locus or (product[:18] + ('…' if len(product) > 18 else '') if product else 'Gene')
        if primary_label.startswith('MMARE11_'):
            primary_label = primary_label.replace('MMARE11_','')
        mid = (g_begin + g_end)/2
        text_x.append(mid)
        # Centreer label exact op lane_y (geen offset)
        text_y.append(lane_y)
        text_meta.append(primary_label)
        # Hover: product als titel + locus + lengte + strand
        length_bp = g_end - g_begin
        hover_product = product if product else '(no product annotation)'
        hover_label = f"<b>{hover_product}</b><br>Locus: {locus if locus else 'NA'}<br>{g_begin:,}-{g_end:,} ({length_bp} bp) Strand: {strand}"
        hover_x.append(mid)
        hover_y.append(lane_y)
        hover_text.append(hover_label)
    if text_x:
        fig.add_trace(go.Scatter(
            x=text_x, y=text_y, mode='text', text=text_meta,
            textfont=dict(size=13, color='white'),  # ~20% groter (11 -> 13)
            hoverinfo='skip', name='gene_labels', showlegend=False,
            textposition='middle center'
        ))
    if hover_x:
        fig.add_trace(go.Scatter(
            x=hover_x, y=hover_y, mode='markers',
            marker=dict(size=8, opacity=0, color='rgba(0,0,0,0)'),
            hovertemplate='%{text}<extra></extra>',
            text=hover_text,
            showlegend=False, name='gene_info'
        ))
    # ...existing code...
    # Vaste coverage bovengrens 250 (waarden erboven buiten beeld) – voorkomen dat genpijlen verdwijnen onderin
    COVERAGE_CAP = 250
    # Herstel naar oorspronkelijke vaste range + annotation uitleg scaling
    init_span = st.session_state.get('initial_view_span', 5_000)
    init_end = min(visible_start + init_span, visible_end)
    fig.update_layout(
        xaxis=dict(range=[visible_start, init_end], showgrid=False, tickformat=',.0f'),
        yaxis=dict(range=[GLOBAL_Y_MIN, GLOBAL_Y_MAX], showticklabels=False, zeroline=True, zerolinecolor='black', fixedrange=True),
        height=550,
        plot_bgcolor='white', margin=dict(l=40,r=10,t=40,b=30), hovermode='closest', dragmode='pan',
        uirevision='genome_view_v1'
    )
    fig.add_annotation(
        x=visible_start + (visible_end-visible_start)*0.01,
        y=GLOBAL_Y_MAX - 1,
        xref='x', yref='y',
        text='Coverage scaling: y = coverage / 10  (top≈28 → ~280 coverage) – values above may clip',
        showarrow=False,
        font=dict(size=11,color='#333'),
        bgcolor='rgba(255,255,255,0.6)',
        bordercolor='#777',
        borderwidth=1,
        align='left'
    )
    st.session_state.fast_metrics.update({'traces': len(fig.data)})
    # Voeg (niet-invasieve) performance annotation (rechts boven) toe
    fm = st.session_state.fast_metrics
    perf_text = f"tr:{fm.get('traces','?')} bins:{fm.get('ta_bins','?')} seg:{fm.get('ta_segments','?')} ds:{fm.get('ta_downsample_step','?')}"
    fig.add_annotation(xref='paper', yref='paper', x=0.99, y=0.98, text=perf_text,
                       showarrow=False, font=dict(size=10,color='#555'), align='right',
                       bordercolor='#ccc', borderwidth=1, bgcolor='rgba(255,255,255,0.5)')
    return fig


def create_genome_browser_page():
    """Create the main genome browser interface."""
    st.set_page_config(
        page_title="M. marinum Genome Browser",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 M. marinum Genome Browser")

    # Read start position from query params if present (server-side sync)
    qp = st.query_params
    
   
       
            
           
def create_genome_browser_page():
    """Create the main genome browser interface."""
    st.set_page_config(
        page_title="M. marinum Genome Browser",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 M. marinum Genome Browser")

    # Read start position from query params if present (server-side sync)
    qp = st.query_params
    if 'start' in qp:
        try:
            requested_start = int(qp['start'])
            if 0 <= requested_start < st.session_state.sequence_length:
                st.session_state.current_start = requested_start
        except:
            pass
    # Metrics container (always collected silently)
    if 'fast_metrics' not in st.session_state:
        st.session_state.fast_metrics = {}
    
    # Initialize session state with caching
    if 'show_codons' not in st.session_state:
        st.session_state.show_codons = False
    if 'show_translations' not in st.session_state:
        st.session_state.show_translations = False
    if 'selected_gene' not in st.session_state:
        st.session_state.selected_gene = None
    if 'current_start' not in st.session_state:
        st.session_state.current_start = 0
    CHUNK_SIZE = 750_000  # Altijd volledige chunk laden
    if 'current_window' not in st.session_state:
        # Zichtbaar venster start op 5k, volledige 750k data wordt toch geladen
        st.session_state.current_window = 5_000
    if 'initial_view_span' not in st.session_state:
        st.session_state.initial_view_span = 5_000
    if 'reading_frame' not in st.session_state:
        st.session_state.reading_frame = 0
    if 'filter_primary_genes' not in st.session_state:
        st.session_state.filter_primary_genes = True
    if 'analysis_gene' not in st.session_state:
        st.session_state.analysis_gene = None
    if 'analysis_extension' not in st.session_state:
        st.session_state.analysis_extension = 100
    if 'show_analysis' not in st.session_state:
        st.session_state.show_analysis = False
    if 'last_plot_params' not in st.session_state:
        st.session_state.last_plot_params = None
    if 'cached_plot' not in st.session_state:
        st.session_state.cached_plot = None
    # CHUNK_SIZE reeds hierboven gedefinieerd voor gebruik bij current_window
    OVERLAP = 250_000
    if 'load_region_start' not in st.session_state:
        st.session_state.load_region_start = 0
    if 'load_region_end' not in st.session_state:
        st.session_state.load_region_end = CHUNK_SIZE
    if 'region_chunk_size' not in st.session_state:
        st.session_state.region_chunk_size = CHUNK_SIZE
    if 'region_overlap' not in st.session_state:
        st.session_state.region_overlap = OVERLAP
    if 'position_log' not in st.session_state:
        st.session_state.position_log = []
    if 'ta_sites' not in st.session_state:
        st.session_state.ta_sites = None
    if 'ta_viability_data' not in st.session_state:
        st.session_state.ta_viability_data = None
    # Remove deprecated fast_mode flag if present
    if 'fast_mode' in st.session_state:
        del st.session_state['fast_mode']
    
    sequence_record = load_genome()
    if not sequence_record:
        return
    
    sequence = sequence_record.seq
    sequence_length = len(sequence)
    
    # Store sequence length in session state for use throughout the app
    st.session_state.sequence_length = sequence_length
    
    # Adjust region end if needed
    if st.session_state.load_region_end > st.session_state.sequence_length:
        st.session_state.load_region_end = st.session_state.sequence_length
    
    annotations = load_optimized_ta_data()

    # Volledige dataset éénmalig laden (genes + TA) op verzoek user
    if annotations is not None and not st.session_state.get('full_genome_loaded', False):
        try:
            # Haal basis stats voor grenzen
            stats = annotations.get_database_stats() if hasattr(annotations, 'get_database_stats') else {'ta_sites_count':0}
            # Load ALL genes
            all_genes = annotations.get_genes_in_region(0, st.session_state.sequence_length)
            # Load ALL TA sites indien nodig
            all_ta = []
            if 'show_tradis_analysis' in st.session_state:
                all_ta = annotations.get_ta_sites_in_region(0, st.session_state.sequence_length) if st.session_state.show_tradis_analysis else []
                # Sla direct raw coverage op; bereken p-waarde buckets
                for site in all_ta:
                    mean = site.get('viability_mean', 0) or 0
                    site['_raw_cov'] = mean
                    std_val = site.get('viability_std', 1.0)
                    mean_val = mean
                    if std_val == 0 and mean_val == 0:
                        site['_p'] = 0.001
                    elif std_val == 0 and mean_val > 0:
                        site['_p'] = 0.01
                    else:
                        site['_p'] = site.get('p_value', 0.5)
            st.session_state.all_genes = all_genes
            st.session_state.all_ta_sites = all_ta
            st.session_state.full_genome_loaded = True
            st.session_state.fast_metrics['full_genome_mode'] = True
        except Exception as e:
            st.warning(f"Kon volledige dataset niet laden: {e}")

    # Lazy load TA sites wanneer toggle van uit->aan gaat en nog niet geladen
    if annotations is not None and st.session_state.get('full_genome_loaded', False):
        if st.session_state.get('show_tradis_analysis', False) and not st.session_state.get('all_ta_sites'):
            try:
                all_ta = annotations.get_ta_sites_in_region(0, st.session_state.sequence_length)
                for site in all_ta:
                    mean = site.get('viability_mean', 0) or 0
                    site['_raw_cov'] = mean
                    std_val = site.get('viability_std', 1.0)
                    mean_val = mean
                    if std_val == 0 and mean_val == 0:
                        site['_p'] = 0.001
                    elif std_val == 0 and mean_val > 0:
                        site['_p'] = 0.01
                    else:
                        site['_p'] = site.get('p_value', 0.5)
                st.session_state.all_ta_sites = all_ta
            except Exception as e:
                st.error(f"TA load error: {e}")
    
    # Gene Search Function
    st.subheader("🔍 Gene Search")
    if annotations is not None:
        col1, col2 = st.columns([3, 1])
        
        with col1:
            # Simple text input for gene search
            search_term = st.text_input(
                "Search for a gene (by name):",
                placeholder="Enter gene name to search...",
                help="Type a gene name to search the database"
            )
        
        with col2:
            col2a, col2b = st.columns(2)
            with col2a:
                if st.button("🎯 Go to Gene") and search_term:
                    # Find the gene using database search
                    gene_match = annotations.search_gene_by_name(search_term) if hasattr(annotations, 'search_gene_by_name') else None
                    
                    if gene_match:
                        gene_center = (gene_match['Begin'] + gene_match['End']) // 2
                    
                    # Check if gene is within current loaded region
                    current_loaded_start = st.session_state.load_region_start
                    current_loaded_end = st.session_state.load_region_end
                    
                    # Always load ±500k around gene center
                    new_region_start = max(0, gene_center - 500_000)
                    new_region_end = min(st.session_state.sequence_length, gene_center + 500_000)
                    st.session_state.load_region_start = new_region_start
                    st.session_state.load_region_end = new_region_end
                    st.session_state.current_start = max(0, gene_center - (st.session_state.current_window // 2))
                    st.session_state.cached_plot = None
                    st.success(f"Loaded region {new_region_start:,}-{new_region_end:,} for {search_term} (center {gene_center:,})")
                    st.rerun()
            
            with col2b:
                if st.button("🔬 Analyze Gene") and search_term:
                    # Set up gene for detailed analysis using SQL database
                    gene_match = annotations.search_gene_by_name(search_term) if hasattr(annotations, 'search_gene_by_name') else None
                    
                    if gene_match:
                        st.session_state.analysis_gene = gene_match
                        st.session_state.show_analysis = True
                        st.success(f"Opening analysis for {search_term}")
                        st.rerun()
                    else:
                        st.error(f"Gene '{search_term}' not found in database")
    
    # Navigatie (verhoogde limiet: grotere vensters toegestaan)
    st.subheader("Navigatie")
    st.session_state.current_window = min(st.session_state.current_window, 1_000_000)
    st.session_state.filter_primary_genes = True
    # Init resolutie voor TA lijn
    if 'ta_line_resolution' not in st.session_state:
        st.session_state.ta_line_resolution = 125  # nieuwe basisresolutie (bp)
    nav_top = st.columns([2,2,2,2,2,2])
    with nav_top[0]:
        # Positie input + validatie
        new_start = st.number_input("Start (bp)", min_value=0, max_value=max(0, st.session_state.sequence_length-100), value=st.session_state.current_start, step=1000)
        if new_start != st.session_state.current_start:
            st.session_state.current_start = new_start
            st.session_state.cached_plot = None
        # Multi-gene search (case-insensitive symbol/name/locus)
        search_term = st.text_input("Gene search (symbol/locus/name)", '')
        matches = []
        if search_term and hasattr(annotations, 'search_genes_multi'):
            try:
                matches = annotations.search_genes_multi(search_term, limit=30)
            except Exception:
                matches = []
        if matches:
            # Build display label
            def fmt(g):
                sym = g.get('symbol') or ''
                loc = g.get('locus_tag') or ''
                prod = g.get('product') or ''
                span = f"{g.get('start_position')}-{g.get('end_position')}"
                primary = sym or loc or g.get('name')
                return f"{primary} | {span} | {prod[:40]}" if prod else f"{primary} | {span}"
            display_options = [fmt(g) for g in matches]
            sel_display = st.selectbox("Matches", display_options, index=0)
            if sel_display:
                sel_idx = display_options.index(sel_display)
                chosen = matches[sel_idx]
                # Center view on chosen gene
                gb = int(chosen.get('start_position',0)); ge=int(chosen.get('end_position',0))
                gene_center = (gb+ge)//2
                st.session_state.current_start = max(0, gene_center - st.session_state.current_window//2)
                st.session_state.selected_gene = (chosen.get('symbol') or chosen.get('locus_tag') or chosen.get('name'))
                st.session_state.cached_plot = None
    # TA resolutie bediening (onder navigatie balk)
    res_col1, res_col2, res_col3 = st.columns([4,1,1])
    with res_col1:
        new_res = st.slider("TA line resolution (bp per step)", min_value=5, max_value=200, step=5, value=int(st.session_state.ta_line_resolution), help="Bepaalt de interpolatie stapgrootte voor de TA lijn over de geladen regio.")
    with res_col2:
        if st.button("-5"):
            st.session_state.ta_line_resolution = max(5, int(st.session_state.ta_line_resolution) - 5)
            new_res = st.session_state.ta_line_resolution
    with res_col3:
        if st.button("+5"):
            st.session_state.ta_line_resolution = min(200, int(st.session_state.ta_line_resolution) + 5)
            new_res = st.session_state.ta_line_resolution
    if new_res != st.session_state.ta_line_resolution:
        st.session_state.ta_line_resolution = new_res
    with nav_top[1]:
        preset_labels = ["5k","10k","20k","30k","100k","300k","500k","700k"]
        preset_map = {"5k":5_000,"10k":10_000,"20k":20_000,"30k":30_000,"100k":100_000,"300k":300_000,"500k":500_000,"700k":700_000}
        current_label = None
        for lbl,val in preset_map.items():
            if abs(st.session_state.current_window - val) < 10:
                current_label = lbl
                break
        default_index = preset_labels.index(current_label) if current_label in preset_labels else 0
        preset = st.selectbox("Window", preset_labels, index=default_index)
        chosen = preset_map[preset]
        if chosen != st.session_state.current_window:
            st.session_state.current_window = chosen
            st.session_state.cached_plot = None
    with nav_top[2]:
        if 'show_tradis_analysis' not in st.session_state:
            st.session_state.show_tradis_analysis = False
        toggle_label = "Hide TA" if st.session_state.show_tradis_analysis else "Show TA"
        if st.button(toggle_label, help="Toon/verberg transposon analyse lijn"):
            st.session_state.show_tradis_analysis = not st.session_state.show_tradis_analysis
            st.session_state.cached_plot = None
        if st.session_state.show_tradis_analysis:
            st.caption("TA line actief (p-value kleur: blauw→geel (0.5)→rood)")
    with nav_top[3]:
        if 'layout_debug' not in st.session_state:
            st.session_state.layout_debug = False
        st.session_state.layout_debug = st.checkbox('Layout debug', value=st.session_state.layout_debug, help='Toon interne y-range en lane informatie')
    with nav_top[4]:
        st.empty()
    with nav_top[5]:
        st.empty()
    # Navigatie knoppen verwijderd op verzoek; alleen directe start & window inputs blijven.

    # Helper functies (inline) - moeten NA gebruik gedefinieerd zijn om lint te vermijden
    # Jump helper
    # (We definieren ze bovenaan voor clarity maar Python voert ze hier ook goed uit.)

    # --- Gene / locus jump invoer ---
    with st.expander("Zoek / Spring naar gen of locus", expanded=False):
        q = st.text_input("Gen-ID of locus (bijv. 123456 of 123000-127000 of genenaam)", key="jump_query")
        colj = st.columns([1,2])
        with colj[0]:
            if st.button("Spring") and q.strip():
                _handle_jump_query(q.strip())
        with colj[1]:
            st.caption("Format: start, start-end, of (prefix van) genenaam. Case-insensitive.")

    # Schaalbalk onder de plot wordt na render toegevoegd; we bewaren keuze dynamisch afhankelijk van venstergrootte.
    
    # Transposon Analysis Status (removed duplicate toggle - use main TradDIS checkbox above)
    st.subheader("🧪 Transposon Analysis (TA Sites)")
    # Extra performance / visual toggles (persist in session)
    perf_col1, perf_col2, perf_col3 = st.columns(3)
    if 'enable_ta_gradient' not in st.session_state:
        st.session_state.enable_ta_gradient = False
    if 'ta_max_points' not in st.session_state:
        st.session_state.ta_max_points = 8000
    with perf_col1:
        st.session_state.enable_ta_gradient = st.checkbox(
            "Gradient lijnen", value=st.session_state.enable_ta_gradient,
            help="Experimenteel: kleurverloop tussen TA punten. Houdt performance in de gaten." )
    with perf_col2:
        st.session_state.ta_max_points = st.number_input(
            "Max points", min_value=1000, max_value=50000, step=1000,
            value=st.session_state.ta_max_points,
            help="Downsampling limiet voor TA sites / performance" )
    with perf_col3:
        if 'last_ta_perf' in st.session_state:
            pinfo = st.session_state.last_ta_perf
            extra = ''
            if pinfo.get('perf_mode'):
                extra = f" | 🚀 Performance mode ({pinfo.get('perf_reason','')})"
            st.caption(f"Rendered {pinfo['rendered_points']}/{pinfo['original_points']} (downsample: {pinfo['downsampled']}){extra} | {pinfo.get('processing_ms',0)} ms | gradient={pinfo.get('gradient_enabled')}" )
   
    # TA analysis data is now loaded directly from SQL database - no pandas needed
    ta_viability_data = None
    ta_statistics_data = None
    ortholog_data = None
    
    # TA analysis parameters with default values
    std_dev_threshold = 0.5
    viability_scale = 1.0

    # Main genome browser with smart caching
    center_position = st.session_state.current_start + (st.session_state.current_window // 2)
    # Log viewer positie
    log = st.session_state.get('position_log', [])
    log.append(int(center_position))
    if len(log) > 500:
        log = log[-500:]
    st.session_state.position_log = log

    # -------- Progressive Loader dicht bij plot --------
    load_row = st.columns([1,4,2])
    with load_row[0]:
        at_end = st.session_state.load_region_end >= sequence_length
        btn_label = "Load next 750k" if not at_end else "End of genome"
        if st.button(btn_label, key="btn_load_next_region_close", disabled=at_end, help="Laad volgende 750k en behoud 250k overlap (sequentieel)"):
            chunk = st.session_state.region_chunk_size
            overlap = st.session_state.region_overlap
            new_start = max(0, st.session_state.load_region_end - overlap)
            new_end = min(sequence_length, new_start + chunk)
            st.session_state.load_region_start = new_start
            st.session_state.load_region_end = new_end
            if not (new_start <= st.session_state.current_start <= new_end):
                st.session_state.current_start = new_start
            st.session_state.cached_plot = None
            st.rerun()
    with load_row[1]:
        span = st.session_state.load_region_end - st.session_state.load_region_start
        st.caption(f"Region {st.session_state.load_region_start:,}-{st.session_state.load_region_end:,} (span {span:,} bp, chunk {st.session_state.region_chunk_size:,}, overlap {st.session_state.region_overlap:,})")
    with load_row[2]:
        st.caption(f"Viewer center: {center_position:,} bp")
    
    # Create cache key for current plot parameters
    current_params = {
        'center': center_position,
        'window': st.session_state.current_window,
        'codons': st.session_state.show_codons,
        'translations': st.session_state.show_translations,
        'frame': st.session_state.reading_frame,
        'region_start': st.session_state.load_region_start,
        'region_end': st.session_state.load_region_end,
        'filter_primary': st.session_state.filter_primary_genes,
        'show_ta': st.session_state.show_tradis_analysis,
        'ta_data_size': len(ta_viability_data) if ta_viability_data else 0,
        'ta_stats_size': len(ta_statistics_data) if ta_statistics_data else 0,
        'ortholog_size': len(ortholog_data) if ortholog_data else 0,
        'std_threshold': std_dev_threshold if st.session_state.show_tradis_analysis else 0,
        'viability_scale': viability_scale if st.session_state.show_tradis_analysis else 1
    }
    
    # Always use unified fast renderer (adaptive tiles + WebGL)
    if (st.session_state.last_plot_params != current_params or 
        st.session_state.cached_plot is None):
        with st.spinner("⚡ Rendering genome view..."):
            try:
                fig = create_fast_genome_browser(
                    sequence,
                    annotations,
                    center_pos=center_position,
                    window_size=st.session_state.current_window
                )
                fm = st.session_state.fast_metrics
                fm['genes_rendered'] = fm.get('genes_rendered', 0)
                fm['ta_rendered'] = fm.get('ta_rendered', 0)
                st.session_state.cached_plot = fig
                st.session_state.last_plot_params = current_params.copy()
            except Exception as e:
                st.error(f"Renderer error: {e}")
                raise
    else:
        fig = st.session_state.cached_plot
    
    # Configure plot for smooth client-side zoom/pan without st.rerun() 
    config = {
        'scrollZoom': True,  # Enable scroll wheel for zoom
        'displayModeBar': True,
        'displaylogo': False,
        'doubleClick': 'reset+autosize',
        'showTips': False,
        'modeBarButtonsToRemove': ['select2d', 'lasso2d', 'autoScale2d'],
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'genome_browser',
            'height': 500,
            'width': 1200,
            'scale': 1
        },
        # Disable interaction throttling for smooth performance
        'responsive': True,
        'editable': False
    }
    
    # Clear navigation instructions  
    st.info("🖱️ **Navigatie:** Muiswiel = Zoom | Sleep = Pan | Dubbelklik = Reset | Shift+Scroll = Horizontaal schuiven | ←/→ = Pan 10%")

    # Display plot with client-side interactions enabled
    plot_container = st.plotly_chart(fig, use_container_width=True, config=config, key="main_genome_plot")

    # --- Dynamische schaalbalk onder de plot ---
    win = st.session_state.current_window
    # Kies mooie eenheid: 1k, 2k, 5k, 10k afhankelijk van window
    import math
    targets = [100,200,500,1000,2000,5000,10000]
    target = targets[0]
    for t in targets:
        if win/5 >= t:
            target = t
    label = f"{target/1000:.0f} kb" if target>=1000 else f"{target} bp"
    # Pixelbreedte approximatie: we gebruiken relatieve div met CSS
    bar_html = f'''<div style="width:100%;display:flex;justify-content:center;margin-top:4px;">
        <div style="position:relative;width:60%;max-width:800px;height:24px;">
            <div style="position:absolute;left:10%;right:10%;top:50%;height:2px;background:#888;"></div>
            <div style="position:absolute;left:35%;width:30%;height:6px;top:50%;margin-top:-3px;background:#444;border:1px solid #222;border-radius:2px;"></div>
            <div style="position:absolute;left:50%;top:50%;transform:translate(-50%, -160%);font-size:11px;color:#333;">{label}</div>
        </div></div>'''
    st.markdown(bar_html, unsafe_allow_html=True)

    # Gene detail panel
    st.subheader("Gene detail")
    sel_gene = st.session_state.get('selected_gene')
    if sel_gene:
        region_genes = get_genes_in_region(annotations, st.session_state.load_region_start, st.session_state.load_region_end)
        gmatch = None
        for g in region_genes:
            sym = g.get('Symbol') or g.get('symbol') or ''
            locus = g.get('Locus tag') or g.get('locus_tag') or ''
            if sel_gene in (sym, locus):
                gmatch = g
                break
        if gmatch:
            gb = gmatch.get('Begin', gmatch.get('start_position',0))
            ge = gmatch.get('End', gmatch.get('end_position',0))
            strand = gmatch.get('strand','+')
            seq_slice = sequence[gb:ge]
            from Bio.Seq import Seq as _Seq
            if strand == '-':
                seq_slice = seq_slice.reverse_complement()
            seq_str = str(seq_slice).upper()
            def translate_frames(s):
                frames = []
                for frame in range(3):
                    codons = [s[i:i+3] for i in range(frame, len(s)-2, 3)]
                    aa = ''.join([str(_Seq(c).translate()) if len(c)==3 else '' for c in codons])
                    frames.append((frame, aa))
                return frames
            frames = translate_frames(seq_str)
            st.markdown(f"**{sel_gene}** Pos: {gb:,}-{ge:,} Len: {ge-gb} bp Strand: {strand}")
            with st.expander("Coding sequence", expanded=False):
                st.code(seq_str, language='text')
            with st.expander("Translations (3 frames)", expanded=False):
                for f, aa in frames:
                    st.markdown(f"Frame {f}:\n``{aa}``")
        else:
            st.caption("Geen details beschikbaar (buiten geladen regio)")
    else:
        st.caption("Geen gen geselecteerd")

    # Unified JS block: relayout debounce -> query param sync + ctrl+wheel horizontal pan
    components.html(
        f"""
        <script>
        (function(){{
            function debounce(fn,ms){{let t;return (...a)=>{{clearTimeout(t);t=setTimeout(()=>fn(...a),ms)}}}}
            const plot = window.parent.document.querySelector('.js-plotly-plot');
            if(!plot) return;
            let lastReloadStart = {st.session_state.current_start};
            const reloadDebounced = debounce(function(range){{
                if(!Array.isArray(range)||range.length<2) return;
                const start = Math.round(range[0]);
                const span = Math.max(1, range[1]-range[0]);
                // Reload alleen als we verder dan 40% van het venster verschoven zijn
                if(Math.abs(start - lastReloadStart) > span * 0.40){{
                    lastReloadStart = start;
                    const url = new URL(window.parent.location.href);
                    url.searchParams.set('start', start.toString());
                    window.parent.location.href = url.toString();
                }}
            }}, 250);
            window.addEventListener('plotly_relayout', function(e){{
                const r0=e['xaxis.range[0]']; const r1=e['xaxis.range[1]'];
                if(r0===undefined||r1===undefined) return; reloadDebounced([r0,r1]);
            }});
            function shiftRange(fracDir){{
                if(!plot || !plot._fullLayout) return;
                const xr = plot._fullLayout.xaxis.range.slice();
                const span = xr[1]-xr[0];
                const delta = span * fracDir;
                const nr = [xr[0]+delta, xr[1]+delta];
                Plotly.relayout(plot, {{"xaxis.range": nr}});
            }}
            // Shift + scroll horizontaal verschuiven (client-side)
            window.parent.document.addEventListener('wheel', function(e){{
                if(!e.shiftKey) return; if(!plot || !plot._fullLayout) return; if(!plot.contains(e.target) && !e.target.closest('.js-plotly-plot')) return;
                e.preventDefault();
                const span = plot._fullLayout.xaxis.range[1]-plot._fullLayout.xaxis.range[0];
                const deltaFrac = 0.12 * (e.deltaY>0?1:-1);
                shiftRange(deltaFrac);
            }}, false);
            // Arrow keys voor 10% pan
            window.parent.document.addEventListener('keydown', function(e){{
                if(e.key === 'ArrowLeft'){{ e.preventDefault(); shiftRange(-0.10); }}
                else if(e.key === 'ArrowRight'){{ e.preventDefault(); shiftRange(0.10); }}
            }});
        }})();
        </script>
        """,
        height=0
    )
    
    # Client-side interactions handled by Plotly - no st.rerun() needed for zoom/pan
    
# Native Plotly zoom enabled via scrollZoom: True in config
    
    # Lightweight single-line status instead of large duplicated panel
    st.caption(f"Center {center_position:,} | Range {st.session_state.current_start:,}-{st.session_state.current_start + st.session_state.current_window:,} | Genome {st.session_state.sequence_length:,} bp")
    
    # Gene translations panel
    if st.session_state.show_translations and annotations is not None:
        st.subheader("Gene Translations")
        
        # Gene selection using database
        if hasattr(annotations, 'get_database_stats'):
            # Get a sample of genes from current visible region for selection
            current_genes = get_genes_in_region(annotations, 
                                              st.session_state.current_start, 
                                              st.session_state.current_start + st.session_state.current_window)
            
            if len(current_genes) > 0:
                # Extract unique gene names from list of dicts
                gene_names = list(set([g.get('Name', '') for g in current_genes if g.get('Name', '').strip()]))[:20]
            else:
                gene_names = []
        else:
            gene_names = []
        
        selected_gene = st.selectbox("Select Gene for Translation", 
                                   options=[''] + list(gene_names))
        
        if selected_gene:
            gene_info = annotations.search_gene_by_name(selected_gene) if hasattr(annotations, 'search_gene_by_name') else None
            
            if gene_info:
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**Gene:** {gene_info.get('Product', gene_info.get('Name', 'Unknown'))}")
                    st.write(f"**Position:** {gene_info['Begin']:,} - {gene_info['End']:,}")
                    st.write(f"**Strand:** {gene_info.get('Strand', '+')}")
                    st.write(f"**Length:** {gene_info['End'] - gene_info['Begin']:,} bp")
                
                with col2:
                    if 'Type' in gene_info:
                        st.write(f"**Type:** {gene_info.get('Type', 'N/A')}")
                    if 'Name' in gene_info:
                        st.write(f"**Name:** {gene_info.get('Name', 'N/A')}")
                    if gene_info.get('ortholog_data'):
                        st.write("**Has ortholog data:** ✅")
            else:
                st.error(f"Gene '{selected_gene}' not found in database")
            
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
            # Get gene count from database stats
            stats = annotations.get_database_stats() if hasattr(annotations, 'get_database_stats') else {'genes_count': 0}
            st.write(f"**Totaal Genen:** {stats['genes_count']:,}")
            
            # Gene type distribution (skip since we can't easily get this from SQL without full query)
            # st.write("**Gen Types:**")
            # Could implement if needed by querying database for gene types
        
        st.write("### Navigatie")
        st.write("🔍 **Zoomen:** Muiswiel op plot")
        st.write("⬅️➡️ **Pannen:** Sleep met muis")
        st.write("🎯 **Reset:** Dubbelklik op plot")
        st.write("⏭️ **Navigeer:** Gebruik knoppen boven plot")
        
      
        
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

    
    # Gene Analysis Window (tijdelijk uitgeschakeld - functie verwijderd)
    if st.session_state.get('show_analysis'):
        st.info("Gene analysis module is tijdelijk verwijderd tijdens opschoning.")

if __name__ == "__main__":
    create_genome_browser_page()
