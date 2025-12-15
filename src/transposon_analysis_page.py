"""
Transposon Analysis Page
------------------------
Dedicated page for TradDIS transposon insertion analysis with statistical visualization.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pickle
from pathlib import Path
import json

def load_optimized_data():
    """Load ALL TA sites from high-performance database - uncoupled from region."""
    try:
        from .database_v2 import HighPerformanceGenomicDB
        
        # Initialize high-performance database
        if not hasattr(st.session_state, 'transposon_hp_db'):
            with st.spinner('🗄️ Loading high-performance genomic database...'):
                st.session_state.transposon_hp_db = HighPerformanceGenomicDB()
        
        return st.session_state.transposon_hp_db
    except Exception as e:
        st.error(f"Error loading database: {e}")
        return None

def get_ta_sites_in_region(database, start_pos, end_pos):
    """Get TA sites within genomic region using SQL database."""
    if not database:
        return []
    
    try:
        return database.get_ta_sites_in_region(start_pos, end_pos)
    except Exception as e:
        st.error(f"Error querying TA sites: {e}")
        return []

def create_smooth_tradis_line_analysis(ta_sites, start_pos, end_pos, smoothing_window=1000):
    """Create smooth continuous TradDIS line with color interpolation for analysis page."""
    if not ta_sites:
        return [], [], []
    
    # Create position array for interpolation
    positions = np.arange(start_pos, end_pos + 1, step=100)  # Every 100bp
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
                values.append(site['viability_mean'])
                p_vals.append(site['p_value'])
            
            # Weighted average
            total_weight = sum(weights)
            if total_weight > 0:
                heights[i] = sum(v * w for v, w in zip(values, weights)) / total_weight
                p_values[i] = sum(p * w for p, w in zip(p_vals, weights)) / total_weight
    
    return positions.tolist(), heights.tolist(), p_values.tolist()

def get_color_from_pvalue(p_value):
    """Convert p-value to color for smooth transitions with interpolation."""
    # Define color points for interpolation
    # Blue -> Orange -> Red -> Dark Red
    if p_value >= 0.05:
        # Not significant - Blue
        return 'rgb(0, 100, 255)'
    elif p_value >= 0.01:
        # Interpolate between blue and orange
        # p_value range: 0.01 to 0.05, map to 0-1
        t = (0.05 - p_value) / (0.05 - 0.01)
        # Blue (0, 100, 255) -> Orange (255, 165, 0)
        r = int(0 + t * 255)
        g = int(100 + t * 65)
        b = int(255 - t * 255)
        return f'rgb({r}, {g}, {b})'
    elif p_value >= 0.001:
        # Interpolate between orange and red
        t = (0.01 - p_value) / (0.01 - 0.001)
        # Orange (255, 165, 0) -> Red (255, 0, 0)
        r = 255
        g = int(165 - t * 165)
        b = 0
        return f'rgb({r}, {g}, {b})'
    else:
        # Interpolate between red and dark red
        t = min(1.0, (0.001 - p_value) / 0.0005)  # Cap at very small p-values
        # Red (255, 0, 0) -> Dark Red (139, 0, 0)
        r = int(255 - t * 116)
        g = 0
        b = 0
        return f'rgb({r}, {g}, {b})'

def calculate_line_segments(ta_sites, window_size=1000):
    """Calculate continuous line segments for visualization."""
    if not ta_sites:
        return [], [], []
    
    # Sort sites by position
    sorted_sites = sorted(ta_sites, key=lambda x: x['position'])
    
    # Group sites into segments
    segments = []
    current_segment = []
    
    for i, site in enumerate(sorted_sites):
        if not current_segment:
            current_segment.append(site)
        else:
            # Check if site is within window of previous site
            if site['position'] - current_segment[-1]['position'] <= window_size:
                current_segment.append(site)
            else:
                # Start new segment
                if len(current_segment) > 1:
                    segments.append(current_segment)
                current_segment = [site]
    
    # Add final segment
    if len(current_segment) > 1:
        segments.append(current_segment)
    
    return segments

def create_integrated_ta_gene_plot(ta_sites, genes, p_value_threshold=0.05, height_scale=1.0, start_pos=None, end_pos=None):
    """Create integrated TA site + gene visualization like genome browser."""
    
    if not ta_sites:
        return go.Figure()
    
    # Determine plot range
    if start_pos is None:
        start_pos = min(site['position'] for site in ta_sites)
    if end_pos is None:
        end_pos = max(site['position'] for site in ta_sites)
    
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
    
    # TA sites smooth line (above baseline)
    smooth_x, smooth_y, smooth_p = create_smooth_tradis_line_analysis(ta_sites, start_pos, end_pos, smoothing_window=2000)
    
    if smooth_x and smooth_y:
        # Scale heights above baseline
        scaled_y = [0.5 + (y * height_scale) for y in smooth_y]
        
        fig.add_trace(go.Scatter(
            x=smooth_x,
            y=scaled_y,
            mode='lines',
            line=dict(color='blue', width=3),
            name='TA Coverage',
            hovertemplate='Position: %{x:,} bp<br>Coverage: %{customdata:.2f}<br><extra></extra>',
            customdata=smooth_y
        ))
    
    # Add individual TA sites as markers  
    ta_positions = [site['position'] for site in ta_sites]
    ta_means = [site.get('viability_mean', 0) for site in ta_sites]
    ta_p_values = [site.get('p_value', 1) for site in ta_sites]
    
    # Color by significance
    colors = ['red' if p < p_value_threshold else 'blue' for p in ta_p_values]
    
    fig.add_trace(go.Scatter(
        x=ta_positions,
        y=[0.1] * len(ta_positions),  # Just above baseline
        mode='markers',
        marker=dict(
            size=6,
            color=colors,
            opacity=0.7
        ),
        name='TA Sites',
        hovertemplate='<b>TA Site</b><br>Position: %{x:,} bp<br>Coverage: %{customdata[0]:.2f}<br>P-value: %{customdata[1]:.3e}<br><extra></extra>',
        customdata=list(zip(ta_means, ta_p_values))
    ))
    
    # Add genes (below baseline like genome browser)
    if genes:
        from .genome_browser import calculate_gene_levels, create_arrow_shape
        
        # Convert genes to proper format
        formatted_genes = []
        for gene in genes:
            formatted_genes.append({
                'Begin': gene['start_position'],
                'End': gene['end_position'],
                'strand': gene['strand'],
                'Name': gene.get('name', gene.get('locus_tag', f"gene_{gene['start_position']}"))
            })
        
        # Calculate gene levels
        gene_levels = calculate_gene_levels(formatted_genes, ta_analysis_active=False)
        
        for i, gene in enumerate(genes):
            gene_start = gene['start_position']
            gene_end = gene['end_position']
            strand = gene['strand']
            product = gene.get('product', 'Unknown')
            locus_tag = gene.get('locus_tag', '')
            
            gene_name = gene.get('name', locus_tag)
            y_pos = gene_levels.get(gene_name, -0.5 - (i * 0.4))
            
            color = "blue" if strand == "+" else "red"
            
            # Create gene arrow
            arrow_x, arrow_y = create_arrow_shape(gene_start, gene_end, y_pos, strand, height=0.3)
            
            fig.add_trace(go.Scatter(
                x=arrow_x,
                y=arrow_y,
                fill="toself",
                fillcolor=color,
                line=dict(color=color, width=1),
                mode="lines",
                name=product,
                showlegend=False,
                hovertemplate=f'<b>{product}</b><br>Locus: {locus_tag}<br>Position: {gene_start:,} - {gene_end:,}<br>Strand: {strand}<br><extra></extra>'
            ))
    
    # Layout
    fig.update_layout(
        title=f"🧬 TradDIS + Gene Analysis - Position {start_pos:,} to {end_pos:,}",
        xaxis=dict(
            title="Genome Position (bp)",
            range=[start_pos, end_pos],
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'
        ),
        yaxis=dict(
            title="",
            range=[-3.0, 4.0],
            showticklabels=False,
            showgrid=False,
            zeroline=True,
            zerolinecolor='black',
            zerolinewidth=2
        ),
        height=600,
        plot_bgcolor="white",
        showlegend=True
    )
    
    return fig

def show_gene_ta_correlations(genes, ta_sites):
    """Show gene-TA site correlation analysis."""
    if not genes or not ta_sites:
        st.info("No genes or TA sites available for correlation analysis.")
        return
    
    correlations = []
    
    for gene in genes:
        gene_start = gene['start_position']
        gene_end = gene['end_position']
        gene_name = gene.get('name', gene.get('locus_tag', f"gene_{gene_start}"))
        
        # Find TA sites within gene
        gene_ta_sites = [site for site in ta_sites if gene_start <= site['position'] <= gene_end]
        
        if gene_ta_sites:
            avg_coverage = np.mean([site['viability_mean'] for site in gene_ta_sites])
            avg_p_value = np.mean([site['p_value'] for site in gene_ta_sites])
            significant_sites = len([site for site in gene_ta_sites if site['p_value'] < 0.05])
            
            correlations.append({
                'Gene': gene_name,
                'Product': gene.get('product', 'Unknown'),
                'Strand': gene['strand'],
                'TA Sites': len(gene_ta_sites),
                'Avg Coverage': avg_coverage,
                'Avg P-value': avg_p_value,
                'Significant Sites': significant_sites,
                'Position': f"{gene_start:,} - {gene_end:,}"
            })
    
    if correlations:
        st.dataframe(correlations, use_container_width=True)
    else:
        st.info("No genes contain TA sites in the current region.")

def create_continuous_line_plot(ta_sites, p_value_threshold=0.05, height_scale=1.0, start_pos=None, end_pos=None):
    """Create smooth continuous line plot with p-value colors and mean-based heights."""
    
    fig = go.Figure()
    
    if not ta_sites:
        st.warning("No TA sites data available")
        return fig
    
    # Determine range
    if start_pos is None:
        start_pos = min(site['position'] for site in ta_sites)
    if end_pos is None:
        end_pos = max(site['position'] for site in ta_sites)
    
    # Create smooth continuous line
    positions, heights, p_values = create_smooth_tradis_line_analysis(
        ta_sites, start_pos, end_pos, smoothing_window=2000
    )
    
    if not positions or not heights:
        st.warning("No valid data for continuous line")
        return fig
    
    # Scale heights
    max_height = max(heights) if heights else 1
    if max_height > 0:
        scaled_heights = [(h / max_height) * height_scale for h in heights]
    else:
        scaled_heights = [0] * len(positions)
    
    # Create smooth color-changing line using markers with color mapping
    # Convert p-values to numeric color scale (0=blue, 1=dark red)
    color_values = []
    for p in p_values:
        if p >= 0.05:
            color_values.append(0.0)  # Blue
        elif p >= 0.01:
            color_values.append(0.33)  # Orange
        elif p >= 0.001:
            color_values.append(0.66)  # Red
        else:
            color_values.append(1.0)  # Dark red
    
    # Add the main continuous line with smooth color gradients
    fig.add_trace(go.Scatter(
        x=positions,
        y=scaled_heights,
        mode='lines+markers',
        line=dict(width=3, color='rgba(0,0,0,0)'),  # Invisible line
        marker=dict(
            size=4,
            color=color_values,
            colorscale=[
                [0.0, 'rgb(0, 100, 255)'],      # Blue (not significant)
                [0.33, 'rgb(255, 165, 0)'],    # Orange (significant)
                [0.66, 'rgb(255, 0, 0)'],      # Red (very significant)
                [1.0, 'rgb(139, 0, 0)']        # Dark red (highly significant)
            ],
            showscale=True,
            colorbar=dict(
                title="Significance<br>Level",
                tickmode="array",
                tickvals=[0.0, 0.33, 0.66, 1.0],
                ticktext=["p≥0.05", "p<0.05", "p<0.01", "p<0.001"],
                x=1.02
            ),
            line=dict(width=1, color='white')  # Marker outlines
        ),
        name='TradDIS Analysis',
        hovertemplate='<b>TradDIS Coverage</b><br>' +
                    'Position: %{x:,} bp<br>' +
                    'Coverage: %{customdata[0]:.2f}<br>' +
                    'P-value: %{customdata[1]:.2e}<br>' +
                    '<extra></extra>',
        customdata=list(zip(heights, p_values)),
        showlegend=True
    ))
    
    # Add connecting lines between markers for true continuous appearance
    for i in range(0, len(positions) - 1, 5):  # Every 5th point for performance
        # Get colors for gradient
        start_color = get_color_from_pvalue(p_values[i])
        end_color = get_color_from_pvalue(p_values[min(i+5, len(p_values)-1)])
        
        fig.add_trace(go.Scatter(
            x=[positions[i], positions[min(i+5, len(positions)-1)]],
            y=[scaled_heights[i], scaled_heights[min(i+5, len(scaled_heights)-1)]],
            mode='lines',
            line=dict(
                width=4,
                color=start_color  # Use start color for segment
            ),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    # Update layout
    fig.update_layout(
        title=f"TradDIS Smooth Continuous Line - Genomic Position {start_pos:,} to {end_pos:,}",
        xaxis_title="Genomic Position (bp)",
        yaxis_title=f"Scaled Mean Coverage (max scale: {height_scale}x)",
        height=600,
        hovermode='closest',
        showlegend=False
    )
    
    return fig

def create_statistical_summary(ta_sites):
    """Create statistical summary plots."""
    if not ta_sites:
        return None
    
    # Extract data
    positions = [site['position'] for site in ta_sites]
    means = [site['viability_mean'] for site in ta_sites]
    p_values = [site['p_value'] for site in ta_sites]
    stds = [site['viability_std'] for site in ta_sites]
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=['P-value Distribution', 'Mean Coverage Distribution', 
                       'Standard Deviation Distribution', 'Significance vs Coverage'],
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # P-value histogram
    fig.add_trace(
        go.Histogram(x=p_values, nbinsx=50, name="P-values", 
                    marker_color='lightblue'),
        row=1, col=1
    )
    
    # Mean coverage histogram  
    fig.add_trace(
        go.Histogram(x=means, nbinsx=50, name="Mean Coverage",
                    marker_color='lightgreen'),
        row=1, col=2
    )
    
    # Standard deviation histogram
    fig.add_trace(
        go.Histogram(x=stds, nbinsx=50, name="Std Deviation",
                    marker_color='lightyellow'),
        row=2, col=1
    )
    
    # Scatter plot: significance vs coverage
    colors = ['red' if p < 0.05 else 'blue' for p in p_values]
    fig.add_trace(
        go.Scatter(x=means, y=[-np.log10(p) for p in p_values],
                  mode='markers', name="Sites",
                  marker=dict(color=colors, size=4, opacity=0.6)),
        row=2, col=2
    )
    
    # Update layout
    fig.update_layout(height=800, showlegend=False)
    fig.update_xaxes(title_text="P-value", row=1, col=1)
    fig.update_xaxes(title_text="Mean Coverage", row=1, col=2)
    fig.update_xaxes(title_text="Standard Deviation", row=2, col=1)
    fig.update_xaxes(title_text="Mean Coverage", row=2, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=2, col=1)
    fig.update_yaxes(title_text="-log10(P-value)", row=2, col=2)
    
    return fig

def show():
    """Main function to display the transposon analysis page."""
    st.title("🧬 TradDIS Transposon Analysis")
    st.markdown("---")
    
    # Load data
    with st.spinner("Loading database..."):
        database = load_optimized_data()
    
    if database is None:
        st.error("Failed to load database. Please ensure the database has been initialized.")
        return
    
    # Load ALL TA sites - uncoupled from position, then filter for display
    with st.spinner("Loading all TA sites from database..."):
        ta_sites = database.get_ta_sites_in_region(0, 10000000)  # Get everything
    
    if not ta_sites:
        st.warning("No TA sites data found in database.")
        return
    
    # DEFAULT: Display first 1M bases like genome browser
    DEFAULT_DISPLAY_END = 1000000
    display_ta_sites = [site for site in ta_sites if site['position'] <= DEFAULT_DISPLAY_END]
    
    st.info(f"📍 **Loaded {len(ta_sites):,} total TA sites, displaying first 1M bases ({len(display_ta_sites):,} sites)**")
    
    # Display summary statistics for both ALL and DISPLAYED data  
    st.subheader("📊 Data Summary")
    
    # Statistics for displayed 1M region
    display_significant = len([site for site in display_ta_sites if site['p_value'] < 0.05])
    display_high_coverage = len([site for site in display_ta_sites if site['viability_mean'] > 10])
    
    # Statistics for complete dataset
    total_sites = len(ta_sites)
    total_significant = len([site for site in ta_sites if site['p_value'] < 0.05])
    total_high_coverage = len([site for site in ta_sites if site['viability_mean'] > 10])
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Displayed Sites (1M)", f"{len(display_ta_sites):,}", 
                 f"Total: {total_sites:,}")
    
    with col2:
        st.metric("Significant (Display)", f"{display_significant:,}", 
                 f"Total: {total_significant:,}")
    
    with col3:
        st.metric("High Coverage (Display)", f"{display_high_coverage:,}",
                 f"Total: {total_high_coverage:,}")
    
    with col4:
        display_mean_p = np.mean([site['p_value'] for site in display_ta_sites]) if display_ta_sites else 0
        st.metric("Mean P-value (Display)", f"{display_mean_p:.2e}")
    
    # Analysis controls
    st.subheader("🔧 Analysis Controls")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        p_value_threshold = st.slider(
            "P-value Threshold", 
            min_value=0.001, 
            max_value=0.1, 
            value=0.05, 
            step=0.001,
            format="%.3f"
        )
    
    with col2:
        height_scale = st.slider(
            "Height Scale Factor",
            min_value=0.1,
            max_value=5.0,
            value=1.0,
            step=0.1
        )
    
    with col3:
        region_filter = st.selectbox(
            "Region Filter",
            options=["All", "Significant Only", "High Coverage Only"],
            index=0
        )
    
    with col4:
        plot_type = st.selectbox(
            "Visualization Type",
            options=["Smooth Continuous Line", "Segmented Regions"],
            index=0
        )
    
    # Region selection for visualization - LINKED TO GENOME BROWSER
    if plot_type == "Smooth Continuous Line":
        st.subheader("📍 Genomic Region Visualization")
        
        # Get genome range from loaded data
        all_positions = [site['position'] for site in ta_sites]
        min_pos = min(all_positions) if all_positions else 0
        max_pos = max(all_positions) if all_positions else DEFAULT_DISPLAY_END
        
        # Auto-link with genome browser OR default to 1M region
        if hasattr(st.session_state, 'current_start') and hasattr(st.session_state, 'current_window'):
            # Use genome browser's current view
            default_start = st.session_state.current_start
            default_end = st.session_state.current_start + st.session_state.current_window
            st.info(f"🔗 **Linked to Gene Browser View:** {default_start:,} - {default_end:,} bp")
        else:
            # Default to first 1M bases like genome browser
            default_start = 0
            default_end = DEFAULT_DISPLAY_END
            st.info(f"🧬 **Default 1M Base View:** {default_start:,} - {default_end:,} bp")
        
        col1, col2 = st.columns(2)
        with col1:
            start_pos = st.number_input(
                "Start Position (bp)",
                min_value=0,  # Always allow from genome start
                max_value=max_pos,
                value=max(0, default_start),  # Ensure valid value
                step=1000
            )
        
        with col2:
            end_pos = st.number_input(
                "End Position (bp)", 
                min_value=1000,  # Minimum meaningful window
                max_value=max_pos,
                value=max(1000, default_end),  # Ensure valid value
                step=1000
            )
        
        if start_pos >= end_pos:
            st.error("Start position must be less than end position!")
            return
    else:
        start_pos = None
        end_pos = None

    # Filter data based on controls
    filtered_sites = ta_sites.copy()
    
    if region_filter == "Significant Only":
        filtered_sites = [site for site in filtered_sites if site['p_value'] < p_value_threshold]
    elif region_filter == "High Coverage Only":
        filtered_sites = [site for site in filtered_sites if site['viability_mean'] > 10]
    
    # Filter by region for visualization (from pre-loaded data - no additional DB queries)
    if plot_type == "Smooth Continuous Line" and start_pos is not None and end_pos is not None:
        # Filter from already loaded data - much faster than DB query
        region_sites = [site for site in ta_sites if start_pos <= site['position'] <= end_pos]
        
        # Apply filtering to region sites
        if region_filter == "Significant Only":
            filtered_sites = [site for site in region_sites if site['p_value'] < p_value_threshold]
        elif region_filter == "High Coverage Only":
            filtered_sites = [site for site in region_sites if site['viability_mean'] > 10]
        else:
            filtered_sites = region_sites
        
        # Load corresponding genes for the region to show gene-TA site links
        region_genes = database.get_genes_in_region(start_pos, end_pos)
        st.info(f"🧬 **Region contains {len(region_genes)} genes and {len(filtered_sites)} TA sites**")
    
    # Main visualization with gene overlay like genome browser
    if plot_type == "Smooth Continuous Line":
        st.subheader("📈 TradDIS Analysis with Gene Context")
        st.markdown(f"**Region {start_pos:,} - {end_pos:,} bp: {len(filtered_sites):,} TA sites, {len(region_genes) if 'region_genes' in locals() else 0} genes**")
        
        if filtered_sites:
            with st.spinner("Creating integrated TA/gene visualization..."):
                # Create enhanced plot with gene overlay like genome browser
                line_plot = create_integrated_ta_gene_plot(
                    filtered_sites,
                    region_genes if 'region_genes' in locals() else [],
                    p_value_threshold=p_value_threshold,
                    height_scale=height_scale,
                    start_pos=start_pos,
                    end_pos=end_pos
                )
            
            st.plotly_chart(line_plot, use_container_width=True)
            
            # Show gene-TA site correlation info
            if 'region_genes' in locals() and region_genes:
                st.subheader("🔗 Gene-TA Site Correlations")
                show_gene_ta_correlations(region_genes, filtered_sites)
                
        else:
            st.warning("No TA sites found in the selected region.")
    
    else:  # Segmented Regions
        st.subheader("📈 Segmented Line Visualization")
        st.markdown(f"**Showing {len(filtered_sites):,} of {total_sites:,} TA sites**")
        
        if filtered_sites:
            # Calculate line segments
            segments = calculate_line_segments(filtered_sites, window_size=5000)
            
            # Use original segmented plotting
            fig = go.Figure()
            
            if segments:
                all_means = [site['viability_mean'] for site in filtered_sites]
                max_mean = max(all_means) if all_means else 1
                
                for seg_idx, segment in enumerate(segments):
                    if len(segment) < 2:
                        continue
                        
                    positions = [site['position'] for site in segment]
                    means = [site['viability_mean'] for site in segment]
                    p_values = [site['p_value'] for site in segment]
                    
                    heights = [(mean / max_mean) * height_scale for mean in means]
                    
                    colors = []
                    for p_val in p_values:
                        if p_val < 0.001:
                            colors.append('darkred')
                        elif p_val < 0.01:
                            colors.append('red')
                        elif p_val < p_value_threshold:
                            colors.append('orange')
                        else:
                            colors.append('blue')
                    
                    fig.add_trace(go.Scatter(
                        x=positions,
                        y=heights,
                        mode='lines+markers',
                        line=dict(width=3, color=colors[0] if colors else 'blue'),
                        marker=dict(size=6, color=colors),
                        text=[f"Position: {pos:,}<br>Mean: {mean:.2f}<br>P-value: {p_val:.2e}" 
                              for pos, mean, p_val in zip(positions, means, p_values)],
                        hovertemplate="%{text}<extra></extra>",
                        name=f"Segment {seg_idx + 1}",
                        showlegend=False
                    ))
                
                fig.update_layout(
                    title="TradDIS Segmented Analysis",
                    xaxis_title="Genomic Position (bp)",
                    yaxis_title=f"Scaled Mean Coverage (max scale: {height_scale}x)",
                    height=600,
                    hovermode='closest'
                )
            
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("No sites match the current filter criteria.")
    
    # Color legend explanation (for both plot types)
    st.markdown("""
    **Color Legend:**
    - 🔵 **Blue**: Not significant (p ≥ 0.05)
    - 🟠 **Orange**: Significant (p < 0.05)
    - 🔴 **Red**: Very significant (p < 0.01)
    - 🔴 **Dark Red**: Highly significant (p < 0.001)
    """)
    
    # Statistical summary plots (for filtered data)
    if filtered_sites:
        st.subheader("📊 Statistical Summary")
        
        with st.spinner("Creating statistical summary..."):
            stats_plot = create_statistical_summary(filtered_sites)
        
        if stats_plot:
            st.plotly_chart(stats_plot, use_container_width=True)
        
        # Data export options
        st.subheader("💾 Data Export")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("📋 Copy Summary Stats"):
                summary_stats = {
                    "total_sites": total_sites,
                    "display_sites": len(display_ta_sites),
                    "total_significant": total_significant,
                    "display_significant": display_significant,
                    "total_high_coverage": total_high_coverage,
                    "display_high_coverage": display_high_coverage,
                    "display_mean_p_value": display_mean_p,
                    "filtered_sites": len(filtered_sites)
                }
                st.code(json.dumps(summary_stats, indent=2))
        
        with col2:
            # Create downloadable CSV
            if st.button("📥 Prepare CSV Download"):
                df_export = pd.DataFrame([
                    {
                        'position': site['position'],
                        'mean_coverage': site['viability_mean'],
                        'std_deviation': site['std'],
                        'p_value': site['p_value'],
                        'significant': site['p_value'] < p_value_threshold
                    }
                    for site in filtered_sites
                ])
                
                csv_data = df_export.to_csv(index=False)
                st.download_button(
                    label="Download Filtered Data",
                    data=csv_data,
                    file_name=f"tradis_analysis_filtered_{len(filtered_sites)}_sites.csv",
                    mime="text/csv"
                )
    
    # Additional information
    st.markdown("---")
    st.markdown("""
    ### About This Analysis
    
    This page provides dedicated analysis of TradDIS (Transposon Directed Insertion-site Sequencing) data:
    
    #### Visualization Types:
    - **Smooth Continuous Line**: True continuous visualization with interpolated values and smooth color transitions
    - **Segmented Regions**: Traditional segment-based visualization with discrete regions
    
    #### Features:
    - **Height**: Based on mean coverage values, scaled by the height factor
    - **Color**: Statistical significance (p-values from t-tests) with smooth transitions
    - **Interactive**: Hover over points for detailed information
    - **Regional Focus**: Select specific genomic regions for detailed analysis
    
    The data has been pre-processed and optimized for fast loading. All statistical analysis 
    (96,259 TA sites, t-tests, p-values) is computed once and cached for performance.
    """)

if __name__ == "__main__":
    show()
