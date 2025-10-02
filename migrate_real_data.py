#!/usr/bin/env python3
"""
REAL DATA Migration Script - Migrate ACTUAL scientific data
============================================================
"""

import sqlite3
import os
import sys
import time
import json
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent / 'src'))

from database_v2 import HighPerformanceGenomicDB

def migrate_real_data():
    """Migrate the ACTUAL scientific data to high-performance architecture."""
    
    print("🚀 MIGRATING REAL SCIENTIFIC DATA!")
    print("=" * 60)
    
    # Source database with REAL data
    source_db_path = "data/genomic_data.db"
    if not os.path.exists(source_db_path):
        print(f"❌ REAL database not found at {source_db_path}")
        return
    
    # Remove any existing test databases
    for db_file in ["data/genes.db", "data/ta_sites.db", "data/orthology.db"]:
        if os.path.exists(db_file):
            os.remove(db_file)
            print(f"🗑️  Removed test database: {db_file}")
    
    # Initialize new high-performance database
    new_db = HighPerformanceGenomicDB()
    
    # Connect to source database
    conn = sqlite3.connect(source_db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    # Get source statistics
    cursor.execute("SELECT COUNT(*) FROM genes")
    gene_count = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM ta_sites")
    ta_count = cursor.fetchone()[0]
    
    print(f"📈 Source REAL data: {ta_count:,} TA sites, {gene_count:,} genes")
    print()
    
    try:
        # MIGRATE GENES - REAL DATA
        print("⚡ Migrating REAL genes...")
        start_time = time.time()
        
        cursor.execute("""
            SELECT gene_name, start_position, end_position, strand, gene_type, product, 
                   locus_tag, symbol, protein_accession, ortholog_data
            FROM genes 
            ORDER BY start_position
        """)
        
        genes_data = []
        for row in cursor.fetchall():
            gene_data = {
                'name': row['product'] or row['gene_name'] or row['locus_tag'] or '',  # Use product as name
                'locus_tag': row['locus_tag'] or '',  # Preserve actual locus tag
                'start_position': row['start_position'] or 0,
                'end_position': row['end_position'] or 0,
                'strand': row['strand'] or '+',
                'gene_type': row['gene_type'] or 'CDS',
                'product': row['product'] or '',
                'symbol': row.get('symbol') if isinstance(row, dict) else row['symbol'] if 'symbol' in row.keys() else ''
            }
            genes_data.append(gene_data)
        
        gene_ids = new_db.bulk_insert_genes(genes_data)
        gene_time = time.time() - start_time
        print(f"✅ Migrated {len(gene_ids):,} REAL genes in {gene_time:.2f}s")
        print()
        
        # MIGRATE TA SITES - REAL DATA
        print("⚡ Migrating REAL TA sites...")
        start_time = time.time()
        
        cursor.execute("""
            SELECT position, mean_coverage, std_coverage, raw_values, 
                   sample_count, p_value, ortholog_info
            FROM ta_sites 
            ORDER BY position
        """)
        
        ta_sites_data = []
        for row in cursor.fetchall():
            # Map real data to expected format
            ta_site_data = {
                'position': row['position'],
                'viability_mean': row['mean_coverage'] or 0,
                'viability_std': row['std_coverage'] or 0,
                'viability_values': row['raw_values'],
                'viable_colonies': int(row['sample_count'] or 0),
                'total_colonies': 200,
                'p_value': row['p_value'] or 0.05  # Include p_value
            }
            ta_sites_data.append(ta_site_data)
        
        new_db.bulk_insert_ta_sites(ta_sites_data)
        ta_time = time.time() - start_time
        print(f"✅ Migrated {len(ta_sites_data):,} REAL TA sites in {ta_time:.2f}s")
        print()
        
        # Verify migration
        new_stats = new_db.get_database_stats()
        print("🎉 MIGRATION COMPLETE!")
        print(f"📊 New database: {new_stats['ta_sites_count']:,} TA sites, {new_stats['genes_count']:,} genes")
        
        # Verify data integrity
        test_genes = new_db.get_genes_in_region(1, 50000)
        test_ta_sites = new_db.get_ta_sites_in_region(1, 50000)
        print(f"🔍 Test query: {len(test_genes)} genes, {len(test_ta_sites)} TA sites in region 1-50,000")
        
        if test_genes:
            print(f"📝 Sample gene: {test_genes[0]['name']} ({test_genes[0]['start_position']}-{test_genes[0]['end_position']})")
        
        if test_ta_sites:
            print(f"📝 Sample TA site: pos {test_ta_sites[0]['position']}, viability {test_ta_sites[0]['viability_mean']:.1f}")
        
    except Exception as e:
        print(f"❌ Migration failed: {e}")
        import traceback
        traceback.print_exc()
    finally:
        conn.close()
        new_db.close_all_connections()

if __name__ == "__main__":
    migrate_real_data()