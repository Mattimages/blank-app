#!/usr/bin/env python3
"""
High-Performance Multi-Database Genomic Database Architecture
============================================================

Normalized database schema with separate databases for maximum speed:
- genes.db: Gene information with optimized indexes
- ta_sites.db: TA site positions and viability data linked to genes
- orthology.db: Gene ortholog relationships

Features:
- Connection pooling for speed
- Optimized indexes on all query columns
- Foreign key constraints for data integrity
- Batch operations for bulk inserts
- Prepared statements for repeated queries
"""

import sqlite3
import os
import logging
from typing import List, Dict, Any, Optional, Tuple
from contextlib import contextmanager
import threading
from collections import defaultdict
import time

class HighPerformanceGenomicDB:
    """Multi-database genomic data manager optimized for speed."""
    
    def __init__(self, data_dir="data"):
        self.data_dir = data_dir
        os.makedirs(data_dir, exist_ok=True)
        
        # Database paths
        self.genes_db_path = os.path.join(data_dir, "genes.db")
        self.ta_sites_db_path = os.path.join(data_dir, "ta_sites.db") 
        self.orthology_db_path = os.path.join(data_dir, "orthology.db")
        
        # Connection pools for speed
        self._genes_connections = {}
        self._ta_sites_connections = {}
        self._orthology_connections = {}
        self._lock = threading.Lock()
        
        # Initialize databases
        self._initialize_databases()

        # Query instrumentation counters
        self.query_stats = {
            'genes': 0,
            'ta_sites': 0,
            'orthologs': 0,
            'last_reset': time.time()
        }

    def reset_query_stats(self):
        """Reset query counters (used for performance diagnostics)."""
        self.query_stats.update({
            'genes': 0,
            'ta_sites': 0,
            'orthologs': 0,
            'last_reset': time.time()
        })

    def get_query_stats(self) -> Dict[str, Any]:
        """Return current query counters and elapsed time since last reset."""
        elapsed = time.time() - self.query_stats.get('last_reset', time.time())
        return {**self.query_stats, 'elapsed_seconds': elapsed}
        
    def _get_connection(self, db_type: str) -> sqlite3.Connection:
        """Get thread-safe database connection with connection pooling."""
        thread_id = threading.get_ident()
        
        with self._lock:
            if db_type == "genes":
                if thread_id not in self._genes_connections:
                    conn = sqlite3.connect(self.genes_db_path)
                    conn.row_factory = sqlite3.Row
                    conn.execute("PRAGMA journal_mode=WAL")  # Speed optimization
                    conn.execute("PRAGMA synchronous=NORMAL")  # Speed vs safety balance
                    conn.execute("PRAGMA cache_size=10000")  # Large cache for speed
                    conn.execute("PRAGMA temp_store=MEMORY")  # Use memory for temp data
                    self._genes_connections[thread_id] = conn
                return self._genes_connections[thread_id]
                
            elif db_type == "ta_sites":
                if thread_id not in self._ta_sites_connections:
                    conn = sqlite3.connect(self.ta_sites_db_path)
                    conn.row_factory = sqlite3.Row
                    conn.execute("PRAGMA journal_mode=WAL")
                    conn.execute("PRAGMA synchronous=NORMAL")
                    conn.execute("PRAGMA cache_size=10000")
                    conn.execute("PRAGMA temp_store=MEMORY")
                    self._ta_sites_connections[thread_id] = conn
                return self._ta_sites_connections[thread_id]
                
            elif db_type == "orthology":
                if thread_id not in self._orthology_connections:
                    conn = sqlite3.connect(self.orthology_db_path)
                    conn.row_factory = sqlite3.Row
                    conn.execute("PRAGMA journal_mode=WAL")
                    conn.execute("PRAGMA synchronous=NORMAL") 
                    conn.execute("PRAGMA cache_size=10000")
                    conn.execute("PRAGMA temp_store=MEMORY")
                    self._orthology_connections[thread_id] = conn
                return self._orthology_connections[thread_id]
    
    def _ensure_column(self, conn: sqlite3.Connection, table: str, column: str, ddl: str):
        """Idempotent kolom toevoeging als deze ontbreekt."""
        cur = conn.cursor()
        cur.execute(f"PRAGMA table_info({table})")
        cols = {r[1] for r in cur.fetchall()}
        if column not in cols:
            logging.warning(f"Migratie: voeg kolom '{column}' toe aan tabel '{table}'.")
            cur.execute(f"ALTER TABLE {table} ADD COLUMN {ddl}")
            conn.commit()

    def _initialize_databases(self):
        """Create optimized database schemas with indexes (en voer eenvoudige migraties uit)."""
        
        # GENES DATABASE - Core gene information
        genes_conn = self._get_connection("genes")
        # Basis schema (zonder assumptie dat symbol al bestaat)
        genes_conn.executescript("""
            CREATE TABLE IF NOT EXISTS genes (
                locus_tag TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                start_position INTEGER NOT NULL,
                end_position INTEGER NOT NULL,
                strand TEXT NOT NULL,
                gene_type TEXT,
                product TEXT,
                symbol TEXT,
                chromosome TEXT DEFAULT 'main',
                UNIQUE(start_position, end_position, chromosome)
            );
        """)
        genes_conn.commit()
        # Migratie: zorg dat symbol kolom aanwezig is
        self._ensure_column(genes_conn, 'genes', 'symbol', "symbol TEXT")
        # Indexen (idempotent)
        genes_conn.executescript("""
            CREATE INDEX IF NOT EXISTS idx_genes_position ON genes(start_position, end_position);
            CREATE INDEX IF NOT EXISTS idx_genes_name ON genes(name);
            CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(symbol);
            CREATE INDEX IF NOT EXISTS idx_genes_region ON genes(chromosome, start_position, end_position);
            CREATE INDEX IF NOT EXISTS idx_genes_strand ON genes(strand);
            CREATE INDEX IF NOT EXISTS idx_genes_locus_tag ON genes(locus_tag);
            CREATE INDEX IF NOT EXISTS idx_genes_position_strand ON genes(start_position, end_position, strand);
        """)
        genes_conn.commit()
        
        # TA SITES DATABASE - Transposon insertion data linked to genes
        ta_conn = self._get_connection("ta_sites")
        ta_conn.executescript("""
            CREATE TABLE IF NOT EXISTS ta_sites (
                ta_site_id INTEGER PRIMARY KEY AUTOINCREMENT,
                position INTEGER NOT NULL UNIQUE,
                locus_tag TEXT,
                viability_mean REAL,
                viability_std REAL,
                viability_values TEXT,  -- JSON array of raw values
                viable_colonies INTEGER,
                total_colonies INTEGER DEFAULT 200,
                p_value REAL DEFAULT 0.05,
                chromosome TEXT DEFAULT 'main',
                FOREIGN KEY (locus_tag) REFERENCES genes(locus_tag)
            );
            
            -- Ultra-fast position-based queries
            CREATE INDEX IF NOT EXISTS idx_ta_sites_position ON ta_sites(position);
            CREATE INDEX IF NOT EXISTS idx_ta_sites_locus_tag ON ta_sites(locus_tag);
            CREATE INDEX IF NOT EXISTS idx_ta_sites_region ON ta_sites(chromosome, position);
            CREATE INDEX IF NOT EXISTS idx_ta_sites_viability ON ta_sites(viability_mean);
            
            -- Compound index for range queries
            CREATE INDEX IF NOT EXISTS idx_ta_sites_position_viability ON ta_sites(position, viability_mean);
        """)
        ta_conn.commit()
        
        # ORTHOLOGY DATABASE - Gene relationships and comparative data
        orth_conn = self._get_connection("orthology")
        orth_conn.executescript("""
            CREATE TABLE IF NOT EXISTS orthologs (
                ortholog_id INTEGER PRIMARY KEY AUTOINCREMENT,
                locus_tag TEXT NOT NULL,
                ortholog_locus_tag TEXT,
                organism TEXT,
                ortholog_name TEXT,
                similarity_score REAL,
                e_value REAL,
                relationship_type TEXT DEFAULT 'ortholog',  -- ortholog, paralog, etc.
                FOREIGN KEY (locus_tag) REFERENCES genes(locus_tag)
            );
            
            -- Indexes for ortholog lookups
            CREATE INDEX IF NOT EXISTS idx_orthologs_locus_tag ON orthologs(locus_tag);
            CREATE INDEX IF NOT EXISTS idx_orthologs_organism ON orthologs(organism);
            CREATE INDEX IF NOT EXISTS idx_orthologs_similarity ON orthologs(similarity_score);
        """)
        orth_conn.commit()
        
        print("✅ High-performance database schema created!")

    def get_database_stats(self) -> Dict[str, int]:
        """Get comprehensive database statistics."""
        stats = {}
        
        # Gene counts
        genes_conn = self._get_connection("genes")
        cursor = genes_conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM genes")
        stats['genes_count'] = cursor.fetchone()[0]
        
        # TA site counts
        ta_conn = self._get_connection("ta_sites")
        cursor = ta_conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM ta_sites")
        stats['ta_sites_count'] = cursor.fetchone()[0]
        
        # Ortholog counts
        orth_conn = self._get_connection("orthology")
        cursor = orth_conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM orthologs")
        stats['orthologs_count'] = cursor.fetchone()[0]
        
        return stats

    def bulk_insert_genes(self, genes_data: List[Dict]) -> List[int]:
        """High-speed bulk insert for genes with duplicate handling."""
        start_time = time.time()
        genes_conn = self._get_connection("genes")
        cursor = genes_conn.cursor()
        
        # Prepare batch insert
        gene_ids = []
        insert_data = []
        
        for gene in genes_data:
            insert_data.append((
                gene.get('locus_tag', ''),  # Primary key
                gene.get('name', ''),
                gene.get('start_position', 0),
                gene.get('end_position', 0),
                gene.get('strand', '+'),
                gene.get('gene_type', ''),
                gene.get('product', ''),
                gene.get('symbol') or gene.get('Symbol') or '',
                gene.get('chromosome', 'main')
            ))
        
        try:
            # Use INSERT OR IGNORE for speed (skip duplicates)
            cursor.executemany("""
                INSERT OR IGNORE INTO genes 
                (locus_tag, name, start_position, end_position, strand, gene_type, product, symbol, chromosome)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, insert_data)
            
            # Get the inserted locus_tags
            for gene in genes_data:
                locus_tag = gene.get('locus_tag', '')
                if locus_tag:
                    gene_ids.append(locus_tag)
            
            genes_conn.commit()
            
        except Exception as e:
            genes_conn.rollback()
            raise e
        
        elapsed = time.time() - start_time
        print(f"⚡ Bulk inserted {len(gene_ids)} genes in {elapsed:.3f}s ({len(gene_ids)/elapsed:.0f} genes/sec)")
        return gene_ids

    def bulk_insert_ta_sites(self, ta_sites_data: List[Dict]):
        """High-speed bulk insert for TA sites with gene linking."""
        start_time = time.time()
        ta_conn = self._get_connection("ta_sites")
        genes_conn = self._get_connection("genes")
        
        cursor = ta_conn.cursor()
        genes_cursor = genes_conn.cursor()
        
        # Prepare batch insert with gene linking
        insert_data = []
        
        for ta_site in ta_sites_data:
            position = ta_site.get('position', 0)
            
            # Find gene that contains this TA site
            genes_cursor.execute("""
                SELECT locus_tag FROM genes 
                WHERE start_position <= ? AND end_position >= ?
                ORDER BY (end_position - start_position) ASC  -- Prefer smaller genes
                LIMIT 1
            """, (position, position))
            
            gene_result = genes_cursor.fetchone()
            locus_tag = gene_result[0] if gene_result else None
            
            insert_data.append((
                position,
                locus_tag,
                ta_site.get('viability_mean'),
                ta_site.get('viability_std'),
                ta_site.get('viability_values'),  # JSON string
                ta_site.get('viable_colonies'),
                ta_site.get('total_colonies', 200),
                ta_site.get('p_value', 0.05),
                ta_site.get('chromosome', 'main')
            ))
        
        try:
            cursor.executemany("""
                INSERT OR REPLACE INTO ta_sites 
                (position, locus_tag, viability_mean, viability_std, viability_values, 
                 viable_colonies, total_colonies, p_value, chromosome)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, insert_data)
            
            ta_conn.commit()
            
        except Exception as e:
            ta_conn.rollback()
            raise e
        
        elapsed = time.time() - start_time
        print(f"⚡ Bulk inserted {len(insert_data)} TA sites in {elapsed:.3f}s ({len(insert_data)/elapsed:.0f} sites/sec)")

    def get_genes_in_region(self, start_pos: int, end_pos: int, chromosome: str = 'main') -> List[Dict]:
        """Ultra-fast gene retrieval with optimized indexing."""
        genes_conn = self._get_connection("genes")
        cursor = genes_conn.cursor()

        # Instrumentation
        try:
            self.query_stats['genes'] += 1
        except Exception:
            pass
        
        cursor.execute("""
            SELECT locus_tag, name, start_position, end_position, strand, gene_type, product, symbol
            FROM genes 
            WHERE chromosome = ? AND start_position <= ? AND end_position >= ?
            ORDER BY start_position
        """, (chromosome, end_pos, start_pos))
        
        return [dict(row) for row in cursor.fetchall()]

    def get_ta_sites_in_region(self, start_pos: int, end_pos: int, chromosome: str = 'main') -> List[Dict]:
        """Ultra-fast TA site retrieval with gene information."""
        ta_conn = self._get_connection("ta_sites")
        genes_conn = self._get_connection("genes")
        
        # Use optimized query with potential JOIN for gene info
        cursor = ta_conn.cursor()

        # Instrumentation
        try:
            self.query_stats['ta_sites'] += 1
        except Exception:
            pass
        cursor.execute("""
            SELECT ta_site_id, position, locus_tag, viability_mean, viability_std, 
                   viability_values, viable_colonies, total_colonies, p_value
            FROM ta_sites 
            WHERE chromosome = ? AND position >= ? AND position <= ?
            ORDER BY position
        """, (chromosome, start_pos, end_pos))
        
        ta_sites = []
        for row in cursor.fetchall():
            ta_site = dict(row)
            
            # Add gene info if linked
            if ta_site['locus_tag']:
                genes_cursor = genes_conn.cursor()
                genes_cursor.execute("SELECT name FROM genes WHERE locus_tag = ?", (ta_site['locus_tag'],))
                gene_result = genes_cursor.fetchone()
                if gene_result:
                    ta_site['gene_name'] = gene_result[0]
            
            ta_sites.append(ta_site)
        
        return ta_sites

    def search_gene_by_name(self, gene_name: str) -> Optional[Dict]:
        """Search gene by (priority): symbol, locus_tag, name, then partial matches.
        Returns dict including symbol if present."""
        genes_conn = self._get_connection("genes")
        cursor = genes_conn.cursor()
        # Normalize query
        q = gene_name.strip()
        if not q:
            return None
        # Exact symbol
        cursor.execute("""
            SELECT locus_tag, name, start_position, end_position, strand, gene_type, product, symbol
            FROM genes WHERE symbol = ? LIMIT 1
        """, (q,))
        row = cursor.fetchone()
        if row:
            return dict(row)
        # Exact locus_tag
        cursor.execute("""
            SELECT locus_tag, name, start_position, end_position, strand, gene_type, product, symbol
            FROM genes WHERE locus_tag = ? LIMIT 1
        """, (q,))
        row = cursor.fetchone()
        if row:
            return dict(row)
        # Exact name
        cursor.execute("""
            SELECT locus_tag, name, start_position, end_position, strand, gene_type, product, symbol
            FROM genes WHERE name = ? LIMIT 1
        """, (q,))
        row = cursor.fetchone()
        if row:
            return dict(row)
        # Partial (symbol or name or locus_tag)
        like = f"%{q}%"
        cursor.execute("""
            SELECT locus_tag, name, start_position, end_position, strand, gene_type, product, symbol
            FROM genes
            WHERE symbol LIKE ? OR name LIKE ? OR locus_tag LIKE ?
            ORDER BY start_position
            LIMIT 1
        """, (like, like, like))
        row = cursor.fetchone()
        return dict(row) if row else None

    def search_genes_multi(self, query: str, limit: int = 20) -> List[Dict]:
        """Case-insensitive multi-result search over symbol, locus_tag, name.
        Returns up to 'limit' matches ordered by start_position."""
        q = (query or '').strip()
        if not q:
            return []
        genes_conn = self._get_connection("genes")
        cursor = genes_conn.cursor()
        like = f"%{q}%"
        cursor.execute(
            """
            SELECT locus_tag, name, start_position, end_position, strand, gene_type, product, symbol
            FROM genes
            WHERE (symbol LIKE ? COLLATE NOCASE)
               OR (name LIKE ? COLLATE NOCASE)
               OR (locus_tag LIKE ? COLLATE NOCASE)
            ORDER BY start_position
            LIMIT ?
            """,
            (like, like, like, limit)
        )
        return [dict(r) for r in cursor.fetchall()]

    def get_gene_orthologs(self, gene_name: str) -> List[Dict]:
        """Get ortholog information for a gene."""
        # First get gene_id
        gene = self.search_gene_by_name(gene_name)
        if not gene:
            return []
        
        orth_conn = self._get_connection("orthology") 
        cursor = orth_conn.cursor()
        
        cursor.execute("""
            SELECT ortholog_gene_id, organism, ortholog_name, similarity_score, 
                   e_value, relationship_type
            FROM orthologs WHERE gene_id = ?
            ORDER BY similarity_score DESC
        """, (gene['gene_id'],))
        
        return [dict(row) for row in cursor.fetchall()]

    def close_all_connections(self):
        """Clean up all database connections."""
        with self._lock:
            for conn_dict in [self._genes_connections, self._ta_sites_connections, self._orthology_connections]:
                for conn in conn_dict.values():
                    conn.close()
                conn_dict.clear()

# Compatibility wrapper for existing code
GenomicDatabase = HighPerformanceGenomicDB