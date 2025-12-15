#!/usr/bin/env python3
"""Backfill script om 'symbol' kolom in genes.db te vullen vanuit TSV.

Gebruik:
  python scripts/backfill_symbols.py --tsv data/ncbi_dataset.tsv --data-dir data

Werking:
  - Leest TSV (kolommen: Accession Begin End Chromosome Orientation Name Symbol Gene ID Gene Type ... Locus tag)
  - Voor iedere rij waar Symbol niet leeg is wordt UPDATE genes SET symbol=? WHERE locus_tag=? uitgevoerd.
  - Rapportage hoeveel records geüpdatet zijn, hoeveel gene rows geen match hadden of al gevuld waren.

Idempotent: herhaald draaien overschrijft alleen lege of verschillende symbol waarden (optioneel --force om altijd te overschrijven).
"""
import csv
import argparse
import sqlite3
from pathlib import Path


def backfill(tsv_path: Path, data_dir: Path, force: bool = False):
    genes_db = data_dir / 'genes.db'
    if not genes_db.exists():
        raise SystemExit(f"genes.db niet gevonden in {data_dir}")
    conn = sqlite3.connect(str(genes_db))
    cur = conn.cursor()
    # Check kolom
    cur.execute("PRAGMA table_info(genes)")
    cols = {r[1] for r in cur.fetchall()}
    if 'symbol' not in cols:
        raise SystemExit("Kolom 'symbol' ontbreekt; run de app zodat migratie draait")

    updated = 0
    skipped_no_symbol = 0
    skipped_already = 0
    not_found = 0

    with open(tsv_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Normaliseer kolomnamen
        field_map = {k.lower(): k for k in reader.fieldnames}
        # Vereiste kolommen
        required = ['locus tag', 'symbol']
        for r in required:
            if r not in field_map:
                raise SystemExit(f"Vereiste kolom ontbreekt in TSV: {r}")
        for row in reader:
            raw_symbol = row.get(field_map['symbol'], '').strip()
            locus_tag = row.get(field_map['locus tag'], '').strip()
            if not locus_tag:
                continue
            if not raw_symbol:
                skipped_no_symbol += 1
                continue
            # Haal huidige waarde
            cur.execute("SELECT symbol FROM genes WHERE locus_tag=?", (locus_tag,))
            res = cur.fetchone()
            if not res:
                not_found += 1
                continue
            current_symbol = res[0]
            if current_symbol and current_symbol.strip() == raw_symbol and not force:
                skipped_already += 1
                continue
            cur.execute("UPDATE genes SET symbol=? WHERE locus_tag=?", (raw_symbol, locus_tag))
            updated += 1
    conn.commit()
    print(f"Backfill klaar: updated={updated}, skipped_no_symbol={skipped_no_symbol}, skipped_already={skipped_already}, not_found={not_found}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tsv', required=True, help='Pad naar ncbi_dataset.tsv')
    ap.add_argument('--data-dir', default='data', help='Directory met genes.db')
    ap.add_argument('--force', action='store_true', help='Overschrijf bestaande symbol waarden')
    args = ap.parse_args()
    backfill(Path(args.tsv), Path(args.data_dir), force=args.force)

if __name__ == '__main__':
    main()
