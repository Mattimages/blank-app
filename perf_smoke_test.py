#!/usr/bin/env python3
"""Lightweight performance & correctness smoke test for genome DB + TA rendering logic.
Run: python perf_smoke_test.py
"""
import time
from src.database_v2 import HighPerformanceGenomicDB


def main():
    db = HighPerformanceGenomicDB(data_dir='data')
    stats = db.get_database_stats()
    print(f"DB Stats: genes={stats['genes_count']} ta_sites={stats['ta_sites_count']} orthologs={stats['orthologs_count']}")

    # Test gene region query timing
    regions = [ (0, 100000), (500000, 600000), (0, 1000000) ]
    for start,end in regions:
        t0=time.time(); genes=db.get_genes_in_region(start,end); dt=time.time()-t0
        print(f"genes_in_region {start}-{end}: {len(genes)} genes in {dt*1000:.1f} ms")

    # Test TA sites region query timing
    for start,end in regions:
        t0=time.time(); ta_sites=db.get_ta_sites_in_region(start,end); dt=time.time()-t0
        print(f"ta_sites_in_region {start}-{end}: {len(ta_sites)} sites in {dt*1000:.1f} ms")

    # Simulate performance mode trigger thresholds used in UI
    MAX_POINTS = 8000
    counts = [2000, 10000, 14000]
    for c in counts:
        perf_mode=False; reason=[]
        if c > MAX_POINTS*1.5:
            perf_mode=True; reason.append(f"points>{MAX_POINTS*1.5}")
        if c > 12000:
            perf_mode=True; reason.append("gradient_off")
        print(f"simulate count={c}: perf_mode={perf_mode} reason={','.join(reason)}")

if __name__ == '__main__':
    main()
