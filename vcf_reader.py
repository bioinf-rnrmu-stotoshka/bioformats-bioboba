import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
from abstract import GenomicDataReader
from record import VariantRecord


class VcfReader(GenomicDataReader):
    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.header_lines: List[str] = []
        self.column_header: List[str] = []
        self._variants_cache = None

    def _parse_header(self):
        pos = self.file.tell()
        for line in self.file:
            if line.startswith("##"):
                self.header_lines.append(line.strip())
            elif line.startswith("#CHROM"):
                self.column_header = line.strip().lstrip("#").split("\t")
                break
        self.file.seek(pos)

    def get_header(self) -> List[str]:
        return self.header_lines

    def read(self) -> Iterator[VariantRecord]:
        if self._variants_cache is not None:
            yield from self._variants_cache
            return

        variants = []
        for line in self.file:
            if line.startswith("#"):
                continue
                
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            qual = float(parts[5]) if parts[5] != "." else 0.0
            filter_ = parts[6]
            info_str = parts[7] if len(parts) > 7 else ""
            
            info = self._parse_info(info_str)
            
            rec = VariantRecord(
                chrom=chrom, 
                pos=pos, 
                ref=ref, 
                alt=alt, 
                info=info
            )
            rec.qual = qual
            rec.filter = filter_
            
            variants.append(rec)
            yield rec

        self._variants_cache = variants

    @staticmethod
    def _parse_info(info_str: str) -> Dict[str, str]:
        d = {}
        if not info_str or info_str == ".":
            return d
        for kv in info_str.split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                d[k] = v
            else:
                d[kv] = "True"
        return d

    def get_chromosomes(self) -> List[str]:
        self._ensure_stats_calculated()
        return sorted(self._stats_cache['chromosomes'])  # ← ИСПРАВЛЕНО ИМЯ

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        return chrom in self.get_chromosomes() and pos > 0

    def get_records_in_region(self, chrom: str, start: int, end: int) -> Iterator[VariantRecord]:
        for rec in self.read():
            if rec.chrom == chrom and start <= rec.pos <= end:
                yield rec

    def filter_records(self, **filters) -> Iterator[VariantRecord]:
        for rec in self.read():
            match = True
            if 'chrom' in filters and rec.chrom != filters['chrom']:
                match = False
            if 'min_qual' in filters and hasattr(rec, 'qual') and rec.qual < filters['min_qual']:
                match = False
            if match:
                yield rec

    def count_variants(self) -> int:
        self._ensure_stats_calculated()
        return self._stats_cache['total_count']  # ← ИСПРАВЛЕНО ИМЯ

    def stats_by_region(self) -> pd.DataFrame:
        self._ensure_stats_calculated()
        counts = self._stats_cache['chromosome_counts']  # ← ИСПРАВЛЕНО ИМЯ
        df = pd.DataFrame(list(counts.items()), columns=["region", "variant_count"])
        return df.sort_values("variant_count", ascending=False)

    def filter_by_region(self, chrom: str, start: int, end: int) -> Iterator[VariantRecord]:
        return self.get_records_in_region(chrom, start, end)

    def filter_by_quality(self, min_qual: float) -> Iterator[VariantRecord]:
        for rec in self.read():
            if hasattr(rec, 'qual') and rec.qual >= min_qual:
                yield rec
