
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

    def _parse_header(self):
        if not self.file:
            self.file = open(self.filepath, "r", encoding="utf-8", errors="replace")

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
        if not self.file:
            self.file = open(self.filepath, "r", encoding="utf-8", errors="replace")

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

            # УБРАН параметр 'id' - используем только те, что есть в VariantRecord
            rec = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, info=info)
            # Добавляем остальные поля как атрибуты
            rec.qual = qual
            rec.filter = filter_
            yield rec

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
        chromosomes = set()
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            chromosomes.add(rec.chrom)

        if self.file:
            self.file.seek(current_pos)
        return sorted(chromosomes)

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        return chrom in self.get_chromosomes() and pos > 0

    def count_variants(self) -> int:
        cnt = 0
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for _ in self.read():
            cnt += 1

        if self.file:
            self.file.seek(current_pos)
        return cnt

    def stats_by_region(self) -> pd.DataFrame:
        from collections import Counter

        c = Counter()
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            c[rec.chrom] += 1

        if self.file:
            self.file.seek(current_pos)

        df = pd.DataFrame(list(c.items()), columns=["region", "variant_count"])
        return df

    def filter_by_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[VariantRecord]:
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if rec.chrom == chrom and start <= rec.pos <= end:
                yield rec

        if self.file:
            self.file.seek(current_pos)

    def filter_by_quality(self, min_qual: float) -> Iterator[VariantRecord]:
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if hasattr(rec, "qual") and rec.qual >= min_qual:
                yield rec

        if self.file:
            self.file.seek(current_pos)

    def get_records_in_region(
        self, chrom: str, start: int, end: int
    ) -> List[VariantRecord]:
        return list(self.filter_by_region(chrom, start, end))

    def filter_records(self, **filters) -> List[VariantRecord]:
        results = []
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            match = True

            if "chrom" in filters and rec.chrom != filters["chrom"]:
                match = False
            if (
                "min_qual" in filters
                and hasattr(rec, "qual")
                and rec.qual < filters["min_qual"]
            ):
                match = False
            if (
                "max_qual" in filters
                and hasattr(rec, "qual")
                and rec.qual > filters["max_qual"]
            ):
                match = False

            if match:
                results.append(rec)

        if self.file:
            self.file.seek(current_pos)
        return results

    def filter_by_quality(self, min_qual: float) -> Iterator[VariantRecord]:
        for rec in self.read():
            if hasattr(rec, 'qual') and rec.qual >= min_qual:
                yield rec

