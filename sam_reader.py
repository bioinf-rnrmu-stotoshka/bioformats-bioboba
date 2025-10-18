import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
import re
from collections import Counter
from abstract import GenomicDataReader
from record import AlignmentRecord


class SamReader(GenomicDataReader):
    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self._header_parsed = False
        self.header: Dict[str, List[str]] = {}

    def _ensure_header_parsed(self):
        if not self._header_parsed:
            self._parse_header()
            self._header_parsed = True

    def _parse_header(self):
        """Читает заголовки SAM-файла и сохраняет их в self.header"""
        self.file.seek(0)
        for line in self.file:
            if not line.startswith("@"):
                break
            self._parse_header_line(line)
        self.file.seek(0)

    def _parse_header_line(self, line: str):
        parts = line.strip().split("\t")
        if parts:
            tag = parts[0]
            self.header.setdefault(tag, []).append("\t".join(parts[1:]))

    def get_header(self) -> Dict[str, List[str]]:
        self._ensure_header_parsed()
        return self.header

    def get_header_group(self, group_tag: str) -> List[str]:
        """Получить конкретную группу заголовков, например '@PG' или '@RG'."""
        self._ensure_header_parsed()
        return self.header.get(group_tag, [])

    # ===== Чтение выравниваний =====
    def read(self) -> Iterator[AlignmentRecord]:
        self.file.seek(0)
        for line in self.file:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11 or fields[2] == "*":
                continue

            qname, flag, rname, pos, mapq, cigar = fields[:6]
            try:
                pos_int = int(pos)
                mapq_int = int(mapq) if mapq != "*" else 0
                flag_int = int(flag)
            except ValueError:
                continue

            aligned_len = self._calc_aligned_length(cigar)
            end_pos = pos_int + aligned_len - 1 if aligned_len > 0 else pos_int

            rec = AlignmentRecord(
                id=qname, chrom=rname, start=pos_int, cigar=cigar, mapq=mapq_int
            )
            rec.flag = flag_int
            rec.end = end_pos
            yield rec

    @staticmethod
    def _calc_aligned_length(cigar: str) -> int:
        if not cigar or cigar == "*":
            return 0
        return sum(
            int(length) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
            if op in ("M", "D", "N", "=", "X")
        )

    # ===== Количество выравниваний =====
    def count_alignments(self) -> int:
        self.file.seek(0)
        count = 0
        for line in self.file:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 11 and fields[2] != "*":
                count += 1
        self.file.seek(0)
        return count

    # ===== Статистика по хромосомам =====
    def stats_by_chromosome(self) -> pd.DataFrame:
        chrom_counts = Counter()
        self.file.seek(0)
        for rec in self.read():
            chrom_counts[rec.chrom] += 1
        self.file.seek(0)
        if not chrom_counts:
            return pd.DataFrame(columns=["chrom", "count"])
        return pd.DataFrame(list(chrom_counts.items()), columns=["chrom", "count"])

    # ===== Фильтр по региону =====
    def filter_by_region(self, chrom: str, start: int, end: int) -> Iterator[AlignmentRecord]:
        if start > end:
            raise ValueError("start must be <= end")
        self.file.seek(0)
        for rec in self.read():
            if rec.chrom == chrom and rec.start <= end and rec.end >= start:
                yield rec
        self.file.seek(0)