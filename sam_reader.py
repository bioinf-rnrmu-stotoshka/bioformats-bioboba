import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
import re
from abstract import GenomicDataReader
from record import AlignmentRecord


class SamReader(GenomicDataReader):
    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.header: Dict[str, List[str]] = {}

    # ===== Заголовки =====
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
        tag = parts[0]
        self.header.setdefault(tag, []).append("\t".join(parts[1:]))

    def get_header(self) -> Dict[str, List[str]]:
        """Возвращает словарь с заголовками SAM-файла"""
        return self.header

    # ===== Чтение выравниваний =====
    def read(self) -> Iterator[AlignmentRecord]:
        for line in self.file:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11 or fields[2] == "*":
                continue

            qname, flag, rname, pos, mapq, cigar = fields[:6]
            aligned_len = self._calc_aligned_length(cigar)
            end_pos = int(pos) + aligned_len - 1 if aligned_len > 0 else int(pos)

            rec = AlignmentRecord(
                id=qname, chrom=rname, start=int(pos), cigar=cigar, mapq=int(mapq)
            )
            rec.flag = int(flag)
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
        current_pos = self.file.tell()
        self.file.seek(0)
        count = sum(1 for _ in self.read())
        self.file.seek(current_pos)
        return count

    # ===== Статистика по хромосомам =====
    def stats_by_chromosome(self) -> pd.DataFrame:
        from collections import Counter

        current_pos = self.file.tell()
        self.file.seek(0)
        cnt = Counter(rec.chrom for rec in self.read())
        self.file.seek(current_pos)
        return pd.DataFrame(list(cnt.items()), columns=["chrom", "count"])

    # ===== Фильтр по региону =====
    def filter_by_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[AlignmentRecord]:
        current_pos = self.file.tell()
        self.file.seek(0)
        for rec in self.read():
            if rec.chrom == chrom and rec.end >= start and rec.start <= end:
                yield rec
        self.file.seek(current_pos)

