import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
import re  # модуль для работы с регулярными выражениями
from abstract import GenomicDataReader
from record import AlignmentRecord


class SamReader(GenomicDataReader):
    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.header: Dict[str, List[str]] = {}

    def _parse_header(self):
        pos = self.file.tell()  # запоминает текущую позицию в файле
        for line in self.file:
            if not line.startswith("@"):
                break
            self._parse_header_line(line)
        self.file.seek(pos)  # возврат к исходной позиции в файле

    def _parse_header_line(self, line: str):
        parts = line.strip().split("\t")
        tag = parts[0]  # первый элемент (тег заголовка)
        self.header.setdefault(tag, []).append("\t".join(parts[1:]))

    def get_header(self) -> Dict[str, List[str]]:
        return self.header

    def read(self) -> Iterator[AlignmentRecord]:
        for line in self.file:

            fields = line.strip().split("\t")
            if len(fields) < 11:  # проверяет, что в строке достаточно полей
                continue

            if fields[2] == "*":  # пропускает записи без хромосомы
                continue

            qname, flag, rname, pos, mapq, cigar = (
                fields[0],
                fields[1],
                fields[2],
                fields[3],
                fields[4],
                fields[5],
            )

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

        total = 0
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            if op in ("M", "D", "N", "=", "X"):
                total += int(length)
        return total

    def get_chromosomes(self) -> List[str]:
        #собирает уникальные хромосомы
        chromosomes = set()
        # СБРАСЫВАЕМ ПОЗИЦИЮ ФАЙЛА ПЕРЕД ЧТЕНИЕМ
        current_pos = self.file.tell()
        self.file.seek(0)

        for rec in self.read():
            if rec.chrom != "*":
                chromosomes.add(rec.chrom)

        # ВОССТАНАВЛИВАЕМ ПОЗИЦИЮ
        self.file.seek(current_pos)
        return sorted(chromosomes)

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Абстрактный метод - проверка координат"""
        return chrom in self.get_chromosomes() and pos > 0

    def calculate_coverage(self, chrom: str) -> Dict[int, int]:
        '''вычисляет покрытие по позициям'''
        coverage = {}
        current_pos = self.file.tell()
        self.file.seek(0)

        for rec in self.read():
            if rec.chrom == chrom and rec.cigar != "*":
                aligned_len = self._calc_aligned_length(rec.cigar)
                for i in range(rec.start, rec.start + aligned_len):
                    coverage[i] = coverage.get(i, 0) + 1

        self.file.seek(current_pos)
        return coverage

    def filter_alignments(self, flag: int) -> Iterator[AlignmentRecord]:
        '''Это фильтр-генератор, который отбирает только те выравнивания, 
        которые соответствуют определённым критериям, заданным через битовые флаги'''
        current_pos = self.file.tell()
        self.file.seek(0)

        for rec in self.read():
            if hasattr(rec, "flag") and rec.flag & flag:
                #проверяет, есть ли у объекта rec атрибут flag
                #побитовая операция И между флагами текущей записи и маской фильтра
                yield rec

        self.file.seek(current_pos)

    def count_alignments(self) -> int:
        '''Считает общее количество выравниваний'''
        cnt = 0
        current_pos = self.file.tell()
        self.file.seek(0)

        for _ in self.read():
            cnt += 1

        self.file.seek(current_pos)
        return cnt

    def stats_by_chromosome(self) -> pd.DataFrame:
        '''Статистика по хромосомам через pandas'''
        from collections import Counter #специальный словарь для подсчета
                                                    
        cnt = Counter()
        current_pos = self.file.tell()
        self.file.seek(0)

        for rec in self.read():
            cnt[rec.chrom] += 1

        self.file.seek(current_pos)
        df = pd.DataFrame(list(cnt.items()), columns=["chrom", "count"])
        return df

    def filter_by_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[AlignmentRecord]:
        '''Фильтрует по геномной области'''
        current_pos = self.file.tell()
        self.file.seek(0)

        for rec in self.read():
            if rec.chrom == chrom and rec.end >= start and rec.start <= end:
                yield rec

        self.file.seek(current_pos)

    def get_records_in_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[AlignmentRecord]:
        """Абстрактный метод - записи в регионе"""
        return list(self.filter_by_region(chrom, start, end))

    def filter_records(self, **filters) -> Iterator[AlignmentRecord]:
        """Абстрактный метод - фильтрация записей"""
        results = []
        current_pos = self.file.tell()
        self.file.seek(0)

        for rec in self.read():
            match = True

            if "chrom" in filters and rec.chrom != filters["chrom"]:
                match = False
            if "min_mapq" in filters and rec.mapq < filters["min_mapq"]:
                match = False
            if (
                "flag" in filters
                and hasattr(rec, "flag")
                and not (rec.flag & filters["flag"])
            ):
                match = False

            if match:
                results.append(rec)

        self.file.seek(current_pos)

