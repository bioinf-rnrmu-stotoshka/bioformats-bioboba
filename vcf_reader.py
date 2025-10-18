import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
from collections import Counter
from abstract import GenomicDataReader
from record import VariantRecord


class VcfReader(GenomicDataReader):
    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.header_lines: List[str] = []
        self.column_header: List[str] = []

    def _parse_header(self) -> None:
        """
        Парсит заголовок VCF-файла. Вызывается автоматически при входе в контекст.
        Предполагается, что self.file уже открыт.
        """
        if self._header_parsed:
            return

        # Сохраняем текущую позицию (хотя в __enter__ мы в начале файла)
        start_pos = self.file.tell()
        self.file.seek(0)

        for line in self.file:
            if line.startswith("##"):
                self.header_lines.append(line.strip())
            elif line.startswith("#CHROM"):
                self.column_header = line.strip().lstrip("#").split("\t")
                break
        else:
            # Если не найден #CHROM — это некорректный VCF
            raise ValueError("Заголовок VCF не содержит строки #CHROM")

        # Возвращаемся в начало файла для последующего чтения записей
        self.file.seek(0)
        self._header_parsed = True

    def get_header(self) -> List[str]:
        """Возвращает все мета-заголовки (##...)."""
        if not self._header_parsed:
            self._parse_header()
        return self.header_lines

    def get_header_group(self, key: str) -> List[str]:
        """
        Возвращает список строк заголовка, начинающихся с `##{key}=`.
        Пример: get_header_group("INFO") → все ##INFO=...
        """
        if not self._header_parsed:
            self._parse_header()
        prefix = f"##{key}="
        return [line for line in self.header_lines if line.startswith(prefix)]

    def read(self) -> Iterator[VariantRecord]:
        """Генератор записей вариантов, пропускающий заголовки."""
        if not self._header_parsed:
            self._parse_header()

        for line in self.file:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            try:
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                filter_ = parts[6] if len(parts) > 6 else "."
                info_str = parts[7] if len(parts) > 7 else ""
                info = self._parse_info(info_str)

                rec = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, info=info)
                rec.filter = filter_
                yield rec
            except (ValueError, IndexError):
                continue  # Пропускаем некорректные строки

    @staticmethod
    def _parse_info(info_str: str) -> Dict[str, str]:
        if not info_str or info_str == ".":
            return {}
        d = {}
        for kv in info_str.split(";"):
            if not kv:
                continue
            if "=" in kv:
                k, v = kv.split("=", 1)
                d[k] = v
            else:
                d[kv] = "True"
        return d

    def count_variants(self) -> int:
        """Возвращает общее количество вариантов в файле."""
        if not self._header_parsed:
            self._parse_header()

        self.file.seek(0)
        count = 0
        for line in self.file:
            if not line.startswith("#"):
                count += 1
        self.file.seek(0)
        return count

    def stats_by_region(self) -> pd.DataFrame:
        """
        Возвращает статистику: количество вариантов по хромосомам (регион = хромосома).
        """
        if not self._header_parsed:
            self._parse_header()

        counter = Counter()
        self.file.seek(0)
        for rec in self.read():
            counter[rec.chrom] += 1
        self.file.seek(0)

        if not counter:
            return pd.DataFrame(columns=["region", "variant_count"])
        df = pd.DataFrame(counter.items(), columns=["region", "variant_count"])
        return df.sort_values("region").reset_index(drop=True)

    def filter_by_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[VariantRecord]:
        """
        Возвращает варианты, попадающие в геномный регион [start, end] (1-based, inclusive).
        Аналог bedtools intersect.
        """
        if start > end:
            raise ValueError("start must be <= end")
        if not self._header_parsed:
            self._parse_header()

        self.file.seek(0)
        for rec in self.read():
            if rec.chrom == chrom and start <= rec.pos <= end:
                yield rec
        self.file.seek(0)