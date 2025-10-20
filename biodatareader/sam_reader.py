import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
import re
from collections import Counter
from abstract import GenomicDataReader
from record import AlignmentRecord


class SamReader(GenomicDataReader):
    """
    Реализация ридера для чтения SAM-файлов (Sequence Alignment/Map).

    Поддерживает:
        - Парсинг заголовков (строки, начинающиеся с '@').
        - Итеративное чтение выравниваний как объектов AlignmentRecord.
        - Подсчёт общего числа выравниваний.
        - Сбор статистики по хромосомам с использованием pandas.
        - Фильтрацию выравниваний по геномному региону.

    Attributes:
        filepath (Path): Путь к SAM-файлу.
        file (file object or None): Открытый файловый дескриптор.
        _header_parsed (bool): Флаг, указывающий, был ли уже распарсен заголовок.
        header (Dict[str, List[str]]): Словарь заголовков, где ключ — тег (например, '@SQ'),
            значение — список строк содержимого без тега.
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализирует SamReader с указанным путём к файлу.

        Args:
            filepath (str | Path): Путь к SAM-файлу.
        """
        super().__init__(filepath)
        self._header_parsed = False
        self.header: Dict[str, List[str]] = {}

    def _ensure_header_parsed(self) -> None:
        """
        Гарантирует, что заголовок был распарсен.

        Если парсинг ещё не выполнялся, вызывает _parse_header() и устанавливает флаг.
        """
        if not self._header_parsed:
            self._parse_header()
            self._header_parsed = True

    def _parse_header(self) -> None:
        """
        Читает и парсит заголовочные строки SAM-файла.

        Перемещает файловый указатель в начало, читает строки до первой
        незаголовочной (не начинающейся с '@'), сохраняя заголовки в self.header.
        После завершения возвращает указатель в начало файла.
        """
        self.file.seek(0)
        for line in self.file:
            if not line.startswith("@"):
                break
            self._parse_header_line(line)
        self.file.seek(0)

    def _parse_header_line(self, line: str) -> None:
        """
        Парсит одну заголовочную строку и сохраняет её в self.header.

        Разделяет строку по табуляции, первый элемент — тег (например, '@SQ'),
        остальные — содержимое. Тег используется как ключ в словаре header.

        Args:
            line (str): Одна строка из SAM-файла, начинающаяся с '@'.
        """
        parts = line.strip().split("\t")
        if parts:
            tag = parts[0]
            self.header.setdefault(tag, []).append("\t".join(parts[1:]))

    def get_header(self) -> Dict[str, List[str]]:
        """
        Возвращает все заголовки SAM-файла в виде словаря.

        Гарантирует, что заголовок был распарсен перед возвратом.

        Returns:
            Dict[str, List[str]]: Словарь, где ключи — теги заголовков (например, '@SQ'),
            значения — списки строк содержимого (без тега).
        """
        self._ensure_header_parsed()
        return self.header

    def get_header_group(self, group_tag: str) -> List[str]:
        """
        Возвращает список записей для указанной группы заголовков.

        Например, '@SQ' для информации о последовательностях или '@RG' для групп ридов.

        Args:
            group_tag (str): Тег заголовочной группы (например, '@SQ', '@PG').

        Returns:
            List[str]: Список строк содержимого для данной группы. Пустой список,
            если группа отсутствует.
        """
        self._ensure_header_parsed()
        return self.header.get(group_tag, [])

    def read(self) -> Iterator[AlignmentRecord]:
        """
        Итеративно читает выравнивания из SAM-файла и возвращает объекты AlignmentRecord.

        Пропускает заголовочные строки и некорректные записи (менее 11 полей или '*'
        в поле RNAME). Вычисляет конечную позицию выравнивания на основе CIGAR.

        Yields:
            AlignmentRecord: Объект выравнивания с заполненными атрибутами.

        Note:
            Метод сбрасывает файловый указатель в начало перед чтением.
        """
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
        """
        Вычисляет длину выравнивания на референсе по CIGAR-строке.

        Учитывает операции, которые "потребляют" референс: M, D, N, =, X.
        Игнорирует операции вставок (I), мягких (S) и жёстких (H) клипов и др.

        Args:
            cigar (str): CIGAR-строка (например, "100M2D30M").

        Returns:
            int: Длина выравнивания в позициях референса.

        Example:
            >>> SamReader._calc_aligned_length("50M10D20M")
            80
            >>> SamReader._calc_aligned_length("30I")  # только вставка — не влияет на референс
            0
        """
        if not cigar or cigar == "*":
            return 0
        return sum(
            int(length) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
            if op in ("M", "D", "N", "=", "X")
        )

    def count_alignments(self) -> int:
        """
        Подсчитывает общее количество валидных выравниваний в SAM-файле.

        Пропускает заголовки и строки с '*' в поле RNAME или менее чем 11 полями.

        Returns:
            int: Количество валидных выравниваний.

        Note:
            Метод сбрасывает файловый указатель в начало после подсчёта.
        """
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

    def stats_by_chromosome(self) -> pd.DataFrame:
        """
        Собирает статистику количества выравниваний по хромосомам.

        Возвращает pandas DataFrame с двумя колонками: 'chrom' и 'count'.

        Returns:
            pd.DataFrame: Таблица с количеством выравниваний на каждую хромосому.
            Если выравниваний нет, возвращает пустой DataFrame с указанными колонками.

        Note:
            Метод сбрасывает файловый указатель в начало после сбора статистики.
        """
        chrom_counts = Counter()
        self.file.seek(0)
        for rec in self.read():
            chrom_counts[rec.chrom] += 1
        self.file.seek(0)
        if not chrom_counts:
            return pd.DataFrame(columns=["chrom", "count"])
        return pd.DataFrame(list(chrom_counts.items()), columns=["chrom", "count"])

    def filter_by_region(self, chrom: str, start: int, end: int) -> Iterator[AlignmentRecord]:
        """
        Фильтрует выравнивания по заданному геномному региону.

        Возвращает итератор по выравниваниям, перекрывающим регион [start, end]
        на указанной хромосоме (включительно, 1-based или 0-based — зависит от данных,
        но сравнение выполняется как числовое перекрытие).

        Args:
            chrom (str): Название хромосомы (например, 'chr1').
            start (int): Начало региона (включительно).
            end (int): Конец региона (включительно).

        Yields:
            AlignmentRecord: Выравнивания, перекрывающие указанный регион.

        Raises:
            ValueError: Если start > end.

        Note:
            Метод сбрасывает файловый указатель в начало после завершения итерации.
        """
        if start > end:
            raise ValueError("start must be <= end")
        self.file.seek(0)
        for rec in self.read():
            if rec.chrom == chrom and rec.start <= end and rec.end >= start:
                yield rec
        self.file.seek(0)