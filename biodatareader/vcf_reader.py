import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
from collections import Counter
from abstract import GenomicDataReader
from record import VariantRecord


class VcfReader(GenomicDataReader):
    """
    Реализация ридера для чтения VCF-файлов (Variant Call Format).

    Поддерживает:
        - Парсинг мета-заголовков (##...) и колонок данных (#CHROM...).
        - Итеративное чтение вариантов как объектов VariantRecord.
        - Подсчёт общего числа вариантов.
        - Сбор статистики по хромосомам (регион = хромосома).
        - Фильтрацию вариантов по геномному региону (аналог bedtools intersect).

    Attributes:
        filepath (Path): Путь к VCF-файлу.
        file (file object or None): Открытый файловый дескриптор.
        _header_parsed (bool): Флаг, указывающий, был ли уже распарсен заголовок.
        header_lines (List[str]): Список строк мета-заголовков (##...).
        column_header (List[str]): Список названий колонок из строки #CHROM.
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализирует VcfReader с указанным путём к файлу.

        Args:
            filepath (str | Path): Путь к VCF-файлу.
        """
        super().__init__(filepath)
        self.header_lines: List[str] = []
        self.column_header: List[str] = []

    def _parse_header(self) -> None:
        """
        Парсит заголовок VCF-файла.

        Читает файл с начала, извлекает все строки, начинающиеся с '##',
        а также строку с колонками '#CHROM'. Вызывается автоматически
        при первом обращении к данным или при входе в контекстный менеджер.

        Raises:
            ValueError: Если строка '#CHROM' не найдена (некорректный VCF).
        """
        if self._header_parsed:
            return

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
        """
        Возвращает все мета-заголовки VCF-файла (строки, начинающиеся с '##').

        Гарантирует, что заголовок был распарсен перед возвратом.

        Returns:
            List[str]: Список строк мета-заголовков.
        """
        if not self._header_parsed:
            self._parse_header()
        return self.header_lines

    def get_header_group(self, key: str) -> List[str]:
        """
        Возвращает список мета-заголовков, относящихся к указанной группе.

        Например, `get_header_group("INFO")` вернёт все строки вида `##INFO=...`.

        Args:
            key (str): Ключ группы (например, "INFO", "FILTER", "FORMAT").

        Returns:
            List[str]: Список строк заголовка, начинающихся с `##{key}=`.
        """
        if not self._header_parsed:
            self._parse_header()
        prefix = f"##{key}="
        return [line for line in self.header_lines if line.startswith(prefix)]

    def read(self) -> Iterator[VariantRecord]:
        """
        Итеративно читает варианты из VCF-файла и возвращает объекты VariantRecord.

        Пропускает все строки заголовка (начинающиеся с '#'). Парсит обязательные
        поля: CHROM, POS, REF, ALT, а также опциональные FILTER и INFO.
        Некорректные строки пропускаются без ошибки.

        Yields:
            VariantRecord: Объект варианта с атрибутами chrom, pos, ref, alt, info.
            Дополнительно устанавливается атрибут `filter` (из колонки FILTER).

        Note:
            Метод автоматически парсит заголовок при первом вызове.
        """
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
        """
        Парсит поле INFO VCF-файла в словарь.

        Поле INFO состоит из пар ключ=значение, разделённых ';'.
        Если элемент не содержит '=', он интерпретируется как флаг (значение "True").

        Args:
            info_str (str): Строка из колонки INFO (например, "DP=30;AF=0.5;DB").

        Returns:
            Dict[str, str]: Словарь, где ключи — идентификаторы аннотаций,
            значения — строки (числа остаются в виде строк).

        Example:
            >>> VcfReader._parse_info("DP=30;AF=0.5;DB")
            {'DP': '30', 'AF': '0.5', 'DB': 'True'}
        """
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
        """
        Подсчитывает общее количество вариантов в VCF-файле.

        Считает все непустые строки, не начинающиеся с '#'.

        Returns:
            int: Общее число вариантов.

        Note:
            Метод сбрасывает файловый указатель в начало после подсчёта.
        """
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
        Собирает статистику количества вариантов по хромосомам.

        Возвращает pandas DataFrame с двумя колонками: 'region' (хромосома)
        и 'variant_count' (число вариантов). Результат отсортирован по региону.

        Returns:
            pd.DataFrame: Таблица с количеством вариантов на каждую хромосому.
            Если вариантов нет, возвращает пустой DataFrame с указанными колонками.

        Note:
            Метод сбрасывает файловый указатель в начало после сбора статистики.
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
        Фильтрует варианты по заданному геномному региону.

        Возвращает итератор по вариантам, позиция которых (POS) попадает
        в диапазон [start, end] (включительно, 1-based, как в спецификации VCF).

        Args:
            chrom (str): Название хромосомы (например, 'chr1').
            start (int): Начало региона (включительно, ≥1).
            end (int): Конец региона (включительно).

        Yields:
            VariantRecord: Варианты, удовлетворяющие условию фильтрации.

        Raises:
            ValueError: Если start > end.

        Note:
            Метод сбрасывает файловый указатель в начало после завершения итерации.
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