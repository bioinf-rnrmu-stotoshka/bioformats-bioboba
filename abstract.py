from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Tuple
from record import Record, SequenceRecord
from pathlib import Path


class Reader(ABC):
    """
    Абстрактный класс, который показывает логику работу ридеров
    """

    def __init__(self, filepath: str | Path):
        self.filepath = Path(filepath)
        self.file = None

    @abstractmethod
    def read(self) -> Iterator[Record]:
        """
        Абстрактный метод
        Должен возвращать итератор по объектам типа Record (последовательности, выравнивания или варианты)
        """
        pass

    def close(self):
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    def __enter__(self):
        self.file = open(self.filepath, "r", encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class SequenceReader(Reader):
    """
    Абстрактный класс, который показывает логику работы ридеров FASTA и FASTQ
    """

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)

    @abstractmethod
    def read(self) -> Iterator[SequenceRecord]:
        """
        Абстрактный метод
        Должен возвращать итератор объектов SequenceRecord
        """
        pass

    def __enter__(self):
        """
        Открывает файл при входе в блок with
        """
        self.file = open(self.filepath, "r", encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Автоматически закрывает файл при выходе из блока 'with'.
        """
        self.close()

    def close(self):
        """
        Закрывает открытый файл (если он был открыт).
        """
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    @abstractmethod
    def get_sequence(self, seq: str, id: str) -> SequenceRecord:
        """
        Создаёт объект SequenceRecord из идентификатора и последовательности.
        Может добавлять качество (в FASTQ) или пропускать его (в FASTA).
        """
        pass

    @abstractmethod
    def validate_sequence(self, seq: str) -> bool:
        """
        Проверяет корректность последовательности
        """
        pass


class GenomicDataReader(Reader):
    """
    Абстрактный класс для ридеров геномных данных (SAM, VCF)
    """

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self._header_parsed = False
        self._chromosomes = []

    def __enter__(self):
        """
        Открывает файл и парсит заголовок при входе в контекст
        """
        try:
            self.file = open(self.filepath, "r", encoding="utf-8")
            self._parse_header()
            return self
        except Exception as e:
            # если ошибка — аккуратно закрываем файл, чтобы не остался висеть
            if self.file and not self.file.closed:
                self.file.close()
            raise RuntimeError(f"Ошибка при открытии или парсинге файла {self.filepath}: {e}")

    @abstractmethod
    def _parse_header(self):
        """
        Абстрактный метод для парсинга заголовка файла
        Должен быть реализован в дочерних классах
        """
        pass

    @abstractmethod
    def read(self) -> Iterator[Record]:
        """
        Чтение записей с пропуском заголовочных строк
        """
        pass

    @abstractmethod
    def get_chromosomes(self) -> List[str]:
        """Возвращает список хромосом"""
        pass

    @abstractmethod
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Проверяет корректность координат"""
        pass

    @abstractmethod
    def get_records_in_region(self, chrom: str, start: int, end: int) -> Iterator[Record]:
        """
        Возвращает записи, попадающие в указанный регион
        """
        pass

    @abstractmethod
    def filter_records(self, **filters) -> Iterator[Record]:
        """
        Фильтрация записей по различным критериям
        filters: словарь с критериями фильтрации
        """
        pass

    def get_chromosome_length(self, chrom: str) -> int:
        """Получить длину хромосомы (если доступно из заголовка)"""
        # Базовая реализация, может быть переопределена в дочерних классах
        pass

    def get_statistics(self) -> Dict[str, any]:
        """
        Базовая статистика по файлу
        Должна быть расширена в дочерних классах
        """
        return {
            "file_path": str(self.filepath),
            "file_size": self.filepath.stat().st_size if self.filepath.exists() else 0,
            "chromosomes": self.get_chromosomes(),
            "chromosome_count": len(self.get_chromosomes()),
        }

    def close(self):
        """Закрывает файл"""
        super().close()
        self._header_parsed = False
