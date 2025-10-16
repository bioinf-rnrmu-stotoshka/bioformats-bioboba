from abc import ABC, abstractmethod
from typing import Iterator
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
        pass


    def __enter__(self):
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