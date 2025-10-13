from abc import ABC, abstractmethod
from typing import Iterator


class Reader(ABC):

    def __init__(self, filename: str):
        self.filename = filename


    @abstractmethod
    def read(self) -> Iterator["Record"]:
        """
        Абстрактный метод
        Должен возвращать итератор по объектам типа Record (последовательности, выравнивания или варианты)
        """
        pass


    @abstractmethod
    def _parse_line(self, line: str) -> "Record":
        """
        Абстрактный метод
        Определяет, как парсить строку файла в объект Record
        """
        pass


    def close(self):
        pass


    def __enter__(self):
        return self


    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
