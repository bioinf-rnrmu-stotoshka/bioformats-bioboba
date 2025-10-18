from pathlib import Path
from typing import Iterator
from abstract import SequenceReader
from record import SequenceRecord


class FastaReader(SequenceReader):
    """
    Реализация ридера для чтения FASTA-файлов.

    Поддерживает итеративное чтение записей в формате FASTA и сбор базовой статистики:
    количество последовательностей и средняя длина. Каждая запись преобразуется
    в объект SequenceRecord.

    Attributes:
        filepath (Path): Путь к FASTA-файлу.
        file (file object or None): Открытый файловый дескриптор.
        _seq_count (int): Количество прочитанных последовательностей.
        _total_length (int): Суммарная длина всех прочитанных последовательностей.
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализирует FastaReader с указанным путём к файлу.

        Args:
            filepath (str | Path): Путь к FASTA-файлу.
        """
        super().__init__(filepath)
        self.file = None
        self._seq_count = 0
        self._total_length = 0

    def __enter__(self):
        """
        Поддержка контекстного менеджера (with-блока).

        Открывает FASTA-файл в режиме чтения и возвращает экземпляр FastaReader.

        Returns:
            FastaReader: Текущий экземпляр после открытия файла.
        """
        self.file = open(self.filepath, "r")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Завершение работы контекстного менеджера.

        Автоматически закрывает файл при выходе из with-блока.

        Args:
            exc_type (type or None): Тип исключения, если оно возникло.
            exc_value (Exception or None): Экземпляр исключения.
            traceback (traceback or None): Объект трассировки стека.
        """
        self.close()

    def close(self):
        """
        Закрывает открытый файл, если он существует и не закрыт.

        Устанавливает self.file в None после закрытия.
        """
        if self.file:
            self.file.close()
            self.file = None

    def read(self) -> Iterator[SequenceRecord]:
        """
        Итеративно читает FASTA-файл и возвращает объекты SequenceRecord.

        Метод обрабатывает многострочные последовательности, объединяя строки
        между заголовками. Поддерживает пропуск пустых строк.

        Yields:
            SequenceRecord: Объект, содержащий идентификатор и последовательность.

        Raises:
            ValueError: Если обнаружена недопустимая последовательность (см. _validate_sequence).
            FileNotFoundError: Если файл не найден (при первом обращении к self.file).
        """
        if not self.file:
            self.file = open(self.filepath, "r")

        seq_id, seq_lines = None, []

        for line in self.file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Если уже есть накопленная последовательность — отдать
                if seq_id is not None:
                    record = self._get_sequence(seq_id, "".join(seq_lines))
                    yield record
                    # учёт статистики
                    self._seq_count += 1
                    self._total_length += len(record.sequence)
                seq_id = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)

        # Последняя запись
        if seq_id:
            record = self._get_sequence(seq_id, "".join(seq_lines))
            yield record
            self._seq_count += 1
            self._total_length += len(record.sequence)

    def _get_sequence(self, seq_id: str, seq: str) -> SequenceRecord:
        """
        Создаёт объект SequenceRecord из идентификатора и строки последовательности.

        Приводит последовательность к верхнему регистру и проверяет её валидность.

        Args:
            seq_id (str): Идентификатор последовательности (без символа '>').
            seq (str): Строка последовательности (может быть многострочной до объединения).

        Returns:
            SequenceRecord: Валидный объект последовательности.

        Raises:
            ValueError: Если последовательность содержит недопустимые символы.
        """
        seq = seq.strip().upper()
        if not self._validate_sequence(seq):
            raise ValueError(f"Некорректная последовательность для {seq_id}")
        return SequenceRecord(id=seq_id, sequence=seq)

    def _validate_sequence(self, seq: str) -> bool:
        """
        Проверяет, что последовательность состоит только из допустимых символов.

        Использует расширенный алфавит IUPAC для нуклеотидов, включая символы
        неопределённости и пробелы (дефис для выравниваний).

        Args:
            seq (str): Последовательность для проверки.

        Returns:
            bool: True, если все символы допустимы; иначе False.
        """
        allowed_chars = set("ACGTURYKMSWBDHVN-")  # стандартные нуклеотидные символы IUPAC
        return all(base in allowed_chars for base in seq)

    def get_seq_score(self) -> int:
        """
        Возвращает количество прочитанных последовательностей.

        Note:
            Счётчик обновляется только при вызове метода read() и проходе по итератору.

        Returns:
            int: Число последовательностей, обработанных на данный момент.
        """
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """
        Возвращает среднюю длину прочитанных последовательностей.

        Если ни одна последовательность не была прочитана, возвращает 0.0.

        Returns:
            float: Средняя длина последовательности в нуклеотидах.
        """
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count