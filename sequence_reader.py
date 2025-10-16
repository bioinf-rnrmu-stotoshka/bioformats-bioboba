from pathlib import Path
from typing import Iterator
from abstract import SequenceReader
from record import SequenceRecord

class FastaReader(SequenceReader):
    """Ридер для чтения FASTA-файлов."""

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.file = None
        self._seq_count = 0
        self._total_length = 0

    def __enter__(self):
        self.file = open(self.filepath, "r")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

    def read(self) -> Iterator[SequenceRecord]:
        """Читает FASTA построчно, возвращая объекты SequenceRecord."""
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
                    record = self.get_sequence(seq_id, "".join(seq_lines))
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
            record = self.get_sequence(seq_id, "".join(seq_lines))
            yield record
            self._seq_count += 1
            self._total_length += len(record.sequence)

    def get_sequence(self, seq_id: str, seq: str) -> SequenceRecord:
        """Создаёт объект SequenceRecord, проверяя валидность последовательности."""
        seq = seq.strip().upper()
        if not self.validate_sequence(seq):
            raise ValueError(f"Некорректная последовательность для {seq_id}")
        return SequenceRecord(id=seq_id, sequence=seq)

    def validate_sequence(self, seq: str) -> bool:
        """Проверяет, что последовательность содержит только корректные символы ДНК/белков."""
        allowed_chars = set("ACGTURYKMSWBDHVN-")  # стандартные нуклеотидные символы IUPAC
        return all(base in allowed_chars for base in seq)

    def get_seq_score(self) -> int:
        """Возвращает количество считанных последовательностей."""
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """Возвращает среднюю длину последовательности."""
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count
