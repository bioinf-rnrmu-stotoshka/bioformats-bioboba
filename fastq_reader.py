from pathlib import Path
from typing import Iterator
import gzip
from abstract import SequenceReader
from record import SequenceRecord


class FastqReader(SequenceReader):
    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.file = None

    def __enter__(self):
        if str(self.filepath).endswith('.gz'):
            self.file = gzip.open(self.filepath, "rt", encoding="ascii")
        else:
            self.file = open(self.filepath, "r", encoding="ascii")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

    def read(self) -> Iterator[SequenceRecord]:
        if not self.file:
            if str(self.filepath).endswith('.gz'):
                self.file = gzip.open(self.filepath, "rt", encoding="ascii")
            else:
                self.file = open(self.filepath, "r", encoding="ascii")

        while True:
            header = self.file.readline()
            if not header:
                break

            sequence = self.file.readline().rstrip('\n')
            plus_line = self.file.readline()
            quality = self.file.readline().rstrip('\n')

            if not (header and sequence and plus_line and quality):
                break

            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ: expected '@', got {header.strip()!r}")
            if not plus_line.startswith("+"):
                raise ValueError(f"Invalid FASTQ: expected '+', got {plus_line.strip()!r}")

            seq_id = header[1:].split(maxsplit=1)[0] if len(header) > 1 else "unknown"

            seq_clean = sequence.upper()
            qual_clean = quality

            if len(seq_clean) != len(qual_clean):
                raise ValueError(f"Sequence and quality length mismatch for {seq_id}")

            if not seq_clean:
                raise ValueError(f"Empty sequence for {seq_id}")

            # Валидация качества (опционально)
            quality_scores = self._parse_quality(qual_clean)
            # Можно добавить предупреждение, если Q < 0 или > 60

            record = SequenceRecord(id=seq_id, sequence=seq_clean, quality=quality_scores)
            yield record

    @staticmethod
    def _parse_quality(quality_str: str) -> list[int]:
        # Phred+33: ASCII 33 = Q0, 126 = Q93
        return [ord(ch) - 33 for ch in quality_str]