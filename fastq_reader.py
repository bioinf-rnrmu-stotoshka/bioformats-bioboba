from pathlib import Path
from typing import Iterator
import gzip
from abstract import SequenceReader
from record import SequenceRecord

class FastqReader(SequenceReader):

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.file = None
        self._seq_count = 0
        self._total_length = 0
        self._total_quality = 0.0

    def __enter__(self):
        if str(self.filepath).endswith('.gz'):
            self.file = gzip.open(self.filepath, "rt")
        else:
            self.file = open(self.filepath, "r")
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
                self.file = gzip.open(self.filepath, "rt")
            else:
                self.file = open(self.filepath, "r")

        while True:
            header = self.file.readline()
            if not header:
                break
                
            sequence = self.file.readline()
            plus_line = self.file.readline()
            quality = self.file.readline()

            if not (header and sequence and plus_line and quality):
                break

            if not header.startswith("@"):
                raise ValueError(f"Неверный формат FASTQ в строке: {header.strip()}")
            
            if not plus_line.startswith("+"):
                raise ValueError(f"Ожидалась строка '+', получено: {plus_line.strip()}")

            seq_id = header[1:].split()[0].strip() if header else "unknown"
            
            seq_clean = "".join(sequence.split()).upper()
            qual_clean = "".join(quality.split())
            
            if len(seq_clean) != len(qual_clean):
                raise ValueError(f"Длины последовательности и качества не совпадают для {seq_id}")
            
            if not seq_clean:
                raise ValueError(f"Пустая последовательность для {seq_id}")

            record = SequenceRecord(id=seq_id, sequence=seq_clean, quality=self._parse_quality(qual_clean))
            
            self._seq_count += 1
            self._total_length += len(seq_clean)
            self._total_quality += self._calculate_mean_quality(record.quality)
            
            yield record

    def _parse_quality(self, quality_str: str) -> list[int]:
        return [ord(ch) - 33 for ch in quality_str] if quality_str else []

    def _calculate_mean_quality(self, quality: list[int]) -> float:
        if not quality:
            return 0.0
        return sum(quality) / len(quality)

    def _validate_sequence(self, seq: str) -> bool:
        allowed_chars = set("ACGTURYKMSWBDHVN-")
        return all(base in allowed_chars for base in seq)

    def get_seq_score(self) -> int:
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count

    def get_mean_quality(self) -> float:
        if self._seq_count == 0:
            return 0.0
        return self._total_quality / self._seq_count