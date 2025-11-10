"""
BioDataReader — библиотека для чтения и анализа геномных файлов
(FASTA, FASTQ, SAM, VCF).
"""
# Readers
from .fasta_reader import FastaReader
from .fastq_reader import FastqReader
from .sam_reader import SamReader
from .vcf_reader import VcfReader

# Базовые классы и записи
from .abstract import Reader, GenomicDataReader, SequenceReader
from .record import Record, SequenceRecord, AlignmentRecord, VariantRecord

# Анализ
from .analyze_fastq import analyze_fastq


# Метаданные пакета
__version__ = "1.0.0"
__author__ = "BioBoba team"
__license__ = "MIT"