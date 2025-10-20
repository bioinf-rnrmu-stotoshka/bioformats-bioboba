# BioDataReader

Библиотека для чтения и анализа геномных файлов в форматах **FASTA**, **FASTQ**, **SAM** и **VCF**.

## Возможности

- Чтение последовательностей из FASTA/FASTQ
- Парсинг SAM-файлов (выравниваний)
- Анализ VCF-файлов (вариантов)
- Валидация данных и статистика
- Готовые CLI-скрипты для быстрого запуска

## Установка

```bash
# Клонирование репозитория
git clone https://github.com/bioinf-rnrmu-stotoshka/bioformats-bioboba.git
cd bioformats-bioboba
```

## Пример использования
```bash
python biodatareader/run_<формат>.py <входной_файл> [дополнительные_параметры]
```
### Класс FastaReader
```bash
python biodatareader/run_fasta.py GCA_000006945.2_ASM694v2_genomic.fna
```
**Пример вывода**
```text
Количество последовательностей: 2
Средняя длина последовательности: 2475691.50
```
### Класс FastqReader
```bash
python biodatareader/run_fastq.py SRR3280893_1.fastq
```
**Пример вывода**

![1q](https://github.com/user-attachments/assets/87ba1d3a-53b9-4bcb-a218-72cba5191b46)
![2q](https://github.com/user-attachments/assets/67da3057-7914-4215-83c3-caccc59ac130)
![3q](https://github.com/user-attachments/assets/2c13dbfa-751e-471e-99cf-673d9840e5e7)

### Класс SamReader
```bash
python biodatareader/run_sam.py Col0_C1.100k.sam
```
**Пример вывода**
```text
=== Заголовки SAM-файла ===
@SQ:
  SN:1  LN:30427671
  SN:2  LN:19698289
  SN:3  LN:23459830
  SN:4  LN:18585056
  SN:5  LN:26975502
  SN:C  LN:154478
  SN:M  LN:366924

=== Общее количество выравниваний: 34,298

=== Статистика по хромосомам ===
chrom  count
    1  34298
```
### Класс VcfReader
```bash
python biodatareader/run_vcf.py HG00098.vcf
```
**Пример вывода**
```text
======================================================================
1. ЗАГОЛОВКИ VCF-ФАЙЛА
======================================================================
Найдено 22 мета-заголовков (##...)
  ##fileformat=VCFv4.0
  ##filedat=20101112
  ##datarelease=20100804
  ##samples=629
  ##description="Where BI calls are present, genotypes and alleles are from BI.  In there absence, UM genotypes are used.  If neither are available, no genotype information is present and the alleles are from the NCBI calls."
  ... и ещё 17 строк

======================================================================
2. ИНФОРМАЦИЯ ПО ГРУППАМ ЗАГОЛОВКОВ
======================================================================

##INFO — найдено записей: 9
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
  ##INFO=<ID=CB,Number=.,Type=String,Description="List of centres that called, UM (University of Michigan), BI (Broad Institute), BC (Boston College), NCBI">
    ... и ещё 7

##FILTER — не найдены

##FORMAT — найдено записей: 7
  ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
    ... и ещё 5

##contig — не найдены

======================================================================
3. КОЛИЧЕСТВО ВАРИАНТОВ
======================================================================
Общее количество вариантов: 46,065

======================================================================
4. СТАТИСТИКА ПО РЕГИОНАМ (ХРОМОСОМАМ)
======================================================================
region  variant_count
    22          46065

Всего регионов (хромосом) с вариантами: 1

Анализ завершён.
```

## Диаграмма классов
```mermaid
classDiagram
    %% Abstract Classes
    class Reader {
        <<abstract>>
        -filepath: Path
        -file: file object or None
        +__init__(filepath: str | Path)
        +read()* Iterator[Record]
        +close()
        +__enter__() Reader
        +__exit__(exc_type, exc_val, exc_tb)
    }

    class SequenceReader {
        <<abstract>>
        +__init__(filepath: str | Path)
        +read()* Iterator[SequenceRecord]
    }

    class GenomicDataReader {
        <<abstract>>
        -_header_parsed: bool
        +__init__(filepath: str | Path)
        +__enter__() GenomicDataReader
        +_parse_header()*
        +read()* Iterator[Record]
        +close()
    }

    %% Record Classes
    class Record {
        <<abstract>>
        -id: str
        +__init__(id: str)
        +__repr__() str
    }

    class SequenceRecord {
        +sequence: str
        +quality: list[int] or None
        +__init__(id: str, sequence: str, quality: list[int] or None)
    }

    class AlignmentRecord {
        +chrom: str
        +start: int
        +end: int
        +cigar: str
        +mapq: int
        +flag: int
        +__init__(id: str, chrom: str, start: int, cigar: str, mapq: int)
        +__repr__() str
    }

    class VariantRecord {
        +chrom: str
        +pos: int
        +ref: str
        +alt: str
        +info: dict
        +filter: str
        +__init__(chrom: str, pos: int, ref: str, alt: str, info: dict)
        +__repr__() str
    }

    %% Concrete Reader Implementations
    class FastaReader {
        -_seq_count: int
        -_total_length: int
        +__init__(filepath: str | Path)
        +__enter__() FastaReader
        +__exit__(exc_type, exc_value, traceback)
        +close()
        +read() Iterator[SequenceRecord]
        -_get_sequence(seq_id: str, seq: str) SequenceRecord
        -_validate_sequence(seq: str) bool
        +get_seq_score() int
        +get_mean_seq_length() float
    }

    class FastqReader {
        +__init__(filepath: str | Path)
        +__enter__() FastqReader
        +__exit__(exc_type, exc_value, traceback)
        +close()
        +read() Iterator[SequenceRecord]
        +_parse_quality(quality_str: str) list[int]$
    }

    class SamReader {
        -header: Dict[str, List[str]]
        +__init__(filepath: str | Path)
        +_ensure_header_parsed()
        +_parse_header()
        -_parse_header_line(line: str)
        +get_header() Dict[str, List[str]]
        +get_header_group(group_tag: str) List[str]
        +read() Iterator[AlignmentRecord]
        +_calc_aligned_length(cigar: str) int$
        +count_alignments() int
        +stats_by_chromosome() pd.DataFrame
        +filter_by_region(chrom: str, start: int, end: int) Iterator[AlignmentRecord]
    }

    class VcfReader {
        -header_lines: List[str]
        -column_header: List[str]
        +__init__(filepath: str | Path)
        +_parse_header()
        +get_header() List[str]
        +get_header_group(key: str) List[str]
        +read() Iterator[VariantRecord]
        +_parse_info(info_str: str) Dict[str, str]$
        +count_variants() int
        +stats_by_region() pd.DataFrame
        +filter_by_region(chrom: str, start: int, end: int) Iterator[VariantRecord]
    }

    %% Inheritance Relationships
    Reader <|-- SequenceReader
    Reader <|-- GenomicDataReader
    SequenceReader <|-- FastaReader
    SequenceReader <|-- FastqReader
    GenomicDataReader <|-- SamReader
    GenomicDataReader <|-- VcfReader

    Record <|-- SequenceRecord
    Record <|-- AlignmentRecord
    Record <|-- VariantRecord

    %% Association Relationships
    Reader --> Record : reads
    SequenceReader --> SequenceRecord : reads
    GenomicDataReader --> Record : reads
    SamReader --> AlignmentRecord : reads
    VcfReader --> VariantRecord : reads
```


## Документация
Полная документация с описанием классов и методов доступна в папке docs/ .
Чтобы собрать локально:

```bash
cd docs
make html # для Windows - ./make html
```


## Вклад в проект
См. CONTRIBUTING.md

## Лицензия

Этот проект распространяется под лицензией MIT.
