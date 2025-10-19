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
git clone https://github.com/bioinf-rnrmu-stotoshka/bioformats-bioboba.git
cd bioformats-bioboba
```


## Пример использования
```bash
python biodatareader/run_fastq.py sample.fastq
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
