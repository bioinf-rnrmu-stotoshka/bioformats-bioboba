"""
Анализ SAM-файла с поддержкой командной строки.

Example:
    Запуск без региона:

    .. code-block:: bash

        python run_sam.py sample.sam

    Запуск с фильтрацией по региону:

    .. code-block:: bash

        python run_sam.py sample.sam chr1 10000 20000
"""

import sys
from pathlib import Path
import argparse
from sam_reader import SamReader


def main() -> None:
    """
    Основная функция CLI-утилиты для анализа SAM-файлов.

    Поддерживает:
        - Чтение и отображение заголовков SAM-файла.
        - Подсчёт общего числа выравниваний.
        - Сбор статистики по хромосомам (количество выравниваний на каждую хромосому).
        - Опциональную фильтрацию выравниваний по геномному региону (хромосома, начало, конец).

    Аргументы командной строки:
        sam_file (str): Обязательный путь к SAM-файлу.
        chrom (str, optional): Название хромосомы для фильтрации.
        start (int, optional): Начало региона (1-based).
        end (int, optional): Конец региона.

    Программа завершается с кодом 1 в следующих случаях:
        - Указанный файл не существует.
        - Задана хромосома, но не указаны обе координаты (start и end).
        - Начало региона больше конца.
        - Произошла ошибка при чтении или обработке SAM-файла.

    Вывод:
        - Заголовки SAM (если есть).
        - Общее количество выравниваний.
        - Таблица статистики по хромосомам.
        - Список выравниваний в заданном регионе (если регион указан).

    Raises:
        SystemExit: При ошибках валидации входных данных или обработки файла.
    """
    parser = argparse.ArgumentParser(
        description="Анализ SAM-файла: заголовки, статистика, фильтрация по региону."
    )
    parser.add_argument(
        "sam_file",
        type=Path,
        help="Путь к SAM-файлу"
    )
    parser.add_argument(
        "chrom",
        nargs="?",
        default=None,
        help="Хромосома для фильтрации (например: '1', 'chr1')"
    )
    parser.add_argument(
        "start",
        nargs="?",
        type=int,
        default=None,
        help="Начало региона (1-based)"
    )
    parser.add_argument(
        "end",
        nargs="?",
        type=int,
        default=None,
        help="Конец региона"
    )
    args = parser.parse_args()

    # Проверка существования файла
    if not args.sam_file.is_file():
        print(f"Ошибка: файл не найден — {args.sam_file}", file=sys.stderr)
        sys.exit(1)

    # Проверка корректности региона
    if args.chrom is not None:
        if args.start is None or args.end is None:
            print("Ошибка: если задана хромосома, должны быть указаны и START, и END.", file=sys.stderr)
            sys.exit(1)
        if args.start > args.end:
            print("Ошибка: START не может быть больше END.", file=sys.stderr)
            sys.exit(1)

    # === Основной анализ ===
    try:
        with SamReader(args.sam_file) as reader:
            # === 1. Заголовки ===
            print("\n=== Заголовки SAM-файла ===")
            header = reader.get_header()
            if header:
                for tag, entries in header.items():
                    print(f"{tag}:")
                    # Убираем дубликаты и сортируем
                    for line in sorted(set(entries)):
                        print(f"  {line}")
            else:
                print("Заголовки отсутствуют.")

            # === 2. Количество выравниваний ===
            total = reader.count_alignments()
            print(f"\n=== Общее количество выравниваний: {total:,}")

            # === 3. Статистика по хромосомам ===
            print("\n=== Статистика по хромосомам ===")
            df_stats = reader.stats_by_chromosome()
            if df_stats.empty:
                print("Нет выравниваний для анализа.")
            else:
                print(df_stats.to_string(index=False))

            # === 4. Фильтрация по региону (если задан) ===
            if args.chrom is not None:
                print(f"\n=== Выравнивания в регионе {args.chrom}:{args.start}-{args.end} ===")
                found = False
                for rec in reader.filter_by_region(args.chrom, args.start, args.end):
                    print(f"{rec.id}\t{rec.chrom}\t{rec.start}\t{rec.end}\t{rec.cigar}")
                    found = True
                if not found:
                    print("Нет выравниваний в указанном регионе.")

    except Exception as e:
        print(f"Ошибка при обработке файла: {e}", file=sys.stderr)
        sys.exit(1)

    print()


if __name__ == "__main__":
    main()