"""
Демонстрационная программа для анализа VCF-файлов.

1. Получение заголовка и информации по отдельным группам заголовков (##INFO, ##FILTER и т.д.)
2. Получение количества вариантов.
3. Получение статистики “количество вариантов — регион” (регион = хромосома) с использованием pandas.
4. Получение вариантов в заданном геномном отрезке (аналог bedtools intersect).


Example:
    Запуск без региона:

    .. code-block:: bash

        python run_vcf.py sample.vcf

    Запуск с фильтрацией по региону:

    .. code-block:: bash

        python run_vcf.py sample.vcf chr1 10000 20000

"""

import sys
import argparse
from pathlib import Path
from vcf_reader import VcfReader


def main() -> None:
    """
    Основная функция CLI-утилиты для анализа VCF-файлов.


    - Выводит мета-заголовки VCF-файла (строки, начинающиеся с '##').
    - Группирует и отображает информацию по ключевым секциям заголовка: ##INFO, ##FILTER, ##FORMAT, ##contig.
    - Подсчитывает общее количество вариантов в файле.
    - Формирует и выводит статистику по хромосомам (регионам) с использованием pandas.
    - (Опционально) Фильтрует и отображает варианты в заданном геномном регионе.

    Аргументы командной строки:

    - vcf_file (str): Обязательный путь к VCF-файлу.
    - chrom (str, optional): Название хромосомы для фильтрации (например, 'chr1').
    - start (int, optional): Начало региона (1-based, включительно).
    - end (int, optional): Конец региона (включительно).

    Программа завершается с кодом 1 в следующих случаях:

    - Указанный файл не существует.
    - Задана хромосома, но не указаны обе координаты (start и end).
    - Нарушены ограничения на координаты (start < 1 или start > end).
    - Произошла ошибка при чтении или обработке VCF-файла.

    Вывод:

    - Краткая сводка по заголовкам.
    - Статистика по группам метаинформации.
    - Общее число вариантов.
    - Таблица распределения вариантов по хромосомам.
    - Список вариантов в указанном регионе (максимум 5 первых, если их больше).

    Note:
        Для корректной работы требуется, чтобы класс VcfReader реализовывал
        следующие методы:
            - get_header() → list[str]
            - get_header_group(group: str) → list[str]
            - count_variants() → int
            - stats_by_region() → pandas.DataFrame
            - filter_by_region(chrom: str, start: int, end: int) → Iterator[VariantRecord]

    Raises:
        SystemExit: При ошибках валидации входных данных или обработки файла.
    """
    parser = argparse.ArgumentParser(
        description="Анализ VCF-файлов: заголовки, статистика, фильтрация по региону."
    )
    parser.add_argument(
        "vcf_file",
        type=Path,
        help="Путь к входному VCF-файлу"
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
        help="Начало региона (1-based, включительно)"
    )
    parser.add_argument(
        "end",
        nargs="?",
        type=int,
        default=None,
        help="Конец региона (включительно)"
    )
    args = parser.parse_args()

    # Проверка существования файла
    if not args.vcf_file.is_file():
        print(f"Ошибка: файл не найден — {args.vcf_file}", file=sys.stderr)
        sys.exit(1)

    # Проверка корректности региона
    if args.chrom is not None:
        if args.start is None or args.end is None:
            print("Ошибка: если указана хромосома, должны быть заданы и START, и END.", file=sys.stderr)
            sys.exit(1)
        if args.start > args.end or args.start < 1:
            print("Ошибка: START должен быть ≥ 1 и ≤ END.", file=sys.stderr)
            sys.exit(1)

    try:
        with VcfReader(args.vcf_file) as reader:
            print("=" * 70)
            print("1. ЗАГОЛОВКИ VCF-ФАЙЛА")
            print("=" * 70)
            header = reader.get_header()
            if header:
                print(f"Найдено {len(header)} мета-заголовков (##...)")
                # Показываем первые 5, чтобы не засорять вывод
                for line in header[:5]:
                    print(f"  {line}")
                if len(header) > 5:
                    print(f"  ... и ещё {len(header) - 5} строк")
            else:
                print("Мета-заголовки не найдены.")

            print("\n" + "=" * 70)
            print("2. ИНФОРМАЦИЯ ПО ГРУППАМ ЗАГОЛОВКОВ")
            print("=" * 70)
            groups_to_show = ["INFO", "FILTER", "FORMAT", "contig"]
            for group in groups_to_show:
                entries = reader.get_header_group(group)
                if entries:
                    print(f"\n##{group} — найдено записей: {len(entries)}")
                    for entry in entries[:2]:  # Показываем максимум 2
                        print(f"  {entry}")
                    if len(entries) > 2:
                        print(f"    ... и ещё {len(entries) - 2}")
                else:
                    print(f"\n##{group} — не найдены")

            print("\n" + "=" * 70)
            print("3. КОЛИЧЕСТВО ВАРИАНТОВ")
            print("=" * 70)
            total = reader.count_variants()
            print(f"Общее количество вариантов: {total:,}")

            print("\n" + "=" * 70)
            print("4. СТАТИСТИКА ПО РЕГИОНАМ (ХРОМОСОМАМ)")
            print("=" * 70)
            stats_df = reader.stats_by_region()
            if stats_df.empty:
                print("Нет вариантов для анализа.")
            else:
                print(stats_df.to_string(index=False))
                print(f"\nВсего регионов (хромосом) с вариантами: {len(stats_df)}")

            # === Фильтрация по региону ===
            if args.chrom is not None:
                print(f"\n" + "=" * 70)
                print(f"5. ВАРИАНТЫ В РЕГИОНЕ: {args.chrom}:{args.start}-{args.end}")
                print("=" * 70)
                variants_in_region = list(reader.filter_by_region(args.chrom, args.start, args.end))
                print(f"Найдено вариантов: {len(variants_in_region)}")
                for i, var in enumerate(variants_in_region[:5], 1):  # Показываем первые 5
                    print(f"  {i}. {var.chrom}:{var.pos} {var.ref}>{var.alt}")
                if len(variants_in_region) > 5:
                    print(f"  ... и ещё {len(variants_in_region) - 5}")

    except Exception as e:
        print(f"Ошибка при обработке VCF-файла: {e}", file=sys.stderr)
        sys.exit(1)

    print("\nАнализ завершён.")


if __name__ == "__main__":
    main()