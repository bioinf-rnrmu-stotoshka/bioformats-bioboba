from pathlib import Path
from sam_reader import SamReader

# === 1. Открываем SAM-файл ===
sam_path = Path("team-project/Col0_C1.100k.sam")
with SamReader(sam_path) as reader:
    # Автоматически разбираем заголовок
    reader._parse_header()

    # === 2. Получаем и выводим заголовок ===
    print("\n=== Заголовки SAM-файла ===")
    header = reader.get_header()
    for tag, entries in header.items():
        print(f"{tag}:")
        # Убираем дубликаты и сортируем
        for line in sorted(set(entries)):
            print(f"  {line}")

    # === 3. Получаем количество выравниваний ===
    total = reader.count_alignments()
    print(f"\n=== Общее количество выравниваний: {total}")

    # === 4. Статистика по хромосомам ===
    print("\n=== Статистика по хромосомам ===")
    df_stats = reader.stats_by_chromosome()
    print(df_stats)

    # === 5. Получаем выравнивания в определённом регионе ===
    # Используем реальные имена хромосом из заголовка (@SQ SN:XXX)
    # Например, если в заголовке есть SN:1, используем chrom="1"
    chrom = "1"
    start = 10000
    end = 20000

    print(f"\n=== Выравнивания в регионе {chrom}:{start}-{end} ===")
    for rec in reader.filter_by_region(chrom, start, end):
        print(f"{rec.id}\t{rec.chrom}\t{rec.start}\t{rec.end}\t{rec.cigar}")
