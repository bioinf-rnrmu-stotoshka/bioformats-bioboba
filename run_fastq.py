import sys
from pathlib import Path
from analyze_fastq import analyze_fastq


def main():
    if len(sys.argv) != 2:
        print("Использование: python demo.py <путь_к_fastq_файлу>")
        sys.exit(1)

    file_path = Path(sys.argv[1])

    if not file_path.exists():
        print(f"Ошибка: файл не найден — {file_path}")
        sys.exit(1)

    print(f"Анализ файла: {file_path}")
    print("-" * 40)
    analyze_fastq(file_path)


if __name__ == "__main__":
    main()