import sys
from pathlib import Path
from analyze_fastq import analyze_fastq


def main() -> None:
    """
    Основная функция демонстрационного скрипта для анализа FASTQ-файлов.

    Скрипт принимает ровно один аргумент командной строки — путь к FASTQ-файлу
    (возможно, сжатому в формате .gz). Выполняет следующие действия:
    - Проверяет, что передан ровно один аргумент.
    - Проверяет существование указанного файла.
    - Запускает функцию analyze_fastq для визуального и статистического анализа.

    В случае ошибки (неверное число аргументов или отсутствие файла)
    выводит сообщение об ошибке и завершает выполнение с кодом 1.

    Example:
        Запуск из командной строки:

        .. code-block:: bash

            python run_fastq.py example.fastq

    Raises:
        SystemExit: При неверном количестве аргументов или если файл не найден.
    """
    if len(sys.argv) != 2:
        print("Использование: python run_fastq.py <путь_к_fastq_файлу>")
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