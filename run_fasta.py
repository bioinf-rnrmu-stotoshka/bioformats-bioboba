import sys
from pathlib import Path
from fasta_reader import FastaReader


def main() -> None:
    """
    Основная функция демонстрационного скрипта для анализа FASTA-файлов.

    Скрипт принимает ровно один аргумент командной строки — путь к FASTA-файлу.
    Выполняет следующие действия:
        - Проверяет корректность аргумента и существование файла.
        - Читает все последовательности из файла с помощью FastaReader.
        - Выводит в консоль количество последовательностей и их среднюю длину.

    В случае ошибки (некорректный ввод, отсутствие файла, недопустимые символы
    в последовательности и т.д.) выводит сообщение об ошибке в stderr и завершает
    выполнение с кодом 1.

    Usage:
        python run_fasta.py <путь_к_файлу.fasta>

    Raises:
        SystemExit: При неверном количестве аргументов, отсутствии файла
                    или ошибках во время чтения/валидации.
    """
    if len(sys.argv) != 2:
        print("Использование: python run_fasta.py <путь_к_файлу.fasta>", file=sys.stderr)
        sys.exit(1)

    fasta_path = Path(sys.argv[1])

    if not fasta_path.is_file():
        print(f"Ошибка: файл '{fasta_path}' не существует или не является файлом.", file=sys.stderr)
        sys.exit(1)

    try:
        with FastaReader(fasta_path) as reader:
            # Прочитать все последовательности (чтобы собрать статистику)
            list(reader.read())  # вызов read() заполняет внутренние счётчики

            count = reader.get_seq_score()
            mean_len = reader.get_mean_seq_length()

            print(f"Количество последовательностей: {count}")
            print(f"Средняя длина последовательности: {mean_len:.2f}")

    except ValueError as e:
        print(f"Ошибка валидации последовательности: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Ошибка при чтении файла: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()