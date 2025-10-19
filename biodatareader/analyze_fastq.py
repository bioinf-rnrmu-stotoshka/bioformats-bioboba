import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List
from fastq_reader import FastqReader


def analyze_fastq(file_path: str | Path) -> None:
    """
    Анализирует FASTQ-файл и визуализирует ключевые метрики качества последовательностей.

    Функция выполняет следующие действия:
        - Проверяет существование и непустоту файла.
        - Собирает статистику по:
            * длине последовательностей,
            * среднему качеству по позициям (Phred score),
            * содержанию нуклеотидов (A, T, G, C) по позициям.
        - Выводит текстовую сводку в консоль.
        - Строит три графика:
            1. Распределение длин последовательностей.
            2. Среднее качество по каждой позиции в риде.
            3. Процентное содержание каждого нуклеотида по позициям.

    Args:
        file_path (str | Path): Путь к FASTQ-файлу для анализа.

    Raises:
        FileNotFoundError: Если указанный файл не существует.
        RuntimeError: Если возникает ошибка при чтении файла (через FastqReader).

    Returns:
        None
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    if file_path.stat().st_size == 0:
        print("File is empty")
        return

    # Сбор статистики
    sequence_lengths: List[int] = []
    quality_per_position: Dict[int, List[int]] = {}
    base_content_per_position: Dict[int, Dict[str, int]] = {}  # pos -> {A: N, T: N, ...}

    total_sequences = 0
    total_length = 0

    with FastqReader(file_path) as reader:
        for record in reader.read():
            seq = record.sequence
            qual = record.quality
            seq_len = len(seq)

            total_sequences += 1
            total_length += seq_len
            sequence_lengths.append(seq_len)

            # Сбор качества по позициям
            for i, q in enumerate(qual):
                if i not in quality_per_position:
                    quality_per_position[i] = []
                quality_per_position[i].append(q)

            # Сбор содержания оснований по позициям
            for i, base in enumerate(seq):
                if i not in base_content_per_position:
                    base_content_per_position[i] = {"A": 0, "T": 0, "G": 0, "C": 0}
                if base in "ATGC":
                    base_content_per_position[i][base] += 1

    if total_sequences == 0:
        print("No valid sequences found")
        return

    mean_length = total_length / total_sequences
    print(f"Total sequences: {total_sequences}")
    print(f"Mean sequence length: {mean_length:.1f} bp")

    # === 1. Sequence Length Distribution ===
    plt.figure(figsize=(10, 6))
    plt.hist(sequence_lengths, bins=min(50, len(set(sequence_lengths))),
             color='lightblue', edgecolor='black', alpha=0.7)
    plt.title("Sequence Length Distribution")
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()

    # === 2. Per Base Sequence Quality ===
    if quality_per_position:
        positions = sorted(quality_per_position.keys())
        mean_qualities = [sum(quality_per_position[p]) / len(quality_per_position[p]) for p in positions]

        plt.figure(figsize=(12, 6))
        plt.plot(positions, mean_qualities, color='blue', linewidth=2, label='Mean Quality')
        plt.axhline(y=20, color='red', linestyle='--', alpha=0.7, label='Q20')
        plt.axhline(y=30, color='green', linestyle='--', alpha=0.7, label='Q30')
        plt.title("Per Base Sequence Quality")
        plt.xlabel("Position in Read (bp)")
        plt.ylabel("Phred Quality Score")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.tight_layout()
        plt.show()

    # === 3. Per Base Sequence Content ===
    if base_content_per_position:
        positions = sorted(base_content_per_position.keys())
        a_counts = []
        t_counts = []
        g_counts = []
        c_counts = []

        for p in positions:
            total_at_pos = sum(base_content_per_position[p].values())
            if total_at_pos == 0:
                a = t = g = c = 0
            else:
                a = base_content_per_position[p]["A"] / total_at_pos * 100
                t = base_content_per_position[p]["T"] / total_at_pos * 100
                g = base_content_per_position[p]["G"] / total_at_pos * 100
                c = base_content_per_position[p]["C"] / total_at_pos * 100
            a_counts.append(a)
            t_counts.append(t)
            g_counts.append(g)
            c_counts.append(c)

        plt.figure(figsize=(12, 6))
        plt.plot(positions, a_counts, label='A', color='green')
        plt.plot(positions, t_counts, label='T', color='red')
        plt.plot(positions, g_counts, label='G', color='orange')
        plt.plot(positions, c_counts, label='C', color='blue')
        plt.axhline(y=25, color='black', linestyle='--', alpha=0.5, linewidth=0.8)
        plt.title("Per Base Sequence Content")
        plt.xlabel("Position in Read (bp)")
        plt.ylabel("Percentage (%)")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.tight_layout()
        plt.show()

    print("Analysis completed.")