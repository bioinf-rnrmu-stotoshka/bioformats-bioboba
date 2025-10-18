from pathlib import Path
from typing import Dict, List
import matplotlib.pyplot as plt
from fastq_reader import FastqReader

def analyze_fastq(file_path: str | Path):
    
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Файл не найден: {file_path}")
    
    if file_path.stat().st_size == 0:
        print("Файл пуст")
        return

    with FastqReader(file_path) as reader:

        sequence_lengths = []
        gc_content_per_read = []
        mean_qualities_per_read = []
        quality_per_position: Dict[int, List[float]] = {}
        
        for record in reader.read():
            seq_len = len(record.sequence)
            sequence_lengths.append(seq_len)

            gc_count = record.sequence.count('G') + record.sequence.count('C')
            gc_percent = (gc_count / seq_len * 100) if seq_len > 0 else 0
            gc_content_per_read.append(gc_percent)

            mean_qual = sum(record.quality) / len(record.quality) if record.quality else 0
            mean_qualities_per_read.append(mean_qual)

            for pos in range(seq_len):
                if pos not in quality_per_position:
                    quality_per_position[pos] = []
                if pos < len(record.quality):
                    quality_per_position[pos].append(record.quality[pos])

        print(f"Прочитано записей: {reader.get_seq_score()}")
        print(f"Средняя длина: {reader.get_mean_seq_length():.1f}")
        print(f"Среднее качество: {reader.get_mean_quality():.1f}")

        if not sequence_lengths:
            print("Нет данных для построения графиков")
            return

        plt.figure(figsize=(10, 6))
        plt.hist(sequence_lengths, bins=30, color='lightblue', edgecolor='black', alpha=0.7)
        plt.title("Sequence Length Distribution")
        plt.xlabel("Sequence Length (bp)")
        plt.ylabel("Frequency")
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.tight_layout()
        plt.show()

        if quality_per_position:
            plt.figure(figsize=(12, 6))
            positions = sorted(quality_per_position.keys())
            mean_qualities = []
            for pos in positions:
                quals = quality_per_position[pos]
                mean_qualities.append(sum(quals) / len(quals))
            
            plt.plot(positions, mean_qualities, color='blue', linewidth=2, label='Mean Quality')
            plt.axhline(y=20, color='red', linestyle='--', alpha=0.7, label='Threshold (Q20)')
            plt.axhline(y=30, color='green', linestyle='--', alpha=0.7, label='Threshold (Q30)')
            
            plt.title("Per Base Sequence Quality")
            plt.xlabel("Position in read (bp)")
            plt.ylabel("Phred Quality Score")
            plt.legend()
            plt.grid(True, linestyle='--', alpha=0.3)
            plt.tight_layout()
            plt.show()

        plt.figure(figsize=(10, 6))
        plt.hist(gc_content_per_read, bins=30, color='lightgreen', edgecolor='black', alpha=0.7)
        plt.axvline(x=50, color='red', linestyle='--', alpha=0.7, label='Expected (50%)')
        
        mean_gc = sum(gc_content_per_read) / len(gc_content_per_read)
        plt.axvline(x=mean_gc, color='blue', linestyle='--', alpha=0.7, 
                   label=f'Mean ({mean_gc:.1f}%)')
        
        plt.title("GC Content Distribution")
        plt.xlabel("GC Content (%)")
        plt.ylabel("Frequency")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(10, 6))
        plt.scatter(gc_content_per_read, mean_qualities_per_read,
                   alpha=0.5, color='purple', s=10)
        plt.xlabel("GC Content (%)")
        plt.ylabel("Mean Read Quality")
        plt.title("Quality vs GC Content")
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.tight_layout()
        plt.show()