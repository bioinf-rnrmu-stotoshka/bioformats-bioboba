# test_run.py
from sam_reader import SamReader
from vcf_reader import VcfReader

def test_sam(file_path):
    """Быстрая проверка SAM ридера"""
    print("Testing SAM reader...")
    try:
        with SamReader(file_path) as reader:
            # Базовые методы
            header = reader.get_header()
            count = reader.count_alignments()
            chroms = reader.get_chromosomes()
            
            print(f"Header: {len(header)} lines")
            print(f"Alignments: {count}")
            print(f"Chromosomes: {chroms[:3]}")  # первые 3
            
            # Первая запись - ИСПРАВЛЕНО: используем rec.start вместо rec.pos
            for rec in reader.read():
                print(f"First record: {rec.chrom}:{rec.start}")  # ← ИСПРАВЛЕНО ЗДЕСЬ
                break
                
    except Exception as e:
        print(f"SAM Error: {e}")

def test_vcf(file_path):
    """Быстрая проверка VCF ридера"""
    print("\nTesting VCF reader...")
    try:
        with VcfReader(file_path) as reader:
            # Базовые методы
            header = reader.get_header()
            count = reader.count_variants()
            chroms = reader.get_chromosomes()
            
            print(f"Header: {len(header)} lines")
            print(f"Variants: {count}")
            print(f"Chromosomes: {chroms[:3]}")  # первые 3
            
            # Первая запись
            for rec in reader.read():
                print(f"First record: {rec.chrom}:{rec.pos} {rec.ref}>{rec.alt}")
                break
                
    except Exception as e:
        print(f"VCF Error: {e}")

if __name__ == "__main__":
    # Укажи пути к своим файлам
    sam_file = "Col0_C1.100k.sam"  # поменяй на свой SAM
    vcf_file = "HG00098.vcf"  # поменяй на свой VCF
    
    test_sam(sam_file)
    test_vcf(vcf_file)