from assembly import Assembler


if __name__ == "__main__":   
    assembly: Assembler = Assembler(filename="data/hemoglobin.fastq", k=11)
    assembly.make_superpath()
    print(assembly.find_eulerian_path())