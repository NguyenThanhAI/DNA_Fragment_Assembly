import os
from typing import List, Dict


class Loader(object):
    
    
    def __init__(self, reads: List[str]) -> None:
        
        self.reads = reads
    
    
    @staticmethod
    def load(filename: str):
        
        if not (filename.endswith(".fq") or filename.endswith(".fastq")):
            raise ValueError("Kiểu file không được hỗ trợ, hãy sử dụng file định dạng .fq hoặc .fastq")
        
        if not os.path.isfile(path=filename):
            raise Exception("File {} không tồn tại".format(filename))
        
        if not os.access(path=filename, mode=os.F_OK):
            raise Exception("File {} không đọc được".format(filename))
        
        reads: List[str] = []
        
        for line in open(file=filename, mode="r"):
            line = line.strip()
            if line[:1] == "@":
                name = line[1:]
                is_seq_next: bool = True
            elif is_seq_next:
                if "N" in line:
                    is_seq_next = False
                    continue
                
                reads.append(line)
                is_seq_next = False
                
        return Loader(reads=reads)
    
    
    def __getitem__(self, n: int) -> str:
        
        return self.reads[n]
    
    
    def __len__(self) -> int:
        
        return len(self.reads)
    
    
    def __str__(self) -> List[str]:
        
        return "".join(str(self.reads))
    
    
#reads = Loader.load(filename="data/hemoglobin.fastq")
#print(reads)