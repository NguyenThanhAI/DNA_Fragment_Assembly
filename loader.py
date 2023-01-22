import os
from typing import List, Dict


class Loader(object):
    """
    Lưu thông tin về gen
    """
    
    def __init__(self, reads: List[str]) -> None:
        """
        Sử dụng phương thức tĩnh Loader.load("filename.fq") để khởi tạo đối tượng lưu các reads
        """
        self.reads = reads
    
    
    @staticmethod
    def load(filename: str):
        
        # Kiểm tra loại file
        if not (filename.endswith(".fq") or filename.endswith(".fastq")):
            raise ValueError("Kiểu file không được hỗ trợ, hãy sử dụng file định dạng .fq hoặc .fastq")
        
        # Kiểm tra file có tồn tại không
        if not os.path.isfile(path=filename):
            raise Exception("File {} không tồn tại".format(filename))
        
        # Kiểm tra file có đọc được không
        if not os.access(path=filename, mode=os.F_OK):
            raise Exception("File {} không đọc được".format(filename))
        
        # Tạo một list các string để lưu các read
        reads: List[str] = []
        
        # Đọc từng dòng của file
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
        """
        Đọc một read tại vị trí n
        """
        
        return self.reads[n]
    
    
    def __len__(self) -> int:
        """
        Số các reads đọc được
        """
        
        return len(self.reads)
    
    
    def __str__(self) -> List[str]:
        """
        Trả về list các string là các reads đã đọc được
        """
        
        return "".join(str(self.reads))
    
    
#reads = Loader.load(filename="data/hemoglobin.fastq")
#print(reads)