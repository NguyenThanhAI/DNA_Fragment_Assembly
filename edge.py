from typing import List, Dict
from vertex import Vertex
from read import Read


class Edge(object):
    
    
    def __init__(self, in_vertex: Vertex, out_vertex: Vertex, sequence: str) -> None:
        """

        Args:
            in_vertex (Vertex): Đỉnh bắt đầu cạnh hiện tại
            out_vertex (Vertex): Đỉnh kết thúc cạnh hiện tại
            sequence (str): Chuỗi đại diện cho cạnh hiện tại
        """
        
        self.sequence: str = sequence
        self.in_vertex: Vertex = in_vertex
        self.out_vertex: Vertex = out_vertex
        self.reads: List[Read] = []
        self.visited: bool = False
        
    
    def __getitem__(self, n: int) -> str:
        """Lấy ký tự thứ n trong chuỗi đại diện cho cạnh

        Args:
            n (int): Vị trí n mà ký tự cần lấy

        Returns:
            str: Ký tự tại vị trí được chọn
        """
        
        return self.sequence[n]
    
    
    def __len__(self) -> int:
        """Trả về độ dài của chuỗi đại diện cho cạnh

        Returns:
            int: Độ dài của chuỗi đại diện cho cạnh
        """
        return len(self.sequence)
    
    
    def __str__(self) -> str:
        """Trả vệ nội dung chuỗi đại diện cho đỉnh

        Returns:
            str: Chuỗi đại diện cho đỉnh
        """
        return self.sequence
    