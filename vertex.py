from typing import List
from edge import Edge


class Vertex(object):
    
    
    def __init__(self, sequence: str) -> None:
        """

        Args:
            sequence (str): Chuỗi đánh dấu đỉnh
        """
        
        self.sequence: str = sequence
        self.in_edges: List[Edge] = []
        self.out_edges: List[Edge] = []
    
    
    def __getitem__(self, n: int) -> str:
        """Lấy ký tự thứ n trong chuỗi đại diện cho đỉnh

        Args:
            n (int): Vị trí n mà ký tự cần lấy

        Returns:
            str: Ký tự tại vị trí được chọn
        """
        
        return self.sequence[n]
    
    
    def __len__(self) -> int:
        """Trả về độ dài của chuỗi đại diện cho đỉnh

        Returns:
            int: Độ dài của chuỗi đại diện cho đỉnh
        """
        
        return len(self.sequence)
    
    
    def __str__(self) -> str:
        """Trả vệ nội dung chuỗi đại diện cho đỉnh

        Returns:
            str: Chuỗi đại diện cho đỉnh
        """
        
        return self.sequence
    
    
    def add_out_edge(self, out_vertex, edge: Edge) -> None:
        """ Thêm một cạnh đi ra từ đỉnh hiện tại

        Args:
            out_vertex (Vertex): Đỉnh đích của cạnh đi ra từ đỉnh hiện tại
            edge (_type_): Cạnh cần được thêm đi ra từ đỉnh hiện tại
        """    
        
        self.out_edges.append(edge)
        out_vertex.in_edges.append(edge)
        