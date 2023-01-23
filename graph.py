import copy
import math
from typing import List, Dict, Optional
from loader import Loader
from vertex import Vertex
from edge import Edge
from read import Read



class Graph(object):
    
    def __init__(self, seqs: Optional[Loader], k: int, threshold: int, error_correct: bool = False) -> None:
        """

        Args:
            seqs (Loader): Các reads đọc được trong loader
            k (int): k-mers, số ký tự trong một chuỗi đại diện cho một cạnh
            threshold (int): ngưỡng để sửa lỗi
            error_correct (bool, optional): Có sử lỗi hay không. Defaults to False.
        """
        
        self.vertex_list: List[Vertex] = [] # Danh sách các đỉnh trong đồ thị
        self.vertex_dict: Dict[str, Vertex] = {} # Danh sách các đỉnh được đánh chỉ mục bởi chuỗi đại diện
        self.edge_list: List[Edge] = [] # Danh sách các cạnh trong đồ thị
        self.edge_dict: Dict[str, Edge] = {} # Danh sách các cạnh trong đồ thị được chỉ mục bởi chuỗi đại diện
        self.read_list: List[Read] = [] # Danh sách các reads trong đồ thị
        self.k: int = k # Độ dài chuỗi đại diện cho một cạnh
        self.seqs: Optional[Loader] = seqs # Các read được đọc từ loader
        self.threshold: int = threshold # Ngưỡng để sửa lỗi
        
        
    def __str__(self) -> str:
        """_summary_

        Returns:
            str: In ra các cạnh và các đỉnh kề
        """
        
        out: str = ""
        for edge in self.edge_list:
            out = out + str(edge) + ": " + str(edge.in_vertex) + ", " + str(edge.out_vertex) + "\n"
            
        return out
    
    
    def new_vertex(self, sequence: str) -> Vertex:
        """Tạo ra một đỉnh mới thêm vào đồ thị

        Args:
            sequence (str): Chuỗi đại diện cho đỉnh

        Returns:
            Vertex: Đỉnh mới được tạo ra
        """
        
        vertex: Vertex = Vertex(sequence=sequence)
        self.vertex_list.append(vertex)
        self.vertex_dict[sequence] = vertex
        
        return vertex
    
    
    def new_edge(self, in_vertex: Vertex, out_vertex: Vertex, sequence: str) -> Edge:
        """Tạo ra một cạnh mới thêm vào đồ thị khi cho biết đỉnh vào, đỉnh ra, chuỗi đại diện cho cạnh

        Args:
            in_vertex (Vertex): Đỉnh vào cạnh mới
            out_vertex (Vertex): Đỉnh ra cạnh mới
            sequence (str): Chuỗi đại diện cho cạnh mới

        Returns:
            Edge: Cạnh mới được tạo ra
        """
        
        edge: Edge = Edge(in_vertex=in_vertex, out_vertex=out_vertex, sequence=sequence)
        self.edge_list.append(edge)
        in_vertex.add_out_edge(out_vertex=out_vertex, edge=edge)
        self.edge_dict[sequence] = edge
        
        return edge
    
    
    def get_unvisited(self, vertex: Vertex) -> Edge:
        """Lấy cạnh ra của một đỉnh mà chưa được đi qua

        Args:
            vertex (Vertex): đỉnh được tìm trong các cạnh ra tương ứng

        Returns:
            Edge: cạnh ra đầu tiên chưa được đi qua của đỉnh là tham số đầu vào của hàm
        """
        
        for edge in vertex.out_edges:
            if not edge.visited:
                return edge
            
        return None
    
    
    def merge(self, x: Edge, y: Edge) -> Edge:
        """Gộp hai cạnh kề nhau x và y

        Args:
            x (Edge): Cạnh liền trước
            y (Edge): Cạnh liền sau

        Returns:
            Edge: Cạnh kết quả là cạnh được gộp hai cạnh x và cạnh y
        """
        
        # Lấy các đỉnh là đỉnh vào của cạnh x, đỉnh ra của cạnh x và đỉnh ra của cạnh y
        in_vertex: Vertex = x.in_vertex
        mid_vertex: Vertex = x.out_vertex
        out_vertex: Vertex = y.out_vertex
        
        assert mid_vertex == y.in_vertex
        
        # Kiểm tra các đỉnh vẫn còn kích hoạt (không có đỉnh nào không có cạnh vào hoặc cạnh ra)
        
        if len(in_vertex.out_edges) == 0 or len(mid_vertex.in_edges) == 0 or \
            len(mid_vertex.out_edges) == 0 or len(out_vertex.in_edges) == 0:
                return None
            
        # Tạo chuỗi đại diện mới cho cạnh mới        
        y_length: int = len(y.sequence)
        seq: str = x.sequence + y.sequence[self.k-1:]
        
        # Tạo một cạnh mới        
        z: Edge = self.new_edge(in_vertex=in_vertex, out_vertex=out_vertex, sequence=seq)
        
        # Cập nhật các đỉnh và đường đi
        if x in in_vertex.out_edges:
            in_vertex.out_edges.remove(x)
        if x in mid_vertex.in_edges:
            mid_vertex.in_edges.remove(x)
        if y in mid_vertex.out_edges:
            mid_vertex.out_edges.remove(y)
        if y in out_vertex.in_edges:
            out_vertex.in_edges.remove(y)
        for read in self.read_list:
            read.update(x=x, y=y, z=z)
            
        return z
        
        
    def is_mergeable(self, p: Edge, x: Edge, y: Edge) -> bool:
        """Xem có thể gộp hai cạnh x và y thành một cạnh được không

        Args:
            p (Edge): Cạnh kề trước cạnh x
            x (Edge): Cạnh x
            y (Edge): Cạnh kề sau cạnh x

        Returns:
            bool: Nếu có thể gộp hai cạnh x và y, trả về True. Nếu không trả về False
        """
        
        is_spanned: bool = False
        # Kiểm tra hai cạnh x và y có cùng một read hay không
        for read in x.reads:
            if read in y.reads:
                is_spanned = True
                
        # Kiểm tra xem có tồn tại một read mà thông tin bị mâu thuẫn không?
        is_conflicted: bool = False
        for read in x.reads:
            # Bắt đầu từ cạnh kề trước, theo thường đi này
            for i in range(len(read.edges)-2):
                # Tìm x
                if read.edges[i] == x:
                    # Bắt đầu khớp
                    if (p == None and i == 0) or p == read.edges[i-1]:
                        # Kiểm tra nếu y không phải ở vị trí sau x
                        if y != read.edges[i+1]:
                            is_conflicted = True
                    else:
                        
                        if y == read.edges[i+1]:
                            is_conflicted = True
                            
            return is_spanned and not is_conflicted # Sai căn lề??????????
        

    def clean(self) -> None:
        """Loại bỏ các cạnh và đỉnh trống từ đồ thị
        """
        
        # Loại bỏ các đỉnh rỗng
        to_remove_vertex: List[Vertex] = []
        for vertex in self.vertex_list:
            if len(vertex.in_edges) == 0 and len(vertex.out_edges) == 0:
                to_remove_vertex.append(vertex)
                
        # Xóa các đỉnh rỗng
        for vertex in to_remove_vertex:
            self.vertex_list.remove(vertex)
            
        # Loại bỏ các cạnh rỗng
        to_remove_edge: List[Edge] = []
        for edge in self.edge_list:
            if len(edge.reads) == 0:
                to_remove_edge.append(edge)
            elif edge.in_vertex not in self.vertex_list:
                to_remove_edge.append(edge)
            elif edge.out_vertex not in self.vertex_list:
                to_remove_edge.append(edge)
                
        # Xóa các cạnh rỗng
        for edge in to_remove_edge:
            self.edge_list.remove(edge)
            