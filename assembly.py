from loader import Loader
from graph import *


class Assembler(object):
    def __init__(self, filename: str, k: int, error_correct: bool=False) -> None:
        """Khởi tạo Assembler

        Args:
            filename (str): File chứa các read
            k (int): Độ dài một k-mer
            error_correct (bool, optional): Có sửa lỗi hay không. Defaults to False.
        """
        
        # Load file
        reads: Loader = Loader.load(filename=filename)
        
        # Khởi tạo đồ thị
        self.graph: Graph = Graph(seqs=reads, k=k, threshold=error_correct)
        self.k: int = k
        
    
    def make_superpath(self) -> None:
        """
        Tạo các superpath bằng cách lặp
        """
        while self.superpath_consider():
            self.graph.clean()
        
        self.graph.clean()
        
        
    def superpath_consider(self) -> bool:
        """Gộp các cạnh nếu có thể

        Returns:
            bool: Có thể gộp được hai cạnh nào đó không
        """
        for read in self.graph.read_list:
            for i in range(len(read.edges)-1):
                # Lấy cạnh phía trước nếu có
                if i > 0:
                    p: Edge = read.edges[i-1]
                else:
                    p = None
                # Lấy hai cạnh x và y kề nhau
                x: Edge = read.edges[i]
                y: Edge = read.edges[i+1]
                # Kiểm tra xem các cạnh có gộp được không
                if self.graph.is_mergeable(p=p, x=x, y=y):
                    if self.graph.merge(x, y):
                        return True
                    
        return False
    
    
    def is_eulerian(self) -> bool:
        """
        Kiểm tra đồ thị có phải là đồ thị Euler hay không

        Returns:
            bool: Nếu là đồ thị Euler trả về True, nếu không là False
        """
        semis: int = 0
        for vertex in self.graph.vertex_list:
            diff = abs(len(vertex.out_edges) - len(vertex.in_edges))
            
            if diff > 1:
                return False
            elif diff == 1:
                semis += 1
        
        # Không phải là đồ thị Euler nếu có hơn 2 đỉnh không cân bằng
        if semis > 2:
            return False
        
        return True
    
    
    def balance(self) -> bool:
        """
        Kết nối các đỉnh bán cân bằng với đỉnh bán cân bằng còn lại.
        Hàm này chỉ thực hiện được khi đồ thị là đồ thị Euler.

        Returns:
            bool: Trả về False nếu không thực hiện được. Trả về True nếu nối thành công
        """
        
        # Tìm các đỉnh không cân bẳng
        semis: List[Vertex] = []
        for vertex in self.graph.vertex_list:
            diff = abs(len(vertex.out_edges) - len(vertex.in_edges))
            if diff == 1:
                semis += [vertex]
                
        if len(semis) != 2:
            print("Quá nhiều đỉnh không cân bằng {}. Đồ thị không phải là đồ thị Euler".format(len(semis)))
            return False
        if len(semis[0].in_edges) > len(semis[0].out_edges):
            self.graph.new_edge(in_vertex=semis[0], out_vertex=semis[1], sequence=semis[0].sequence)
        else:
            self.graph.new_edge(in_vertex=semis[1], out_vertex=semis[0], sequence=semis[1].sequence)
            
        return True
    
    
    def find_eulerian_path(self) -> str:
        """
        Xây dựng đường đi Euler trên đồ thị sử dụng thuật toán Heirholzer

        Returns:
            str: Đường đi Euler kết quả
        """
        
        # Khởi tạo
        current_path: List[Edge] = []
        final_path: List[Edge] = []
        
        # Bắt đầu với đỉnh bán cân bẳng với số cạnh vào nhỏ hơn số cạnh ra
        edge: Edge = None
        for vertex in self.graph.vertex_list:
            diff = abs(len(vertex.out_edges) - len(vertex.in_edges))
            if diff == 1 and len(vertex.in_edges) < len(vertex.out_edges):
                edge = self.graph.get_unvisited(vertex=vertex)
        
        # Nếu không có đỉnh nào như vậy thì chọn cạnh ra của đỉnh đầu tiên
        if not edge:
            edge = self.graph.get_unvisited(self.graph.vertex_list[0])
        
        # Thêm các cạnh vào ngăn xếp
        while edge != None:
            edge.visited = True
            current_path.append(edge)
            edge = self.graph.get_unvisited(edge.out_vertex)
        
        # Lấy tất cả các cạnh khác chưa được đi qua và xây dựng đường đi Euler
        
        while len(current_path) > 0:
            edge = current_path.pop()
            final_path.append(edge)
            edge = self.graph.get_unvisited(edge.in_vertex)
            
            # Lặp tất cả các cạnh chưa được đi qua
            while edge != None:
                edge.visited = True
                current_path.append(edge.out_vertex)
                
        # In kết quả bằng cách thêm vào phía trước danh sách
        sequence: str = ""
        while len(final_path):
            edge = final_path.pop()
            if len(final_path) == 0: # Cạnh cuối cùng
                sequence += edge.sequence
            else:
                s_len: int = len(edge.sequence)
                sequence += edge.sequence[:s_len-self.k+1] # Chỉ lấy các ký tự đầu
                
        return sequence
        


        