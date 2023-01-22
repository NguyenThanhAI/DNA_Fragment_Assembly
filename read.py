from typing import List, Dict
from vertex import Vertex
from edge import Edge


class Read(object):
    
    def __init__(self, sequence: str, read_id: int) -> None:
        """

        Args:
            sequence (str): Chuỗi đại diện cho read
            read_id (int): id của read
        """
        self.sequence: str = sequence
        self.read_id: int = read_id
        self.edges: List[Edge] = []
        
        
    def __getitem__(self, n: int) -> Edge:
        """Lấy cạnh thứ n trong danh sách các cạnh của read hiện tại

        Args:
            n (int): Vị trí của cạnh cần được lấy

        Returns:
            Edge: Cạnh ở vị trí cần được lấy
        """
        return self.edges[n]
    
    
    def change_x(self, x: Edge, z: Edge) -> bool:
        """Thay đổi cạnh cuối cùng trong tập các cạnh thành cạnh mới z nếu x là cạnh cuối

        Args:
            x (Edge): Cạnh cũ ở cuối
            z (Edge): Cạnh mới

        Returns:
            bool: True nếu cạnh cũ được thay đổi, nếu không là False
        """
        
        if self.edges[-1] == x:
            self.edges[-1] = z
            
            # Thêm read hiện tại và tập các read của cạnh mới
            z.reads.append(self)
            
            # Xóa read hiện tại từ cạnh cũ
            if self in x.reads:
                x.reads.remove(self)
                
            return True
        
        return False
    
    
    def change_y(self, y: Edge, z: Edge) -> bool:
        """Thay đổi cạnh đầu tiên trong tập các cạnh của read hiện tại thành cạnh mới z nếu cạnh đầu tiên là y

        Args:
            y (Edge): Cạnh cũ ở đầu
            z (Edge): Cạnh mới

        Returns:
            bool: True nếu cạnh cũ được thay đổi, nếu không là False
        """
        if self.edges[0] == y:
            self.edges[0] = z
            
            # Thêm read hiện tại và tập các read của cạnh mới
            z.reads.append(self)
            
            # Xóa read hiện tại từ cạnh cũ
            if self in y.reads:
                y.reads.remove(self)
                
            return True
        
        return False
    
    
    def change_xy(self, x: Edge, y: Edge, z: Edge) -> bool:
        """Thay đổi hai cạnh liên tiếp x, y trong đường dẫn
    

        Args:
            x (Edge): Cạnh đầu tiên liền trước
            y (Edge): Cạnh liền sau cạnh x
            z (Edge): Cạnh mới thay thế hai cạnh x và cạnh y

        Returns:
            bool: True nếu hai cạnh bị thay đổi bởi z, nếu không là False
        """
        
        found_xy: bool = False
        
        # Xét từng cặp cạnh liền nhau có phải là cạnh x là cạnh đầu tiên, y là cạnh liền cạnh x hay không:
        for i in range(len(self.edges)-1):
            
            # Nếu tìm thấy x là cạnh liền trước, y là cạnh liền sau cạnh x
            if self.edges[i] == x and self.edges[i+1] == y:
                found_xy = True
                
                # Thêm cạnh z, xóa bỏ các cạnh x và cạnh y
                self.edges[i] = z
                self.edges[i+1] = None
                
                # Thêm read hiện tại vào danh sách các read của cạnh mới z
                z.reads.append(self)
                
                # Xóa read hiện tại từ danh sách các read của các cạnh bị xóa x và y
                if self in x.reads:
                    x.reads.remove(self)
                if self in y.reads:
                    y.reads.remove(self)
        
        # Xóa tất cả các cạnh mà vị trí đó là None
        self.edges = [e for e in self.edges if e != None]
        
        return found_xy
    
    
    def update(self, x: Edge, y: Edge, z: Edge) -> bool:
        """Chạy tất cả các phương thức cập nhật cạnh mới

        Args:
            x (Edge): Cạnh nếu là cạnh cuối cùng trong danh sách các cạnh sẽ bị thay thế bởi cạnh z
            y (Edge): Cạnh nếu là cạnh đầu tiên trong danh sách các cạnh sẽ bị thay thế bởi cạnh z
            z (Edge): Cạnh mới

        Returns:
            bool: True nếu có thay đổi, nếu không là False
        """
        
        # Xét từng trường hợp
        changed: List[bool] = [self.change_xy(x=x, y=y, z=z), self.change_x(x=x, z=z), self.change_y(y=y, z=z)]
        
        return any(changed)
    