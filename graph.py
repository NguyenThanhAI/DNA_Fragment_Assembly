import copy
import math
from typing import List, Dict, Optional, Tuple
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
        
        # k_mers là danh sách của các danh sách của các k-mer, corrected_seqs là danh sách của các danh sách của các k-mers đã được chỉnh sửa
        self.corrected_seqs, self.k_mers = self.error_correction(threshold=self.threshold)
            
        # Ghép lại các k-mer khi đã sửa lỗi
        for n in range(len(self.corrected_seqs)):
            new_read: str = ""
            for i in self.corrected_seqs[n]:
                if new_read == "":
                    new_read += i
                else:
                    new_read += i[-1]
            self.corrected_seqs[n] = new_read
            
        
        if error_correct:
            self.seqs: List[str] = self.corrected_seqs
            
        for s in range(len(seqs)):
            # Lấy các read
            seq: str = seqs[s]
            # Tạo object Read
            read: Read = Read(sequence=seq, read_id=s)
            
            # Tạo các đỉnh và các cạnh
            for i in range(len(seq)-k+1):
                # Tạo k-mer
                k_mer: str = seq[i:i+k]
                prefix: str = k_mer[:k-1]
                suffix: str = k_mer[1:]
                
                # Tạo đỉnh tiền tố
                if prefix in self.vertex_dict:
                    p_vertex: Vertex = self.vertex_dict[prefix]
                else:
                    p_vertex: Vertex = self.new_vertex(sequence=prefix)
                    
                # Tạo đỉnh hậu tố
                if suffix in self.vertex_dict:
                    s_vertex: Vertex = self.vertex_dict[suffix]
                else:
                    s_vertex: Vertex = self.new_vertex(sequence=suffix)
                    
                # Tạo cạnh
                if k_mer in self.edge_dict:
                    edge: Edge = self.edge_dict[k_mer]
                else:
                    edge: Edge = self.new_edge(in_vertex=p_vertex, out_vertex=s_vertex, sequence=k_mer)
                    
                # Thêm cạnh vào danh sách cạnh của read
                read.edges.append(edge)
                # Thêm read vào danh sách read của cạnh
                edge.reads.append(read)                
                
        
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
    
    
    def error_correction(self, threshold: int) -> Tuple[List[str], List[str]]:
        """Sửa lỗi các reads và lưu vào đồ thị

        Args:
            threshold (int): Ngưỡng sửa lỗi để một k-mers được gọi là "đặc"

        Returns:
            Tuplle[List[str], List[str]]: Danh sách các k-mers được sửa lỗi và danh sách các k-mers gốc
        """
        
        """
        Là danh sách của từng danh sách của các k-mers, 
        mỗi danh sách con là danh sách của các k-mers từ một read
        """
        list_sequences: List[List[str]] = []
        freq_dict: Dict[str, int] = {}
        
        for read in self.seqs:
            one_sequence_kmers: List[str] = [] # Danh sách các k-mers từ một read
            for i in range(len(read)-(self.k)+1):
                one_sequence_kmers.append(read[i:i+(self.k)])
                # Tính số lần xuất hiện theo k-mers
                if read[i:i+(self.k)] not in freq_dict:
                    freq_dict[read[i:i+(self.k)]] = 1
                else:
                    freq_dict[read[i:i+(self.k)]] += 1
            list_sequences.append(one_sequence_kmers)
        
        seq_list: List[List[str]] = copy.deepcopy(list_sequences)
        
        # Lặp qua từng read của list sequence
        for read in seq_list:
            # Lặp qua từng k-mer trong từng danh sách con (đại diện cho từng read)
            for i in range(len(read)):
                k_mer: str = read[i]
                
                """
                Tra số lần xuất hiện của k-mers trong bảng tần số và 
                kiểm tra số lần xuất hiện có lớn hơn hoặc bằng ngưỡng hay không?
                """
                ocurrences: int = freq_dict[k_mer]
                if ocurrences < threshold:
                    """
                    Nếu số lần xuất hiện trong bảng tần số ít hơn ngưỡng, ta kiểm tra thêm k-1 k-mer tiếp theo và 
                    nếu số lần xuất hiện cũng những k-mer này cũng nhỏ hơn ngưỡng.
                    Nếu k/2 trong số k-mer này số lần xuất hiện cũng nhỏ hơn ngưỡng, 
                    ta cần phải thay đổi, nếu không ta giữ nguyên
                    """
                    cutoff = math.ceil((len(k_mer)/2.0))
                    # Lặp thêm k-1 k-mer sau k-mer hiện tại
                    for j in range(len(k_mer)): # range(1, len(k_mer))?????
                        if i + j >= len(read):
                            break
                        if freq_dict[read[i+j]] < threshold:
                            cutoff -= 1

                        # Nếu quá nhiều k-mer ngay sau dưới ngưỡng, ta cần phải sửa lỗi
                        if cutoff <= 0: # Sai lề???????
                            choices: List[str] = ["A", "C", "G", "T"]
                            best_score: int = 0
                            k_mer_to_change: int = -1 # k-mer đầu tiên cần được thay đổi
                            letter: int = -1 # Chỉ số vị trí của ký tự cần bị thay đổi trong tập các lựa chọn A, C, G, T

                            # Chỉ lặp qua một nửa số ký tự trong k-mer
                            for x in range(int(math.ceil(len(k_mer)/2.0))):
                                index: int = -1 - x
                                # Theo dõi score của từng ký tự
                                scores: List[int] = [0, 0, 0, 0]
                                for y in range(len(k_mer)-x):
                                    # Lặp qua từng k-mer cũng bao gồm ký tự để so sánh score
                                    k_mer_index: int = (i) + y
                                    if k_mer_index >= len(read) or k_mer_index < 0:
                                        continue
                                    for n in range(len(choices)):
                                        option: str = str(read[k_mer_index])
                                        option_list: List[str] = list(option)
                                        new_index: int = index - y
                                        if new_index < -len(k_mer):
                                            new_index += len(k_mer)
                                        option_list[new_index] = choices[n]
                                        option_back_to_string: str = "".join(option_list)
                                        option = option_back_to_string
                                        if option in freq_dict:
                                            # Thêm vào chênh lệch về số lần xuất hiện của k-mer thay thế và k-mer ban đầu
                                            scores[n] += (freq_dict[option] - freq_dict[read[k_mer_index]])
                                        else:
                                            # Trừ đi nếu không trong bảng tần số
                                            scores[n] -= 10

                                for n in range(len(scores)):
                                    if scores[n] > best_score:
                                        best_score = scores[n]
                                        k_mer_to_change = (i - x)
                                        letter = n

                            if best_score == 0:
                                # Giữ nguyên ban đầu
                                continue
                            else:
                                # Thay đổi k-mer
                                for p in range(len(k_mer)): # Thay đổi thành -x
                                    if (p + k_mer_to_change) < len(read):
                                        freq_dict[read[k_mer_to_change + p]] -= 1
                                        list_of_k_mer: List[str] = list(read[k_mer_to_change + p])
                                        list_of_k_mer[-(p+1)] = choices[letter]
                                        back_to_string = "".join(list_of_k_mer)
                                        read[k_mer_to_change + p] = back_to_string
                                        if read[k_mer_to_change + p] in freq_dict:
                                            freq_dict[read[k_mer_to_change + p]] += 1
                                        else:
                                            freq_dict[read[k_mer_to_change + p]] = 1
                            
        return seq_list, list_sequences
                    