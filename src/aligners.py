'''
@author: ahmed allam <ahmed.allam@yale.edu>
'''
import numpy as np
import scipy.spatial.distance as scpydist

def transform_scale(score, xref_min, xref_max):
    """transforms score from [0,1] scale  to [-1,1]"""
    x_new = 2*(score-xref_min)/(xref_max-xref_min) - 1
#     x_new = score-0.5
    return(x_new)

class GenericSeqNode():
    """Generalizing elements in Sequences to nodes
    
       This class generalizes the character elements in sequences where each element 
       is represented by a character from some alphabet. It represents elements in the sequence
       as nodes.
       
       Args:
           label:
           trans_time:
           
       Attributes:
           label:
           trans_time:
    """
    __slots__ = ['label', 'trans_time']
    def __init__(self, label, trans_time=None):
        self.label = label
        self.trans_time = trans_time
    def compute_similarity(self, seq_node, **kwargs):
        """compute similarity between two sequence nodes
        
          ..note::
          
              to be implemented in the child class
        """
        pass
    
    # to make it a class method .. 
    def compute_temporal_penalty(self, accum_tseq_1, accum_tseq_2, **kwargs):
        """implement temporal penalty function -- in case of using temporal alignment
        """
        # implement the penalty function
        T_penalty = kwargs.get("T_penalty")
        if(not T_penalty):
            T_penalty = 0.25
        s = 0
        if(accum_tseq_1 or accum_tseq_2):
            s = T_penalty * np.abs(accum_tseq_1-accum_tseq_2)/max(accum_tseq_1, accum_tseq_2)
        return(s)
    
class CharSeqNode(GenericSeqNode):
    """Character node in a sequence
    """
    def __init__(self, label, trans_time=None):
        super().__init__(label, trans_time)
        
    def compute_similarity(self, seq_node, **kwargs):
        # get the score based on the scoring matrix/substitution matrix
        score_matrix = kwargs.get('score_matrix')
        return(score_matrix[self.label, seq_node.label])
    
class SeqNode(GenericSeqNode):
    """Node in a sequence
    
       This class generalizes the character elements in sequences where usually each element 
       is represented by a character from some alphabet. 
       The aim of this class is to represent elements in the sequence as nodes where each node 
       is represented by certain properties in a feature vector.
       
       Args:
           label:
           trans_time:
           feat_vec:
       Attributes:
           label:
           trans_time:
           feat_vec:
    """
    __slots__ = ['feat_vec']
    def __init__(self, label, feat_vec, trans_time=None):
        super().__init__(label, trans_time)
        self.feat_vec = np.array(feat_vec)
        
#     def compute_similarity(self, seq_node, **kwargs):
#         # parse the method for computing the score
#         if(not kwargs.get('score_type')):
#             score_type = 'euclidean'
#         # compute a score
#         featvec_a = self.feat_vec
#         featvec_b = seq_node.feat_vec
#         if(score_type == 'euclidean'):
#             distance = np.linalg.norm(featvec_a - featvec_b)
#             score = 1/(1+distance)
#         elif(score_type == 'manhatten'):
#             distance = np.abs(featvec_a - featvec_b).sum()
#             score = 1/(1+distance)
#         elif(score_type == 'cosine'):
#             try:
#                 score = np.dot(featvec_a, featvec_b)/(np.linalg.norm(featvec_a)*np.linalg.norm(featvec_b))
#             except ZeroDivisionError:
#                 print('dividing by zero..')
#             finally:
#                 return(0)
#         elif(score_type == 'jaccard'):
#             # get activated features
#             set_a_activef = set(np.where(featvec_a!=0)[0])
#             set_b_activef = set(np.where(featvec_b!=0)[0])
#             n_activef = len(set_a_activef.intersection(set_b_activef))
#             # get non active features
#             set_a_nonactivef = set(np.where(featvec_a==0)[0])
#             set_b_nonactivef = set(np.where(featvec_b==0)[0])
#             n_nonactivef = len(set_a_nonactivef.intersection(set_b_nonactivef))
#             score = (n_activef+n_nonactivef)/(len(self.feat_vec))
#         return(score)

    def compute_similarity(self, seq_node, **kwargs):
        """compute similarity between two feature vector representation using :module:``scipy.spatial.distance``
           Args:
               seq_node: instance of :class:`GenericSeqNode` or :class:`CharSeqNode` or :class:`SeqNode`
           Keyword args:
               score_type: string, scoring option from {'euclidean', 'manhatten', 'cosine', 'jaccard', 'correlation'}
               normalize_vec: bool, normalize feature vector components
               transform_score: bool, transform score to -1,1 range
               
        """
        # compute a score
        featvec_a = np.asarray(self.feat_vec)
        featvec_b = np.asarray(seq_node.feat_vec)
        # parse the method for computing the score
        score_type = kwargs.get('score_type')
        if(not score_type):
            score_type = 'euclidean'
            
        normalize_vec = kwargs.get('normalize_vec')
        if(normalize_vec):
            featvec_a = featvec_a/np.sum(featvec_a)
            featvec_b = featvec_b/np.sum(featvec_b)
            
        if(score_type == 'euclidean'):
            distance = scpydist.euclidean(featvec_a, featvec_b)
            if(normalize_vec):
                score = 1-distance
            else:
                score = 1/(1+distance)
        elif(score_type == 'manhatten'):
            distance = scpydist.minkowski(featvec_a, featvec_b, 1)
            if(normalize_vec):
                score = 1-distance
            else:
                score = 1/(1+distance)
        elif(score_type == 'cosine'):
            try:
                distance = scpydist.cosine(featvec_a, featvec_b)
                score = 1-distance
            except ZeroDivisionError:
                print('dividing by zero..')
            finally:
                return(0)
        elif(score_type == 'jaccard'):
            # get activated features
            distance = scpydist.jaccard(featvec_a, featvec_b)
            score = 1-distance
        elif(score_type == 'correlation'):
            distance = scpydist.correlation(featvec_a, featvec_b)
            score = 1-distance
        change_scale = kwargs.get('transform_score')
        if(score_type != 'correlation' and change_scale):
            score = transform_scale(score, 0, 1)
        return(score)

class ScoringMatrices():
    __slots__=[]
    def __init__(self):
        pass
    @staticmethod
    def get_common_scoring_matrix(matrix_name):
        """retrieve common substitution matrices (i.e. BlOSUM, PAM ..)
        
           check these relevant links:
               -`wikipedia page <https://en.wikipedia.org/wiki/Substitution_matrix>`__
               - http://biopython.org/DIST/docs/api/Bio.SubsMat.MatrixInfo-module.html
               - https://web.archive.org/web/19991014010917/http://www.embl-heidelberg.de/%7Evogt/matrices/mlist1.html
           --
           TODO: generate common substitution matrices
           
        """
        pass
    
class GapsParams():
    """ defines gaps parameters required for alignment
    
        Args:
            gap_open: number, penalty for opening a gap 
            gap_ext: number or None (default), penalty for extending an existing gap.
                     if gap_ext is None, the scoring assumes a linear gap penalty,
                     else the scoring assumes affine gap penalty              
    """
    __slots__ = ['gap_open', 'gap_ext']
    def __init__(self, gap_open, gap_ext=None):
        # to add method to verify that gap penalties are negative numbers
        self.gap_open = gap_open
        self.gap_ext = gap_ext
        
    
class Aligner():
    """ Generic sequence aligner
    
        Args:
            align_type: string, defining type of alignment.
                        There are four alignment types:
                            - `global` (i.e. Needlman-Wunch or more precisely a variant of David Sankoff algorithm
                                             see  Sankoff, D. (1972). "Matching sequences under deletion-insertion constraints". 
                                                  Proceedings of the National Academy of Sciences of the United States of America. 69 (1): 4–6
                                        )
                            - `local` (i.e. Smith-Waterman algorithm)
                            - `semi-global` alignment
                            - `end-gap-free` alignment
        TODO:
            to implement overlap alignment such as the following pattern:
                +++++++******--------
                -------******========
                - is the gap and +,*,= are the symbols/characters
            which is roughly performed by assigning 0 to the first column and take the max score in the last row of the dynamic table
                            
        ..note ::
        
            After reading multiple articles and implementations, I came to a conclusion that there were inconsistencies
            in many of the implementations and/or the formulation of the dynamic programming recursion -- especially
            for the affine gap penalty case. Interestingly -- after figuring out that -- I found a paper discussing/reporting
            the same frustration I witnessed. 
                    
                    ` Are all global alignment algorithms and implementations correct? <http://www.biorxiv.org/content/biorxiv/early/2015/11/12/031500.full.pdf>`__
            
            This article should be THE reference paper for debugging/verifying alignment implementation -- especially for affine gap penalty. 
            
    """
    def __init__(self, align_type):
        self.align_type = align_type
        
    def align(self, seq_1, seq_2, gaps_params, **kwargs):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        align_type = self.align_type
        # parse the gaps -- assume that gaps are negative numbers
        gap_open = gaps_params.gap_open
        gap_ext = gaps_params.gap_ext
        l_seq1 = len(seq_1)
        l_seq2 = len(seq_2)
#         print("kwargs ", kwargs)
        if(gap_ext == None): # case of linear gap applied
            # backtrack matrix
            Vp = self._init_pointer_table(l_seq1, l_seq2)
            # initialize V matrix -- holding the best alignment score for seq1[1:i] and seq2[1:j]
            V = np.zeros((l_seq1+1, l_seq2+1), dtype='float32')
            if(align_type == 'global'):
                # fill the top row  with the linear gap
                V[0,1:] = gap_open * np.arange(1, l_seq2+1, dtype='float32')
                # fill the first left column with the linear gap
                V[1:,0] = gap_open * np.arange(1, l_seq1+1, dtype='float32')
               
            elif(align_type == 'semi-global'):
                # in this implementation, seq_2 is always considered the pattern sequence (i.e. generally the pattern sequence is the shorter sequence)
                # we penalize gaps in the pattern sequence (seq_2) but not in the long sequence (seq_1)
                # fill the top row with the linear gap
                V[0, 1:] = gap_open * np.arange(1, l_seq2+1, dtype='float32')    
            for i in range(1, l_seq1+1):
                for j in range(1, l_seq2+1):
                    # compute score
                    node_seq1 = seq_1[i-1]
                    node_seq2 = seq_2[j-1]
                    sim_score = node_seq1.compute_similarity(node_seq2, **kwargs)
                    case_1 = V[i, j-1] + gap_open # align node_j from seq_2 to a gap
                    case_2 = V[i-1, j-1] + sim_score # align node_i from seq_1 to node_j from seq_2
                    case_3 = V[i-1, j] + gap_open # align node_i from seq_1 to a gap 
                    cases_vec = [case_1, case_2, case_3]
#                     print('i,j = ', (i,j))
#                     print('cases_vec ', cases_vec)
                    self._fill_dynamic_and_pointer_tables(cases_vec, V, Vp, (i,j), main_dtable = True)

            return(V, Vp, i, j)

        else: # case of affine gap applied
            # backtrack matrices
            Vp = self._init_pointer_table(l_seq1, l_seq2)
            Fp = self._init_pointer_table(l_seq1, l_seq2)
            Ep = self._init_pointer_table(l_seq1, l_seq2)
            # initialize V matrix -- holding the best alignment score for seq1[1:i] and seq2[1:j]
            V = np.zeros((l_seq1+1, l_seq2+1), dtype='float32')
            # best alignment score where node_i in seq_1 is aligned to a space
            F = np.ones((l_seq1+1, l_seq2+1), dtype='float32') * -np.inf 
            # best alignment score where node_j in seq_2 is aligned to a space
            E = np.ones((l_seq1+1, l_seq2+1), dtype='float32') * -np.inf
              
            if(align_type == 'global'):
                # the affine gap function used
                # g(k) = -gap_open-gap_ext*k where k is the gap length
                # fill the top row  with the affine gap
                V[0, 1:] = gap_open + gap_ext * np.arange(1, l_seq2+1, dtype='float32')
                # fill the first left column with the affine gap
                V[1:, 0] = gap_open + gap_ext * np.arange(1, l_seq1+1, dtype='float32')
                  
            elif(align_type == 'semi-global'):
                # seq_2 is always considered the pattern sequence (i.e. pattern sequence should be the shorter sequence)
                # we penalize gaps in the pattern sequence (seq_2) but not in the long sequence (seq_1)
                # fill the top row with the linear gap
                V[0, 1:] = gap_open * np.arange(1, l_seq2+1, dtype='float32')
            for i in range(1, l_seq1+1):
                for j in range(1, l_seq2+1):
                    # compute score
                    node_seq1 = seq_1[i-1]
                    node_seq2 = seq_2[j-1]
                    sim_score = node_seq1.compute_similarity(node_seq2, **kwargs)

                    e_vec = [E[i, j-1] + gap_ext, 
                             V[i, j-1] + gap_open + gap_ext, 
                             F[i, j-1] + gap_open + gap_ext]
                    self._fill_dynamic_and_pointer_tables(e_vec, E, Ep, (i,j))
                    
                    f_vec = [E[i-1, j] + gap_open + gap_ext,
                             V[i-1, j] + gap_open + gap_ext, 
                             F[i-1, j] + gap_ext]
                    self._fill_dynamic_and_pointer_tables(f_vec, F, Fp, (i,j))

                    v_vec = [E[i-1, j-1] + sim_score,
                             V[i-1, j-1] + sim_score, 
                             F[i-1, j-1] + sim_score] 
                    self._fill_dynamic_and_pointer_tables(v_vec, V, Vp, (i,j), main_dtable = True)
           
            print("E ", E)
            print("Ep ", Ep)
            print("-"*50)
            print("F ", F)
            print("Fp ", Fp)
            print("-"*50)
            print("V ", V)
            print("Vp ", Vp)
            print("-"*50)
            return(V, Vp, E, Ep, F, Fp, i, j)

    def _fill_dynamic_and_pointer_tables(self, vec, dynamic_table, pointer_table, pos, main_dtable = False):
        direc_flags = np.zeros(3)
        vec_max = np.max(vec)
        if(main_dtable and self.align_type == 'local'): # case of V matrix and local alignment is required
            vec_max = max(vec_max, 0)
        direc_flags[np.where(vec == vec_max)[0]] = 1
        i, j = pos
        pointer_table[i,j] = tuple(direc_flags)
        dynamic_table[i,j] = vec_max
        
    def _init_pointer_table(self, l_seq1, l_seq2):
        P = {}
        for i in range(1, l_seq1+1):
            P[i, 0] = (0,0,1)
        for j in range(1, l_seq2+1):
            P[0,j] = (1,0,0)
        return(P)

    def retrieve_alignments(self, *args, **kwargs):
        res = None
        if(len(args) == 4):
            res = self._retrieve_alignments_linear(*args, **kwargs)
        else:
            res = self._retrieve_alignments_affine(*args, **kwargs)
        return(res)
            
    def _retrieve_alignments_linear(self, V, Vp, end_i, end_j, num_paths=1):
        align_type = self.align_type
        # get the starting nodes with max align score corresponding to the specified alignment type
        root_nodes = []
        if(align_type == 'global'):   
            root_nodes.append((end_i, end_j))
        elif(align_type == 'semi-global'):
            # maximum in the last column
            max_lcol = np.max(V[:,-1])
            for indx_row in np.where(V[:,-1] == max_lcol)[0]:
                root_nodes.append((indx_row, end_j))
        elif(align_type == 'end-gap-free'):
            max_lrow = np.max(V[-1,:]) # max last row
            max_lcol = np.max(V[:,-1]) # max last column
            if(max_lrow == max_lcol):
                for indx_col in np.where(V[-1,:]==max_lrow)[0]:
                    root_nodes.append((end_i, indx_col))
                for indx_row in np.where(V[:,-1]==max_lcol)[0]:
                    root_nodes.append((indx_row, end_j))
            elif(max_lrow > max_lcol):
                for indx_col in np.where(V[-1,:]==max_lrow)[0]:
                    root_nodes.append((end_i, indx_col))
            elif(max_lrow < max_lcol):
                for indx_row in np.where(V[:,-1]==max_lcol)[0]:
                    root_nodes.append((indx_row, end_j))
        elif(align_type == 'local'):
            max_matrix = np.max(V)
            for coord in np.argwhere(V==max_matrix):
                root_nodes.append(tuple(coord))
#         print("root_nodes ", root_nodes)
        score = V[root_nodes[0]]
        if(num_paths == 0): # returning only the score
            return((score, ))
        alignments_path = []
        # building the graph
        for root_node in root_nodes:
            # building a graph from the nodes (i.e. coordinates (i,j))
            coord_graph = {}
            terminal_nodes = set()
            self._build_coord_graph_linear(root_node, V, Vp, coord_graph, terminal_nodes)
#             print("coord_graph ", coord_graph)
#             print("terminal_nodes ", terminal_nodes)
            # prepare for traversing the graph using depth-first 
            all_paths = self._traverse_path(root_node, coord_graph, gap_type='linear')
            alignments_path += self._build_alignments(all_paths, terminal_nodes)
            del all_paths
            if(len(alignments_path)>=num_paths):
                return((score, alignments_path))
        return((score, alignments_path))
                        
    def _build_coord_graph_linear(self, curr_node, V, Vp, coord_graph, terminal_nodes):
        i, j = curr_node
        conds = {'global' : (i==0 and j==0),
                 'local' : V[curr_node]==0,
                 'semi-global' : (i==0 and j==0),
                 'end-gap-free' : (i==0 and j==0)}
        
        if(conds[self.align_type]):
            coord_graph[curr_node] = []
            terminal_nodes.add(curr_node)
            return
        if(curr_node not in coord_graph): # explore the children of parent node
            # get children nodes
            children = []
            for pos, direc in enumerate(Vp[curr_node]):
                if(direc != 0):
                    if(pos==0): # case of horizontal direction (i.e. left -- (i, j-1))
                        children.append((i, j-1))
                    elif(pos==1): # case of diagonal direction (i-1, j-1)
                        children.append((i-1, j-1))
                    else: # case of vertical direction (i-1, j)
                        children.append((i-1, j))
            coord_graph[curr_node] = children
            for child in children:
                self._build_coord_graph_linear(child, V, Vp, coord_graph, terminal_nodes)
                
    def _retrieve_alignments_affine(self, V, Vp, E, Ep, F, Fp, end_i, end_j, num_paths=1):
        align_type = self.align_type
        # get the starting nodes with max align score corresponding to the specified alignment type
        root_nodes = []
        if(align_type == 'global'):   
            vec = [np.max(E[end_i, end_j]), np.max(V[end_i, end_j]), np.max(F[end_i, end_j])]
            v_max = np.max(vec)
            dtable_name = None
            for indx_pos in np.where(vec == v_max)[0]:
                if(indx_pos==0):
                    dtable_name = 'E'
                    score = E[end_i, end_j]
                elif(indx_pos==1):
                    dtable_name = 'V'
                    score = V[end_i, end_j]
                elif(indx_pos==2):
                    dtable_name = 'F'
                    score = F[end_i, end_j]
                root_nodes.append(((end_i, end_j), dtable_name))
        elif(align_type == 'semi-global'):
            max_lcol = np.max(V[:,-1])
            for indx_row in np.where(V[:,-1] == max_lcol)[0]:
                root_nodes.append((indx_row, end_j), 'V')
            score = V[root_nodes[0]]
        elif(align_type == 'end-gap-free'):
            max_lrow = np.max(V[-1,:])
            max_lcol = np.max(V[:,-1])
            if(max_lrow == max_lcol):
                for indx_col in np.where(V[-1,:]==max_lrow)[0]:
                    root_nodes.append((end_i, indx_col), 'V')
                for indx_row in np.where(V[:,-1]==max_lcol)[0]:
                    root_nodes.append((indx_row, end_j), 'V')
            elif(max_lrow > max_lcol):
                for indx_col in np.where(V[-1,:]==max_lrow)[0]:
                    root_nodes.append((end_i, indx_col), 'V')
            elif(max_lrow < max_lcol):
                for indx_row in np.where(V[:,-1]==max_lcol)[0]:
                    root_nodes.append((indx_row, end_j), 'V')
            score = V[root_nodes[0]]
        elif(align_type == 'local'):
            max_matrix = np.max(V)
            for coord in np.argwhere(V==max_matrix):
                root_nodes.append((tuple(coord),'V'))
            score = V[root_nodes[0]]

#         print("root_nodes ", root_nodes)
        
        if(num_paths == 0): # returning only the score
            return((score, ))
            
        # building the graph
        alignments_path = []
        for root_node in root_nodes:
            # building a graph from the nodes (i.e. coordinates (i,j))
            coord_graph = {}
            terminal_nodes = set()
            self._build_coord_graph_affine(root_node, V, Vp, Ep, Fp, coord_graph, terminal_nodes)
#             print("coord_graph ", coord_graph)
#             print("terminal_nodes ", terminal_nodes)
            # prepare for traversing the graph using depth-first 
            all_paths = self._traverse_path(root_node, coord_graph, gap_type='affine')
            alignments_path += self._build_alignments(all_paths, terminal_nodes)
            del all_paths
#             print("num_paths ", num_paths)
            if(len(alignments_path)>=num_paths):
                return((score, alignments_path))
        return((score, alignments_path))
    
    def _build_coord_graph_affine(self, curr_node, V, Vp, Ep, Fp, coord_graph, terminal_nodes):
        # parse node
        curr_pos, pointer = curr_node
        i, j = curr_pos
        conds = {'global' : (i==0 and j==0),
                 'local' : V[curr_pos]==0,
                 'semi-global' : (i==0 and j==0),
                 'end-gap-free' : (i==0 and j==0)}
        
        if(conds[self.align_type]):
            coord_graph[curr_node] = []
            terminal_nodes.add(curr_node)
            return
        if(curr_node not in coord_graph): # explore the children of parent node
            # get children nodes
            children = []
            pointer_table = None
            if(pointer == 'V'):
                pointer_table = Vp
                if(i==0 and j!=0):
                    new_pos = (i, j-1)
                elif(i!=0 and j==0):
                    new_pos = (i-1, j) 
                else:
                    new_pos = (i-1, j-1)
            elif(pointer == 'E'):
                pointer_table = Ep
                if(i==0 and j!=0):
                    new_pos = (i, j-1)
                elif(i!=0 and j==0):
                    new_pos = (i-1, j) 
                else:
                    new_pos = (i, j-1)  
            elif(pointer == 'F'):
                pointer_table = Fp
                if(i==0 and j!=0):
                    new_pos = (i, j-1)
                elif(i!=0 and j==0):
                    new_pos = (i-1, j) 
                else:
                    new_pos = (i-1, j)                
            for indx_pos, direc in enumerate(pointer_table[curr_pos]):
                if(direc != 0):
                    if(indx_pos==0): # case of horizontal direction (i.e. left -- (i, j-1))
                        children.append((new_pos, 'E'))
                    elif(indx_pos==1): # case of diagonal direction (i-1, j-1)
                        children.append((new_pos, 'V'))
                    else: # case of vertical direction (i-1, j)
                        children.append((new_pos, 'F'))
            coord_graph[curr_node] = children
            for child in children:
                self._build_coord_graph_affine(child, V, Vp, Ep, Fp, coord_graph, terminal_nodes)
       
    def _build_alignments(self, paths, terminal_nodes):
        alignments_path = []
#         print("inside _build_alignments")
#         print("terminal_nodes ", terminal_nodes)
#         print("paths ", paths)
        for node in terminal_nodes:
            for i in range(len(paths[node])):
                path = paths[node][i][:]
                path.reverse()
                alignments_path.append(path)
        return(alignments_path)
             
    def _traverse_path(self, root_node, coord_graph, gap_type = 'linear'):
        # breadth first search
        q = [] # start a queue 
        q.append(root_node) # add the root node
        q_parents = [{}] # queue for the parents
        visited_path = {}
        track_paths = {} # track paths
        while(q):
            curr_node = q.pop(0)
            curr_parent = q_parents.pop(0)
#             print("current node ", curr_node)
            if(curr_node in curr_parent): # to protect against root node
#                 print("parent ", curr_parent)
                parent = curr_parent[curr_node]
                if((curr_node, parent) not in visited_path):
                    align_repr = self._get_align_repr(parent, curr_node, gap_type)
#                     print("align_rep ", align_repr)
                    if(parent in track_paths):
                        for path in track_paths[parent]:
#                             print("path ", path)
                            tmp = path[:]
                            tmp.append(align_repr)
                            if(curr_node in track_paths):
                                track_paths[curr_node].append(tmp)
                            else:
                                track_paths[curr_node] = [tmp]
#                             print("track_path[{}] = {}".format(curr_node, track_paths[curr_node]))
                    else:
                        track_paths[curr_node] = [[align_repr]]
#                         print("track_path[{}] = {}".format(curr_node, track_paths[curr_node]))
                    
                    visited_path[(curr_node, parent)] = True
                
            for child in coord_graph[curr_node]:
                q.append(child)
                q_parents.append({child:curr_node})
#             print("q ", q)
#             print("q_parents ", q_parents)
        return(track_paths)

    def _get_align_repr(self, parent_node, child_node, gap_type = 'linear'):
#         print('parent_node ', parent_node)
#         print('child_node ', child_node)
        if(gap_type == 'linear'):
            p_i, p_j = parent_node
            c_i, c_j = child_node
        elif(gap_type == 'affine'):
            p_i, p_j = parent_node[0]
            c_i, c_j = child_node[0]
        seq_1  = self.seq_1
        seq_2 = self.seq_2
        align_repr = ()
        if((p_i-1 == c_i) and (p_j-1 == c_j)): # diagonal case
            align_repr = (seq_1[p_i-1].label, seq_2[p_j-1].label) 
        elif((p_i-1 == c_i) and (p_j == c_j)): # vertical case
            align_repr = (seq_1[p_i-1].label, '-') 
        elif((p_i == c_i) and (p_j-1 == c_j)): # horizontal case
            align_repr = ('-', seq_2[p_j-1].label)
        return(align_repr) 
    
class TemporalAligner(Aligner):
    """ Generic temporal sequence aligner
    
        Implements the method reported in:
            `Syed et al. Temporal Needleman-Wunsch <http://ieeexplore.ieee.org/document/7344785/>`__
            
        Currently, it supports linear gap penalty but extends/adapts the notation in the original paper to support all alignment paths
        ---
        
        TODO: 
            extend it to affine gap penalty
    
        Args:
            align_type: string, defining type of alignment.
                        There are four alignment types:
                            - `global` (i.e. Needleman-Wunsch or more precisely a variant of David Sankoff algorithm
                                             see  Sankoff, D. (1972). "Matching sequences under deletion-insertion constraints". 
                                                  Proceedings of the National Academy of Sciences of the United States of America. 69 (1): 4–6
                                        )
                            - `local` (i.e. Smith-Waterman algorithm)
                            - `semi-global` alignment
                            - `end-gap-free` alignment
                            
    """
    def __init__(self, align_type):
        super().__init__(align_type)
        
    def align(self, seq_1, seq_2, gaps_params, **kwargs):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        align_type = self.align_type
        # parse the gaps -- assume that gaps are negative numbers
        gap_open = gaps_params.gap_open
        gap_ext = gaps_params.gap_ext
        l_seq1 = len(seq_1)
        l_seq2 = len(seq_2)

        if(gap_ext == None): # case of linear gap applied
            # backtrack matrix
            Vp = self._init_pointer_table(l_seq1, l_seq2)
            # initialize V matrix -- holding the best alignment score for seq1[1:i] and seq2[1:j]
            V = np.zeros((l_seq1+1, l_seq2+1), dtype='float32')
            # initialize the transition time tables
            TR, TC = self._init_transition_table(seq_1, l_seq1, seq_2, l_seq2)
            if(align_type == 'global'):
                # fill the top row  with the linear gap
                V[0,1:] = gap_open * np.arange(1, l_seq2+1, dtype='float32')
                # fill the first left column with the linear gap
                V[1:,0] = gap_open * np.arange(1, l_seq1+1, dtype='float32')
               
            elif(align_type == 'semi-global'):
                # in this implementation, seq_2 is always considered the pattern sequence (i.e. pattern sequence should be the shorter sequence)
                # we penalize gaps in the pattern sequence (seq_2) but not in the long sequence (seq_1)
                # fill the top row with the linear gap
                V[0, 1:] = gap_open * np.arange(1, l_seq2+1, dtype='float32')    
            for i in range(1, l_seq1+1):
                for j in range(1, l_seq2+1):
                    # compute score
                    node_seq1 = seq_1[i-1]
                    node_seq2 = seq_2[j-1]
                    sim_score = node_seq1.compute_similarity(node_seq2, **kwargs)
                    max_s= -np.inf
                    if(TR[i-1, j-1] == TC[i-1, j-1] == (0,0,0)): # case of TR and TC are (0,0,0) at i-1 and j-1 position
                        temporal_penalty = node_seq1.compute_temporal_penalty(node_seq1.trans_time, 
                                                                              node_seq2.trans_time, 
                                                                              **kwargs)
                        max_s = V[i-1, j-1] + sim_score - temporal_penalty
                    else:
                        for pos in range(3):
                            if(TR[i-1, j-1][pos] or TC[i-1, j-1][pos]):
                                accum_tseq_1 = TR[i-1, j-1][pos] + node_seq1.trans_time
                                accum_tseq_2 = TC[i-1, j-1][pos] + node_seq2.trans_time
                                temporal_penalty = node_seq1.compute_temporal_penalty(accum_tseq_1, 
                                                                                      accum_tseq_2, 
                                                                                      **kwargs)
                                tmp = V[i-1, j-1] + sim_score - temporal_penalty
                                if(tmp>max_s):
                                    max_s = tmp
                        
                    case_1 = V[i, j-1] + gap_open # align node_j from seq_2 to a gap
                    case_2 = max_s # align node_i from seq_1 to node_j from seq_2
                    case_3 = V[i-1, j] + gap_open # align node_i from seq_1 to a gap 
                    cases_vec = [case_1, case_2, case_3]
#                     print('i,j = ', (i,j))
#                     print('cases_vec ', cases_vec)
                    self._fill_dynamic_and_pointer_tables(cases_vec, V, Vp, (i,j), main_dtable = True)
                    # updated the TR and TC tables
                    TR_l = [0, 0, 0]
                    TC_l = [0, 0, 0]                        
                    for indx, flag in enumerate(Vp[i,j]):
                        if(flag):
                            if(indx == 0):
                                TR_l[indx] = TR[i, j-1][-1]
                                TC_l[indx] = TC[i, j-1][0] + node_seq2.trans_time
                            elif(indx == 2):
                                TR_l[indx] = TR[i-1, j][-1] + node_seq1.trans_time
                                TC_l[indx] = TC[i-1, j][0] 
                    TR[i, j] = tuple(TR_l)
                    TC[i, j] = tuple(TC_l)
            print("TR ", TR)
            print("TC ", TC)
            print("V ", V)
            print("Vp ", Vp)    
            return(V, Vp, i, j)

        else: # case of affine gap applied
            raise(Exception("Affine gap is not supported yet!!..."))

    
    def _init_transition_table(self, seq_1, l_seq1, seq_2, l_seq2):
        # TR table
        TR = {(0,0): (0,0,0)}
        for i in range(1, l_seq1+1):
            trans_time = seq_1[i-1].trans_time
            TR[i, 0] = (0, 0, TR[i-1, 0][-1] + trans_time)
        for j in range(1, l_seq2+1):
            TR[0, j] = (0, 0, TR[0, j-1][-1])      
        
        # TC table
        TC = {(0,0): (0,0,0)}
        for j in range(1, l_seq2+1):
            trans_time = seq_2[j-1].trans_time
            TC[0, j] = (TC[0, j-1][0]+trans_time, 0, 0)     
        for i in range(1, l_seq1+1):
            TC[i, 0] = (TC[i-1, 0][0], 0, 0)
        return(TR, TC)

    
    def retrieve_alignments(self, *args, **kwargs):
        res = None
        if(len(args) == 4):
            res = self._retrieve_alignments_linear(*args, **kwargs)
        return(res)
        