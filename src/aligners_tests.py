'''
@author: ahmed allam <ahmed.allam@yale.edu>
'''
import numpy as np
from aligners import TemporalAligner, Aligner, GapsParams
from aligners_utilities import generate_char_seqnodes, custom_substitution_table

def test_tempaligner_1():
    """ this example comes from the original paper of temporal aligner `Syed et al. Temporal Needleman-Wunsch <http://ieeexplore.ieee.org/document/7344785/>`__
    """
    s1 = "ABCD"
    s2 = "AD" 
    print("s1: ", s1)
    print("s2: ", s2)
    uniform_t = 2
    ls1 = generate_char_seqnodes(s1, [uniform_t]*3)
    ls2 = generate_char_seqnodes(s2, [uniform_t*3])
    score_matrix = custom_substitution_table(set(s1), 1, -1.1)
    gaps_param = GapsParams(-0.5)
    print("score_matrix \n", score_matrix)
    print("gap_open: ", -0.5)
    for align_type in ('global', 'local'):
        print("{} alignment ...".format(align_type))
        aligner = TemporalAligner(align_type)
        res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix, T_penalty=0.25)
        all_paths = aligner.retrieve_alignments(*res, num_paths=np.inf)
        print("alignment paths:")
        print(all_paths)
        print("-"*50)
     
def test_tempaligner_2():
    """ variation on :func:`test_tempaligner_1()` function
    """
    s1 = "ABCD"
    s2 = "AD" 
    print("s1: ", s1)
    print("s2: ", s2)
    uniform_t = 2
    ls1 = generate_char_seqnodes(s1, [uniform_t]*3)
    ls2 = generate_char_seqnodes(s2, [uniform_t*4])
    score_matrix = custom_substitution_table(set(s1), 1, -1.1)
    gaps_param = GapsParams(-0.5)
    print("score_matrix \n", score_matrix)
    print("gap_open: ", -0.5)
    for align_type in ('global', 'local'):
        print("{} alignment ...".format(align_type))
        aligner = TemporalAligner(align_type)
        res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix, T_penalty=0.25)
        all_paths = aligner.retrieve_alignments(*res, num_paths=np.inf)
        print("alignment paths:")
        print(all_paths)
        print("-"*50)
        
     
def test_tempaligner_3():
    """ using non-temporal aligner to contrast it with :func:`test_tempaligner_2()` function
    """
    s1 = "ABCD"
    s2 = "AD" 
    print("s1: ", s1)
    print("s2: ", s2)
    uniform_t = 2
    ls1 = generate_char_seqnodes(s1, [uniform_t]*3)
    ls2 = generate_char_seqnodes(s2, [uniform_t*4])
    score_matrix = custom_substitution_table(set(s1), 1, -1.1)
    gaps_param = GapsParams(-0.5)
    print("score_matrix \n", score_matrix)
    print("gap_open: ", -0.5)
    for align_type in ('global', 'local'):
        print("{} alignment ...".format(align_type))
        aligner = Aligner(align_type)
        res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
        all_paths = aligner.retrieve_alignments(*res, num_paths=np.inf)
        print("alignment paths:")
        print(all_paths)
        print("-"*50)
""" 
    The following tests are aimed at checking the correctness of the implementation against
    cases that shown to be challenging for other alignment implementations.
    These tests come from the following paper:                   
        `Are all global alignment algorithms and implementations correct? <http://www.biorxiv.org/content/biorxiv/early/2015/11/12/031500.full.pdf>`__

"""

def test_aligner_1():
    s1 = 'AAAGGG' 
    s2 = 'TTAAAAGGGGTT'
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1)
    ls2 = generate_char_seqnodes(s2)
    aligner = Aligner('global')
    score_matrix = custom_substitution_table(set(s1).union(set(s2)), 0, -1)
    gaps_param = GapsParams(-5, gap_ext=-1)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
    all_paths = aligner.retrieve_alignments(*res)
    print("alignment_type: global")
    print("score_matrix \n", score_matrix)
    print("gap_open:{} , gap_ext:{}".format(-5, -1))
    print("alignment paths:")
    print(all_paths)
    
def test_aligner_2():
    s1 = 'CGCCTTAC'
    s2 = 'AAATTTGC' 
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1)
    ls2 = generate_char_seqnodes(s2)
    aligner = Aligner('global')
    score_matrix = custom_substitution_table(set(s1).union(set(s2)), 10, -30)
    gaps_param = GapsParams(-40, gap_ext=-1)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
    all_paths = aligner.retrieve_alignments(*res)
    print("alignment_type: global")
    print("score_matrix \n", score_matrix)
    print("gap_open:{} , gap_ext:{}".format(-40, -1))
    print("alignment paths:")
    print(all_paths)

def test_aligner_3():
    s1 = 'TAAATTTGC'
    s2 = 'TCGCCTTAC' 
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1)
    ls2 = generate_char_seqnodes(s2)
    aligner = Aligner('global')
    score_matrix = custom_substitution_table(set(s1).union(set(s2)), 10, -30)
    gaps_param = GapsParams(-40, gap_ext=-1)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
    all_paths = aligner.retrieve_alignments(*res)
    print("alignment_type: global")
    print("score_matrix \n", score_matrix)
    print("gap_open:{} , gap_ext:{}".format(-40, -1))
    print("alignment paths:")
    print(all_paths)

def test_aligner_4():
    s1 = 'AGAT'
    s2 = 'CTCT' 
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1)
    ls2 = generate_char_seqnodes(s2)
    aligner = Aligner('global')
    score_matrix = custom_substitution_table(set(s1).union(set(s2)), 10, -30)
    gaps_param = GapsParams(-25, gap_ext=-1)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
    all_paths = aligner.retrieve_alignments(*res)
    print("alignment_type: global")
    print("score_matrix \n", score_matrix)
    print("gap_open:{} , gap_ext:{}".format(-25, -1))
    print("alignment paths:")
    print(all_paths)
    print("-"*50)
    gaps_param = GapsParams(-30, gap_ext=-1)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
    all_paths = aligner.retrieve_alignments(*res)
    print("gap_open:{} , gap_ext:{}".format(-30, -1))
    print("alignment paths:")
    print(all_paths)
    
def test_aligner_5():
    s1 = 'CTGTTCTCCCTACATCAAATGTCTATCCCCGCACCAAGTGGAGATTCCATGAGGATGAGG'
    s2 = 'CTGCCAAGTGGAGATTCCATGAGGATGAGG'
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1)
    ls2 = generate_char_seqnodes(s2)
    aligner = Aligner('global')
    score_matrix = custom_substitution_table(set(s1).union(set(s2)), 100, 0)
    gaps_param = GapsParams(-10, gap_ext=-2)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
    all_paths = aligner.retrieve_alignments(*res)
    print("alignment_type: global")
    print("score_matrix \n", score_matrix)
    print("gap_open:{} , gap_ext:{}".format(-10, -2))
    print("alignment paths:")
    print(all_paths)

def test_aligner_6():
    """ testing semi-global alignment 
        example coming from http://www.comp.nus.edu.sg/~ksung/algo_in_bioinfo/slides/Ch2_sequence_similarity.pdf
    """
    #'semi-global':('ATCCGAACATCCAATCGAAGC', 'AGCATGCAAT'),
    align_options = {'semi-global':('TCAACGATCACCGCA', 'ACCTCACGATCCGA'),'end-gap-free':('TCAACGATCACCGCA', 'ACCTCACGATCCGA')}
    gap_open = -1
    match = 2
    mismatch = -1
    for align_type, seq_tup in align_options.items():
        aligner = Aligner(align_type)
        s1, s2 = seq_tup
        print("s1: ", s1)
        print("s2: ", s2)
        ls1 = generate_char_seqnodes(s1)
        ls2 = generate_char_seqnodes(s2)
        score_matrix = custom_substitution_table(set(s1).union(set(s2)), match, mismatch)
        gaps_param = GapsParams(gap_open, gap_ext=None)
        res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
        all_paths = aligner.retrieve_alignments(*res)
        print("alignment_type: ", align_type)
        print("score_matrix \n", score_matrix)
        print("gap_open:{} , gap_ext:{}".format(gap_open, None))
        print("alignment paths:")
        print(all_paths)
        for path in all_paths[-1]:
            o1=[]
            o2=[]
            for alignment in path:
                o1.append(alignment[0])
                o2.append(alignment[-1])
            print("".join(o1))
            print("".join(o2))
            print("~"*40)
        print("-"*40)
    
def test_aligner_7():
    """ testing local alignment 
        example coming from http://www.comp.nus.edu.sg/~ksung/algo_in_bioinfo/slides/Ch2_sequence_similarity.pdf
    """
    #'semi-global':('ATCCGAACATCCAATCGAAGC', 'AGCATGCAAT'),
    gap_open = -1
    match = 2
    mismatch = -1
    align_type = 'local'
    aligner = TemporalAligner(align_type)
    s1 = 'ACAATCG'
    s2 =  'CTCATGC'
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1, [0]*(len(s1)-1))
    ls2 = generate_char_seqnodes(s2, [0]*(len(s2)-1))
    score_matrix = custom_substitution_table(set(s1).union(set(s2)), match, mismatch)
    gaps_param = GapsParams(gap_open, gap_ext=None)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_matrix)
    print("V: ", res[0])
    print("V: ", res[1])
    all_paths = aligner.retrieve_alignments(*res, num_paths=np.inf)
    print("alignment_type: ", align_type)
    print("score_matrix \n", score_matrix)
    print("gap_open:{} , gap_ext:{}".format(gap_open, None))
    print("alignment paths:")
    print(all_paths)
    for path in all_paths[-1]:
        o1=[]
        o2=[]
        for alignment in path:
            o1.append(alignment[0])
            o2.append(alignment[-1])
        print("".join(o1))
        print("".join(o2))
        print("~"*40)
    print("-"*40)
    
def scoremat_example1():
    # define score mat
    score_mat = {}
    score_mat['A', 'B'] = 1/2
    score_mat['A', 'C'] = 1/4
    score_mat['A', 'D'] = 0
    score_mat['B', 'C'] = 1/2
    score_mat['B', 'D'] = 1/4
    score_mat['C', 'D'] = 1/2
    tmp = {}
    for elem, score in score_mat.items():
        a, b = elem
        tmp[b,a] = score
    score_mat.update(tmp)
    del tmp
    for char in ('A', 'B', 'C', 'D'):
        score_mat[char, char]= 1
    return(score_mat)

def scoremat_example2():
    # define score mat
    score_mat = {}
    score_mat['A', 'B'] = 1/2
    score_mat['A', 'C'] = 1/4
    score_mat['A', 'D'] = 0
    score_mat['B', 'C'] = 1/2
    score_mat['B', 'D'] = 1/4
    score_mat['C', 'D'] = 1/2
    tmp = {}
    for elem, score in score_mat.items():
        a, b = elem
        tmp[b,a] = score
    score_mat.update(tmp)
    del tmp
    for char in ('A', 'B', 'C', 'D'):
        score_mat[char, char]= 1
    for char_i in ('A', 'B', 'C', 'D'):
        for char_j in ('E', 'F', 'G', 'H', 'I', 'J', 'K'):
            score_mat[char_i, char_j] = 0
            score_mat[char_j, char_i] = 0
    return(score_mat)

def test_alignment_paramchoice_1(s1, s2, gap_open):

    score_mat = scoremat_example1()
        
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1)
    ls2 = generate_char_seqnodes(s2)
    aligner = Aligner('global')
    gaps_param = GapsParams(gap_open)
    res = aligner.align(ls1, ls2, gaps_param, score_matrix = score_mat)
    all_paths = aligner.retrieve_alignments(*res)
    print("alignment_type: global")
    print("score_matrix \n", score_mat)
    print("gap_open: ", gap_open)
    print("alignment paths:")
    print(all_paths)

def test_alignment_paramchoice_2(s1, s2, gap_open, transtime_s1=[], transtime_s2=[], align_type='global'):

    score_mat = scoremat_example2()
        
    print("s1: ", s1)
    print("s2: ", s2)
    ls1 = generate_char_seqnodes(s1, transtime_s1)
    ls2 = generate_char_seqnodes(s2, transtime_s2)
    aligner = TemporalAligner(align_type)
    gaps_param = GapsParams(gap_open)
    res = aligner.align(ls2, ls1, gaps_param, score_matrix = score_mat)
    score, all_paths = aligner.retrieve_alignments(*res)
    print("alignment_type: ", align_type)
    print("score_matrix \n", score_mat)
    print("gap_open: ", gap_open)
    print("original_score: ", score)
    len_ls1=len(ls1)
    len_ls2=len(ls2)
    if(len_ls1 < len_ls2):
        max_score = len_ls1 + gap_open*(len_ls2-len_ls1)
    else:
        max_score = len_ls2 + gap_open*(len_ls1-len_ls2)
    print("max score: ", max_score)
    min_score = gap_open*(len_ls1 + len_ls2)
    print("min score: ", min_score)
    transformed_score = (score - min_score)/(max_score - min_score)
    print("transformed score: ", transformed_score)
    print("transformed score 2: ", score/max_score)
    print("alignment paths:")
    print(all_paths)
    
    
def alignment_paramchoice_1():
    s1='A'
    s2='B'
    for gap_open in (0.24, 0.25, 0.26):
        test_alignment_paramchoice_1(s1, s2, gap_open)
        print("-"*40)
        
def alignment_paramchoice_2():
    s1='BC'
    s2='AD'
    for gap_open in (0.24, 0.25, 0.26):
        test_alignment_paramchoice_1(s1, s2, gap_open)
        print("-"*40)
             
def alignment_paramchoice_3():
    s1='BC'
    s2='AD'
    for gap_open in (0.24, 0.25, 0.26):
        test_alignment_paramchoice_1(s1, s2, gap_open)
        print("-"*40)

def alignment_paramchoice_4():
    s1='ABB'
    s2='DCC'
    for gap_open in (0.1, 0.24, 0.26):
        test_alignment_paramchoice_1(s1, s2, gap_open)
        print("-"*40)

def alignment_paramchoice_5():
    s1='AB'
    s2='AB'
    transtime_s1 = [100]
    transtime_s2 = [10]
    print("using temporal aligner")
    gap_open = 0.24
    test_alignment_paramchoice_2(s1, s2, gap_open, transtime_s1, transtime_s2)
    print("-"*40)                        

    print("using non-temporal aligner")
    gap_open = 0.24
    test_alignment_paramchoice_1(s1, s2, gap_open)
    print("-"*40)                 
        
        
def alignment_paramchoice_6():
    """ using temporal penalty as currently implemented, the t_p value specifies the maximum penalty to subtract from 
        any two matching events. Hence, if we have two events that have perfect match score 1, but have huge difference
        in the transition time, the similarity of these two events will be equal to (max_score - t_p)
    """
    s1='AB'
    s2='AB'
    transtime_s1 = [100]
    transtime_s2 = [10]
    for gap_open in (0.24, 0.3, 0.4):
        test_alignment_paramchoice_2(s1, s2, gap_open, transtime_s1, transtime_s2)
    print("-"*40)                        

    s1='AB'
    s2='AC'
    transtime_s1 = [0]
    transtime_s2 = [0]
    gap_open = 0.24
    test_alignment_paramchoice_2(s1, s2, gap_open, transtime_s1, transtime_s2)
    print("-"*40)  
    
def alignment_paramchoice_7(align_type):
    """ using temporal penalty as currently implemented, the t_p value specifies the maximum penalty to subtract from 
        any two matching events. Hence, if we have two events that have perfect match score 1, but have huge difference
        in the transition time, the similarity of these two events will be max_score - t_p
    """
    char_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']                    
    s1='AB'
    gap_open = 0.24
    transtime_s1 = [0]
    for i in range(3, 11):
        transtime_s2 = [0]*(i-1)
        s2="".join(char_list[:i])
        test_alignment_paramchoice_2(s1, s2, gap_open, transtime_s1, transtime_s2, align_type)
        print("-"*40)
        
def alignment_paramspace(align_type):
    pass



