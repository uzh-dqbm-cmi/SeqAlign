'''
@author: ahmed allam <ahmed.allam@yale.edu>
'''

import numpy as np
from .aligners import CharSeqNode, SeqNode


def generate_char_seqnodes(s, t = []):
    """ generate sequence nodes from character elements in a sequence
    
       The generated nodes are instances of :class:`CharSeqNode`
       
       Args:
           s: list, sequence of character elements
        
       Keyword args:
           t: list, (empty default) representing the transition time between elements in the sequence.
              We can specify the transition time between elements in case we want to use :class:`TemporalAligner`
    
    """
    s_l = []
    if(t):
        trans_time = [0] + t[:]
        for indx, char in enumerate(s):
            s_l.append(CharSeqNode(char, trans_time = trans_time[indx]))
    else:
        for char in s:
            s_l.append(CharSeqNode(char))
    return(s_l)

def custom_substitution_table(alphabet, match, mismatch):
    table = {}
    for a in alphabet:
        for b in alphabet:
            if(a == b):
                table[a,b] = match
            else:
                table[a,b] = mismatch
    return(table)

# def generate_p_seqnodes(num_events, fv_size):
#     """ generate synthetic sequence nodes
#     
#        The generated nodes are instances of :class:`SeqNode`
#        
#        Args:
#            num_events: number of events in the sequence
#            fv_size: the size of the feature vector representing the event
#     """
#     # generate time of events raw
#     t = np.random.rand(num_events)
#     t= np.round(t*100)
#     t = np.sort(t)
#     # get the transition time between time events
#     t = np.diff(t)
#     t = np.insert(t, 0, 0) # insert 0 at position 0 since the first element has no transition time
#     # generate feature vector representation
#     seq_l = []
#     for i in range(num_events):
#         label = "event_{}".format(i)
#         fv = np.random.randint(0,2,fv_size)
#         seq_l.append(SeqNode(label, fv, t[i]))
#     return(seq_l)

def generate_p_seqnodes(num_events, fv_size):
    """ generate synthetic sequence nodes
    
       The generated nodes are instances of :class:`SeqNode`
       
       Args:
           num_events: number of events in the sequence
           fv_size: the size of the feature vector representing the event
    """
    # generate time of events raw for the 12 months period
    t = np.random.randint(1, 13, num_events)
    t = np.sort(t)
    # generate feature vector representation
    seq_l = []
    for i in range(num_events):
        label = "event_{}".format(i)
        fv = np.random.randint(0,2,fv_size)
        seq_l.append(SeqNode(label, fv, t[i]))
    return(seq_l)