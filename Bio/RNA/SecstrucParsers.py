#!/usr/bin/env python

from Secstruc import Secstruc, SecstrucError

#
# Parsers for common secondary structure formats
#

from Secstruc import Secstruc, PseudoknotSecstruc

def parse_vienna(data): 
    """Returns a (sequence, secstruc) tuple."""
    seq, secstruc = None, None
    for line in data:
        if line.startswith('>'): continue
        if not seq: seq = line.strip()
        elif not secstruc: 
            secstruc = Secstruc(line.strip())
            return seq, secstruc


def parse_bpseq(data):
    """Returns a list of base pairs parsed from a BPSeq format line iterator."""
    # remove header
    data = [d for d in data][1:] 
    seq = ""
    ss = ""
    for line in data:
        t = line.strip().split()
        if len(t)== 3:
            seq += t[1]
            index = int(t[0])-1
            pair = int(t[2])-1
            if pair == -1: ss += '.'
            elif pair > index: ss += '('
            elif pair < index: ss += ')'
            else:
                raise SecstrucError("Unable to parse this line in BPSEQ:%s"%line)
    return seq, PseudoknotSecstruc(ss)
        

def parse_ct(data):
    """Returns a list of base pairs parsed from a CT format line iterator."""
    # remove header
    data = [d for d in data][1:] 
    seq = ""
    ss = ""
    for line in data:
        t = line.strip().split()
        if len(t)== 6:
            seq += t[1]
            index = int(t[2])
            pair = int(t[4])-1
            if pair == -1: ss += '.'
            elif pair > index: ss += '('
            elif pair < index: ss += ')'
            else:
                raise PseudoknotSecstrucError("Unable to parse this line in BPSEQ:%s"%line)
    return seq, Secstruc(ss)
        
#TODO: implement
def write_vienna(seq, secstruc):
    pass
    
def write_bpseq(seq, secstruc):
    pass
    
def write_ct(seq, secstruc):
    pass
    

