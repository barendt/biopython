
import re

class SecstrucError(Exception): pass
class UnknownElementError(Exception): pass

class Secstruc:
    """
    Secondary structure string with indices for each position.
    Can handle pieces of secondary structure which are interrupted.
    """
    def __init__(self,secstr,indices=None):
        self.Secstruc = secstr
        if indices == None:
            indices = range(len(self))
        if len(self.Secstruc) != len(indices):
            raise SecstrucError("Cannot create Secstruc object (%s %s)"%(secstr,str(indices)))
        self.indices = indices

    def find_overhang5(self):
        """Returns overhang on the 5' end as Secstruc objects."""
        n = len(self)
        for i in range(2,n):
            overhang = self[:i]
            if overhang.get_type() == "5'-overhang":
                return overhang

    def find_overhang3(self):
        """Returns overhang on the 3' end as Secstruc objects."""
        n = len(self)
        for i in range(n-1,0,-1):
            overhang = self[i:]
            if overhang.get_type() == "3'-overhang":
                return overhang

    def find_helices(self):
        # return self.extract_elements('helix')
        helices = []
        begin = 0
        end = 0
        length = 0

        for i,j in self.get_base_pair_indices():
            # check if old helix continues
            if i == begin+length and j==end-length:
                length += 1
            else:
                # save old helix
                if length >0:
                    helices.append(self[begin:begin+length] + self[end-length+1:end+1])
                # new helix starts
                begin = i
                end = j
                length = 1
        if length > 0:
            helices.append(self[begin:begin+length] + self[end-length+1:end+1])
        return helices

    def check_bulge(self, i,j):
        if i>len(self)-3 or j==0: return None, None
        ele = self[i] + self[i+1] + self[i+2] + self[j-1] + self[j]
        if ele.get_type() == 'bulge':
            return ele, i+1
        ele = self[i] + self[i+1] + self[j-2] + self[j-1] + self[j]
        if ele.get_type() == 'bulge':
            return ele, j-1
        return None,None

    def get_base_pair_indices(self):
        """generates (i,j) tuples iterating through base pairs."""
        n = len(self)
        i = 0
        while i<n-1:
            # find next open pair
            while i<n-1 and self.Secstruc[i]!='(': i += 1
            # find corresponding pair
            j = i+1
            bps = 1
            while j<n and bps>0:
                if self.Secstruc[j]=='(': bps += 1
                if self.Secstruc[j]==')': bps -= 1
                if bps == 0: yield i,j
                j += 1
            i += 1
 
    def find_bulges(self):        
        bulges = []
        for i,j in self.get_base_pair_indices():           
            bulge, ind = self.check_bulge(i,j)
            if bulge:
                bulges.append(bulge)
        return bulges

    def find_loops(self):        
        return self.find_elements('loop')
    
    def find_junctions(self):
        """Returns a list of all junctions.
        Any region between two paired bases with more than 1 nested base pair
        Inside is a junction.
        """
        junctions = []
        for i,j in self.get_base_pair_indices():
            # check if i,j enclose a junction.
            nested_bps = 0
            skip_helix = 0
            junction = self[i]
            k = i+1
            while k < j:
                if self.Secstruc[k] == ')':
                    skip_helix -= 1
                    if skip_helix == 0:
                        nested_bps += 1
                if skip_helix == 0:
                    junction += self[k]
                if self.Secstruc[k] == '(':
                    skip_helix += 1
                k += 1
            junction += self[j]
            if nested_bps >= 2:
                junctions.append(junction)
        return junctions
    
    def find_elements(self, element_type):
        """Finds elements of a given type."""
        result = []
        n = len(self)
        i = 0
        while i < n-1:
            j = n
            while j > i:
                ele = self[i:j]
                if ele.get_type() == element_type:
                    result.append(ele)
                j -= 1
            i += 1
        return result
                            

    def find_Secstruc_elements(self):
        """
        Takes a secondary structure, and returns a list of
        component secondary structure elements.
        """
        result = []
        overhang = self.find_overhang5()
        if overhang: result.append(overhang)
        overhang = self.find_overhang3()
        if overhang: result.append(overhang)
        loops = self.find_loops()
        helices = self.find_helices()
        bulges = self.find_bulges()
        junctions = self.find_junctions()
        result += loops + bulges + helices + junctions
        return result

    def get_type(self):
        """Returns a type classification of this secondary structure as a string,
        if it is some basic type."""
        if re.search('^\(\(+\)\)+$',self.Secstruc):
            # check both stems have equal length
            half = len(self.Secstruc)/2
            if len(self.Secstruc)%2 == 0 \
               and re.search('^\(+$',self.Secstruc[:half]) \
               and re.search('^\)+$',self.Secstruc[half:]):
                return 'helix'
            else: return None
        if re.search('^\.+\($',self.Secstruc): return "5'-overhang"
        if re.search('^\)\.+$',self.Secstruc): return "3'-overhang"
        if re.search('^\(\.\.+\)$',self.Secstruc): return "loop"
        if re.search('^\(\(\)\.+\)|\(\.+\(\)\)$',self.Secstruc): return "bulge"
        if re.search('^\((.*\(\)\.*)+\)$',self.Secstruc): return 'junction'
        return None


    def __eq__(self, other):
        if not isinstance(other,Secstruc):
            other = Secstruc(other)
        if self.Secstruc == other.Secstruc \
           and self.indices == other.indices:
            return True
        
    def __add__(self, other):
        newsec = self.Secstruc + other.Secstruc
        newind = self.indices + other.indices
        return Secstruc(newsec, newind)

    def __getitem__(self,index):
        if type(index) == slice:
            newsec = self.Secstruc[index.start:index.stop]
            newind = self.indices[index.start:index.stop]
            return Secstruc(newsec, newind)
        else:
            return Secstruc(self.Secstruc[index],[self.indices[index]])

    def __repr__(self):
        return self.Secstruc + ';' + str(self.indices)

    def __len__(self):
        return len(self.Secstruc)

    def __str__(self):
        return self.Secstruc

def contains_pseudoknot(bplist):
    """Checks if a list of base pairs contains a pseudoknot. Returns boolean"""
    for a1, a2 in bplist:
        for b1, b2 in bplist:
            if (a1 < b1 < a2 < b2) or (a2 < b1 < a1 < b2) or (a2 < b2 < a1 < b1):
                return True
    return False
    
    
def get_pseudoknot_pairs(bplist):
    """Returns a list of base pairs that participate in pseudoknots."""
    result = []
    for a1, a2 in bplist:
        for b1, b2 in bplist:
            if (a1 < b1 < a2 < b2) or (a2 < b1 < a1 < b2) or (a2 < b2 < a1 < b1) \
            or (b1 < a1 < b2 < a2) :
                if (a1, a2) not in result:
                    result.append((a1, a2))
    return result
                
def make_dotbracket_from_bplist(length, bplist):    
    """
    Takes a structure as a list of basepairs
     [(2,7),(4,5),(10,15),(17,20)]
     and converts it to vienna format: 
     "((...).(((.((....))).)))"
    """
    if contains_pseudoknot(bplist):
        raise PseudoknotError("Base pairs contain a pseudoknot: %s"%str(bplist))
    Secstruc = ['.']*length
    for bp in bplist:
        Secstruc[bp[0]-1] = '(' 
        Secstruc[bp[1]-1] = ')'    
    return ''.join(Secstruc)
    
def make_bplist_from_dotbracket(Secstruc):
    """Reads a dot-bracket secondary structure and returns a list of base pair indices"""
    result = []
    for bp in PseudoknotSecstruc(Secstruc).get_base_pair_indices():
        result.append((bp[0]+1, bp[1]+1))
    return result

def has_overlapping_basepairs(bplist):
    indices = {}
    for bp in bplist:
        if indices.has_key(bp[0]): return bp[0]
        if indices.has_key(bp[1]): return bp[1]
        indices[bp[0]]=True
        indices[bp[1]]=True
    return False
    
"""
This modules contains parses for output formats of secondary structure prediction programs
for single sequences or groups of similar sequences (comparative prediction).
They all can be parsed into lists of basepairs, 
and lists of base pairs can be converted to dot-bracket sequences

Example of structure without pseudoknots, takes one line:
    ((..(..)..).((...).(((...)))).)
Example of structure with pseudoknots, takes more than one line:
    ((..(..)..)..((...).(((...)))).)
    .......((........).)............
    ...........((..........)).......
"""
 
class PseudoknotError(Exception): pass

class PseudoknotSecstruc(Secstruc):
    symbols = {
        '(':')', '[':']', '{':'}', 
        'a':'A', 'b':'B', 'c':'C', 'd':'D', 'e':'E', 'f':'F', 
        'g':'G', 'h':'H', 'i':'I', 'j':'J', 'k':'K', 'l':'L', 
        'm':'M', 'n':'N', 'o':'O', 'p':'P', 'q':'Q', 'r':'R', 
        's':'S', 't':'T', 'u':'U', 'v':'V', 'w':'W', 'x':'X', 'y':'Y', 'z':'Z', 
        }
    open_symbols = symbols.keys()
    
    def get_base_pair_indices(self):
        """generates (i,j) tuples iterating through base pairs."""
        n = len(self)
        i = 0
        while i<n-1:
            # find next open pair
            while i<n-1 and self.Secstruc[i] not in self.open_symbols: i += 1
            # find corresponding pair
            open_symbol = self.Secstruc[i]
            close_symbol = self.symbols.get(open_symbol)
            j = i+1
            bps = 1
            while j<n and bps>0:
                if self.Secstruc[j]==open_symbol: bps += 1
                if self.Secstruc[j]==close_symbol: bps -= 1
                if bps == 0: yield i,j
                j += 1
            i += 1
            
    def get_pseudoknot_indices(self):
        """Returns a list of base pair index tuples participating in pseudoknots."""
        bplist = [bp for bp in self.get_base_pair_indices()]
        pknot = get_pseudoknot_pairs(bplist)
        return pknot

