
# coding: utf-8

# In[32]:

from Bio import Seq, SeqIO
codons = ['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG',           'TAG', 'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA', 'CAA', 'CGA',           'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA',           'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG', 'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC',           'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG']
class Heuristics: pass


# In[63]:

class Sequence(Heuristics):
    def __init__(self,filename):
        self.seqs = [seq.seq.upper() for seq in SeqIO.parse(filename,'fasta')]
        self.NumSeq = len(self.seqs)
        self.counts = {Letter:map(lambda seq: seq.count(Letter),self.seqs) for Letter in 'ATGC'}
        
    def get_seqs(self):return self.seqs
    
    def get_counts(self,Letter=None): 
        if Letter == None: return self.counts
        else: return self.counts[Letter]
    
    def get_length(self): return sum([ sum(self.counts[Letter]) for Letter in 'ATGC'])
    
    def get_lengths(self,i): return sum([self.counts[Letter][i] for Letter in 'ATGC']) 


# In[70]:

class GenomeSeq(Sequence):
    def __init__(self,filename):
        Sequence.__init__(self,filename)
        print 'processed',filename
        
    def get_GC(self): return (sum(self.counts['G'])+sum(self.counts['C']))/float(self.get_length())
    def get_GCs(self): return [(self.counts['G'][i]+self.counts['C'][i])/float(self.get_lengths(i)) for i in range(self.NumSeq)]
    
    
    def get_freq(self): return {Letter:float(sum(self.counts[Letter]))/self.length for Letter in 'ATGC'}
    
    def get_freqs(self,Letter=None):
        if Letter == None: return {Letter: [1.0*self.counts[Letter][i]/self.lengths[i] for i in range(self.NumSeq)]                                    for Letter in 'ATGC'}
        else: return [float(self.counts[Letter][i])/self.get_lengths(i) for i in range(self.NumSeq)]


# In[71]:

class CodingSeq(Sequence):
    def __init__(self,filename):
        print 'processing',filename
        Sequence.__init__(self,filename)
        self.seqs_filtered = [seq for seq in self.seqs if 
                             self.check_start(seq) and 
                             self.check_stop(seq) and 
                             self.check_length(seq) and 
                             self.check_shift(seq) and 
                             len(seq)>10]
        
        self.NumSeq_filtered = len(self.seqs_filtered)
        self.short_seqs = [seq[3:-3] for seq in self.seqs_filtered] #with head codon and tail codon removed for GC calc 
        self.short_lengths = map(len,self.short_seqs)
        self.short_counts = {Letter:map(lambda seq: seq.count(Letter),self.short_seqs) for Letter in 'ATGC'}
        self.GC = (sum(self.short_counts['G'])+sum(self.short_counts['C']))/float(sum(self.short_lengths))
        self.GCs = [(self.short_counts['G'][i]+self.short_counts['C'][i])/float(self.short_lengths[i])                     for i in range(self.NumSeq_filtered)]
        
        self.positional_seqs = [[seq[i::3] for seq in self.short_seqs] for i in range(3)]
        self.positional_counts = [{Letter:[self.positional_seqs[i][j].count(Letter)                                            for j in range(len(self.positional_seqs[i]))]                                      for Letter in 'ATGC'} for i in range(3)]
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    ##  constructor helper functions  
    def check_start(self,seq): return str(seq[:3])=='ATG'
    def check_stop(self,seq): return str(seq[-3:]) in ['TAA', 'TGA', 'TAG']
    def check_length(self,seq): return len(seq)%3 == 0
    def check_shift(self,seq):
        seq = str(seq[3:-3])
        for i in range(3,len(seq)-2,3): 
            if seq[i:i+3] in ['TAA', 'TGA', 'TAG']: return False
        return True
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    
    def get_NumCDS(self): return len(self.seqs_filtered)
    def get_GC(self): return self.GC
    def get_GCs(self): return self.GCs
    
    def get_freq(self):
        TotalLength = sum(self.short_lengths)
        return {Letter:float(sum(self.short_counts[Letter]))/TotalLength for Letter in 'ATGC'}
    
    def get_freqs(self,Letter=None):
        if Letter == None: return {Letter: [1.0*self.short_counts[Letter][i]/self.short_lengths[i]                                             for i in range(self.NumSeq_filtered)] for Letter in 'ATGC'}
        else: return [float(self.short_counts[Letter][i])/self.short_lengths[i] for i in range(self.NumSeq_filtered)]    
    
    def get_positional_freq(self):
        return [{Letter:1.0*sum(self.positional_counts[i][Letter])/(sum(map(len,self.positional_seqs[i])))                  for Letter in 'ATGC'} for i in range(3)]
    
    def get_codon_counts(self):
        result = {codon:0 for codon in codons}
        for seq in map(str,self.short_seqs):
            for i in range(0,len(seq)-2,3):
                result[seq[i:i+3]]+=1
        return result
    
    def get_codon_freqs(self):
        counts = self.get_codon_counts()
        total = sum(counts.values())
        return {codon:1.0*counts[codon]/total for codon in codons}

