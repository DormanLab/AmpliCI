import re
import numpy as np


"""
This is the part of importing data from the fastq file.

Read data in fastq form from the file. The data should be in the form of the example showed below without blank line:
------------------------------------------
@SRR2241783.56.2_AAACTTTGA_1
GTCGTAATTTCTAGGTCCCCTCCTGTTGCATAGAG....
+
8CCFGFGGGGGGGGGCFGGGEGFFGGDFD@FCEF......
-----------------------------------------------
Construct three lists containing the name of the reads,
the nucleotide sequences of the reads and the quality sequences of the reads.
Use class ndarray in numpy package to store all nucleotide sequences and each quality sequences.
In each nucleotide sequence, use 0,1,2,3 to represent A,C,T,G.
"""

class Error(Exception):
    pass

class Line(str):
    """A line of text with associated filename and line number."""
    def error(self, message):
        """Return an error relating to this line."""
        return Error("{0}({1}): {2}\n{3}"
                     .format(self.filename, self.lineno, message, self))

class Lines(object):
    """Lines(filename, iterator) wraps 'iterator' so that it yields Line
    objects, with line numbers starting from 0. 'filename' is used in
    error messages.
    """
    def __init__(self, filename, iterator):
        self.filename = filename
        self.lines = enumerate(iterator, start=0)

    def __iter__(self):
        return self

    def __next__(self):
        lineno, s = next(self.lines)
        line = Line(s)
        line.filename = self.filename
        line.lineno = lineno
        return line

def read_fastq(filename, iterator):
    """Read FASTQ data from 'iterator' (which may be a file object or any
    other iterator that yields strings) and generate tuples (sequence
    name, sequence data, quality data). 'filename' is used in error
    messages."""
    at_seqname_re = re.compile(r'@(.+)$')
    sequence_re = re.compile(r'[!-*,-~]*$')
    plus_seqname_re = re.compile(r'\+(.*)$')
    quality_re = re.compile(r'[!-~]*$')

    lines = Lines(filename, iterator)
    for line in lines:
        # First line of block is @<seqname>.
        m = at_seqname_re.match(line)
        if not m:
            raise line.error("Expected @<seqname> but found:")
        seqname = m.group(1)
        try:
            # One line of sequence data.
            line = next(lines)
            m = sequence_re.match(line)
            sequence=m.group(0)
            if not sequence:
                raise line.error("Expected <sequence> but found:")

            # The line following the sequence data consists of a plus sign
            line = next(lines)
            m = plus_seqname_re.match(line)
            if not m:
                raise line.error("Expected +[<seqname>] but found:")

            # One line quality data, containing the same number of characters as the sequence data.
            line = next(lines)
            n_s = len(sequence)
            m = quality_re.match(line)
            if not m:
                raise line.error("Expected <quality> but found:")
            n_q = len(m.group(0))
            if n_s != n_q:
                raise line.error("<quality> is not equal to <sequence>:")
            quality = np.fromstring(m.group(0), dtype=np.int8) - 33
            yield seqname, ''.join(sequence), quality

        except StopIteration:
            raise line.error("End of input before sequence was complete:")


def read_fasta(filename, iterator):
    """Read FASTA data from 'iterator' (which may be a file object or any
    other iterator that yields strings) and generate tuples (sequence
    name, sequence data). 'filename' is used in error
    messages."""
    at_seqname_re = re.compile(r'>(.+)$')
    sequence_re = re.compile(r'[!-*,-~]*$')

    lines = Lines(filename, iterator)
    for line in lines:
        # First line of block is @<seqname>.
        m = at_seqname_re.match(line)
        if not m:
            raise line.error("Expected @<seqname> but found:")
        seqname = m.group(1)
        try:
            # One line of sequence data.
            line = next(lines)
            m = sequence_re.match(line)
            sequence=m.group(0)
            if not sequence:
                raise line.error("Expected <sequence> but found:")

            yield seqname, ''.join(sequence)

        except StopIteration:
            raise line.error("End of input before sequence was complete:")

def change_sequence_format(sequence):
    complement = {'A': '0', 'C': '1', 'T': '2', 'G':'3'}
    new_sequence=' '.join(complement.get(base, base) for base in sequence)
    return np.fromstring(new_sequence, dtype=np.int8, sep=' ')

def change_to_sequence_format(sequence_of_number):
    omega=['A','C','T','G']
    sequence = [omega[e] for e in sequence_of_number]
    sequence=''.join(sequence)
    return sequence

def print_fastq(f,name, seq, qual):
   print('@'+name,file = f)
   print(seq,file = f)
   print('+',file = f)
   print(qual,file = f)
        