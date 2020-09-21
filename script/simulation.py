import fastq
import numpy as np
from scipy.stats import powerlaw
from scipy.stats import multinomial
from scipy.stats import binom
import argparse

parser = argparse.ArgumentParser(description='Simulate umi tagged sequence data')
parser.add_argument('-b', required=True, type=str, help="input barcodes fasta file")
parser.add_argument('-i', required=True, type=str, help="input reads fasta file")
parser.add_argument('-o', required=True, type=str, help="output file base name")
parser.add_argument('-k', required=False, default= 25, type=int, help="number of haplotypes")
parser.add_argument('-c', required=False, default= 7, type=int, help="number of PCR cycles")
parser.add_argument('-e', required=False, default= 0.9, type= float, help="PCR efficiency")
parser.add_argument('-n', required=False, default= 400, type= int, help="number of simulated moleculars")
parser.add_argument('-a', required=False, default= 0.5, type= float, help="parameter for power-law distribution")

params = parser.parse_args()

barcode_fasta = params.b
reads_fasta = params.i
output_file = params.o
K = params.k   ## number of haplotypes
n_cycle = params.c  # PCR amplification
p_amp = params.e
num_uniq = params.n   ## number of sequences before amplification
alpha = params.a ## controling the abundance distribution 

class Error(Exception):
    pass

### error rates from Table 2 in the paper https://link.springer.com/article/10.1186/1471-2105-12-38/tables/2
error_prob = np.array([[0.9995,7.2e-6,7.7e-6,5.1e-4],[1.1e-5,0.9996,4.1e-4,2.1e-6],[9.0e-6,5.7e-4,0.9994,1.4e-5],[3.5e-4,3.2e-6,2.1e-5,0.9996]])
error_prob[0,0] = 1.0 - error_prob[0,1] - error_prob[0,2] - error_prob[0,3]
error_prob[1,1] = 1.0 - error_prob[1,0] - error_prob[1,2] - error_prob[1,3]
error_prob[2,2] = 1.0 - error_prob[2,0] - error_prob[2,1] - error_prob[2,3]
error_prob[3,3] = 1.0 - error_prob[3,0] - error_prob[3,1] - error_prob[3,2]


def change_sequence_format(sequence):
    complement = {'A': '0', 'C': '1', 'T': '2', 'G':'3'}
    new_sequence=' '.join(complement.get(base, base) for base in sequence)
    return np.fromstring(new_sequence, dtype=np.int8, sep=' ')

def change_to_sequence_format(sequence_of_number):
    omega=['A','C','T','G']
    sequence = [omega[e] for e in sequence_of_number]
    sequence=''.join(sequence)
    return sequence

### [TODO] add codes to simulate sequence. 
def simu_seq(hap):
    sequence= []
    nuc_list = change_sequence_format(hap)
    l = len(hap)
    for i in range(0,l):
        p_nu = error_prob[nuc_list[i]]
        #print(p_nu,nuc_list[i])
        nu = np.random.choice(4, size = 1, p = p_nu)[0]
        #print(nu)
        sequence.append(nu)
    return(change_to_sequence_format(sequence))   


def print_fasta(f,name, seq):
   print('>'+name,file = f)
   print(seq,file = f)


def sto_amplif(seq_set, n_j, p_amp):
    n_j_1 = n_j
    for i in range(0,n_j):
        exist = binom.rvs(1,p_amp,size = 1)[0]
        #print(exist)
        if exist == 1:
            n_j_1 += 1   ### simulate a seq and append it to the end of the list seq_set
            seq_set.append(simu_seq(seq_set[i]))
    return(n_j_1)

## n( j + 1) = n( j ) + B(n( j ), Pamp)
## From the paper: Sources of PCR-induced distortions in high-throughput sequencing data sets
def pcr_copy(true_seq,p_amp,n_cycle,ori_abun):   
    abun = np.zeros((n_cycle+1),dtype=np.int32)
    abun[0] = int(ori_abun)
    seq_set = [true_seq]
    for i in range(0,n_cycle):
        # print(abun[i])
        abun[i+1] = sto_amplif(seq_set, abun[i],p_amp)
    return(abun[n_cycle],seq_set)

def simulate_fasta():

    ## Maybe using a powerlaw distritbuion to fit the abundance distribution
    abun_prob = powerlaw.rvs(alpha, size=K)    
    abun_prob = abun_prob / sum(abun_prob)
    print(abun_prob)
    abundance = multinomial.rvs(num_uniq, abun_prob)   ## multivariate normal to simulate count of the original unique sequences
    print(abundance)     ### abundance of true unique sequence

    ### use real barcodes and fasta for simulation
    real_barcodes=[]
    try:
        barcodefile = open(barcode_fasta, 'r')
        barcodes = fastq.read_fasta(barcode_fasta,barcodefile)
        for sequence in barcodes:
            real_barcodes.append(sequence[1])
            nb = len(real_barcodes)
        barcodefile.close()
    except IOError:
        raise Error("No barcodes fasta file")
    
    ### reads
    haplotypes = []
    try:
        readsfile = open(reads_fasta,'r')
        reads = fastq.read_fasta(reads_fasta,readsfile)
        for sequence in reads:
            haplotypes.append(sequence[1])
            nh = len(haplotypes)
        readsfile.close()
    except IOError:
        raise Error("No haplotype fasta file") 

    ### randomly select real_barcodes and haplotypes
    real_barcodes = np.random.choice(real_barcodes,num_uniq,replace = False)
    haplotypes = np.random.choice(haplotypes,K,replace = False)


    list_seq = []
    list_name = []
    list_seq_unamplified = []
    list_true_hap = []  ## indexed
    list_true_UMI = []  ## indexed
    count = 0   ## count number of unique seqeunces
    print(error_prob)

    for k in range(0,K):
        for i in range(0,abundance[k]):
            barcode_id = count    ## may change it to randomly chosen a barcode
            true_seq = real_barcodes[barcode_id]+ haplotypes[k]
            list_seq_unamplified.append(true_seq)  
            (abun, seq_set) = pcr_copy(true_seq, p_amp,n_cycle,1)
            # print(amp)
            list_seq = list_seq + seq_set
            for a in range(0,abun):
                list_true_hap.append(k)
                list_true_UMI.append(barcode_id)
                list_name.append(str(k)+"_"+str(barcode_id))
            count = count + 1
    
    with open(output_file + ".fasta","w") as f:
        for i in range(0,len(list_seq)):
           print_fasta(f,list_name[i],list_seq[i])

    with open(output_file+ ".para.txt","w") as f:
        print("abundance:",file = f)
        print(abundance,file = f)

        print("haploptyes:",file = f)
        for i in range(0,K):
           print_fasta(f,str(i),haplotypes[i])
        
        print("barcodes:",file = f)      
        for i in range(0,count):
           print_fasta(f,str(i),real_barcodes[i])


        print("true hap id:",file = f)
        print(list_true_hap,file = f)
        print("true UMI id:",file = f)
        print(list_true_UMI,file = f)

simulate_fasta()