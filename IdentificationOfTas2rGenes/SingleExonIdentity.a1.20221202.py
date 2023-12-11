#!/usr/bin/env python
# coding=utf-8
import sys
import os
import subprocess

fasta_file = sys.argv[2]
sequence = {}

#Type0=Find start codon
#Type1=no start codon
#Type2=no end codon
#Type3=no start and end codon
#Type4=sequence too short,maybe termination

#detect any start codon ATG site
def ATG_dect(seq) :
    Kmer = len(seq) - 3 + 1
    seq_list = []
    kmer_list = []
    for step in range(Kmer) :
        ATG_seq = {}
        codon_seq = seq[step:step+3]
        kmer_list.append(codon_seq)
        if codon_seq == "ATG" :
            new_seq = seq[step:len(value)]
            seq_start = step+1
            seq_end = len(seq)
            seq_pos = str(seq_start)+"-"+str(seq_end)
            ATG_seq[seq_pos] = new_seq
            seq_list.append(ATG_seq)
        else :
            pass
    if "ATG" in kmer_list :
        return ["Type0",seq_list]
    else :
        if "TAG" or "TAA" or "TGA" in kmer_list :
            return ["Type1",seq]
        else :
            return ["Type3",seq]

###detect partial sequence
def stop_dect(seq) :
    Kmer = len(seq) - 3 + 1
    dict_seq = {}
    for step in range(Kmer) :
        codon_seq = seq[step:step+3]
        if codon_seq in ["TAG","TAA","TGA"] :
            new_seq1 = seq[0:step+3]
            if len(new_seq1) % 3 == 0 :
                i = 0
                for x in range(int(len(new_seq1)/3)) :
                    new_seq2 = new_seq1[x*3:x*3+3]
                    if new_seq2 in ["TAG","TAA","TGA"] :
                        i += 1
                    else :
                        pass
                if i == 1 :
                    pos1 = "1-"+str(len(new_seq1))
                    dict_seq[pos1] = new_seq1
                else :
                    pass
            elif (len(new_seq1)-1) % 3 == 0 :
                i = 0
                for x in range(int(len(new_seq1)/3)) :
                    new_seq2 = new_seq1[x*3+1:x*3+4]
                    if new_seq2 in ["TAG","TAA","TGA"] :
                        i += 1
                    else :
                        pass
                if i == 1 :
                    pos1 = "2-"+str(len(new_seq1))
                    dict_seq[pos1] = new_seq1[1:]
                else :
                    pass
            elif (len(new_seq1)-2) % 3 == 0 :
                i = 0
                for x in range(int(len(new_seq1)/3)) :
                    new_seq2 = new_seq1[x*3+2:x*3+5]
                    if new_seq2 in ["TAG","TAA","TGA"] :
                        i += 1
                    else :
                        pass
                if i == 1 :
                    pos1 = "3-"+str(len(new_seq1))
                    dict_seq[pos1] = new_seq1[2:]
                else :
                    pass
    if not bool(dict_seq) :
        flag = "No any sequence"
        return flag
    else :
        new_dict = {}
        maximum_key = max(dict_seq, key=lambda k: sum(len(v) for v in dict_seq[k]))
        longest_seq = dict_seq[maximum_key]
        new_dict[maximum_key] = longest_seq
        return new_dict

def end_dect(seq) :
    step = int(len(seq)/3)
    codon_list = []
    full_seq = ""
    for i in range(step) :
        new_codon = seq[i*3:i*3+3]
        codon_list.append(new_codon)
        if new_codon in ["TAG","TAA","TGA"] :
            full_seq += seq[0:i*3+3]
            break
        else :
            pass
    if "TAG" or "TAA" or "TGA" in codon_list :
        if len(full_seq) >= 810 :
            return [len(full_seq),full_seq]
        else:
            return ["Type4",i*3+3,full_seq]
    else :
        return ["Type2",seq]

with open(fasta_file, "r") as fin :
    for line in fin :
        line = line.strip()
        if line.startswith(">") :
            gene = line.split()[0].strip(">")
            seqs = ""
        else:
            seqs += line
            sequence[gene] = len(seqs)

exc_sequence = {}
with open(sys.argv[1], "r") as fin1 :
    for line1 in fin1 :
        line1 = line1.strip().split("\t")
        chrom = line1[0]
        up = int(line1[1])
        down = int(line1[2])
        if up < down :
            if up - 500 > 0 and down + 500 <= int(sequence[chrom]) :
                seq_pos = chrom + ":" + str(up - 500) + "-" + str(down + 500)
                seq_id = chrom + "-" + str(up - 500) + "-" + str(down + 500) + "-pc\t501\t"+str(501+down-up)+"\t"+str(down-up+1+1000)+"\t"+str(up)+"\t"+str(down)
                #chrom_up-extend500+down-extend5000-strand \t new_start_in_new_exct_sequence \t new_end_in_new_exct_sequence \t full_exct_sequence_length \t raw_alignment_sequence_start \t raw_alignment_sequence_end
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()
            elif up - 500 <= 0 and down + 500 <= int(sequence[chrom]):
                seq_pos = chrom + ":1-" + str(down + 500)
                seq_id = chrom + "-1-" + str(down + 500) + "-pc\t"+str(up)+"\t"+str(down)+"\t"+str(down+500)+"\t"+str(up)+"\t"+str(down)
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()
            elif up - 500 > 0 and down + 500 > int(sequence[chrom]):
                seq_pos = chrom + ":" + str(up - 500) + "-" + str(int(sequence[chrom]))
                seq_id = chrom + "-" + str(up - 500) + "-" + str(int(sequence[chrom])) + "-pc\t501\t"+str(501+down-up)+"\t"+str(int(sequence[chrom])-up+501)+"\t"+str(up)+"\t"+str(down)
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()
            elif up - 500 <= 0 and down + 500 > int(sequence[chrom]):
                seq_pos = chrom + ":1-" + str(sequence[chrom])
                seq_id = chrom + "-1-" + str(sequence[chrom]) + "-pc\t"+str(up)+"\t"+str(down)+"\t"+str(sequence[chrom])+"\t"+str(up)+"\t"+str(down)
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()
        elif up > down :
            if down - 500 > 0 and up + 500 <= int(sequence[chrom]) :
                seq_pos = chrom + ":" + str(down - 500) + "-" + str(up + 500)
                seq_id = chrom + "-" + str(down - 500) + "-" + str(up + 500) + "-rc\t501\t"+str(501+up-down)+"\t"+str(up-down+1+1000)+"\t"+str(down)+"\t"+str(up)
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -i -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()
            elif down - 500 <= 0 and up + 500 <= int(sequence[chrom]) :
                seq_pos = chrom + ":1-" + str(up + 500)
                seq_id = chrom + "-1-" + str(up + 500) + "-rc\t501\t"+str(501+up-down)+"\t"+str(up+500)+"\t"+str(down)+"\t"+str(up)
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -i -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()
            elif down - 500 > 0 and up + 500 > int(sequence[chrom]) :
                seq_pos = chrom + ":" + str(down - 500) + "-" + str(int(sequence[chrom]))
                seq_id = chrom + "-" + str(down - 500) + "-" + str(int(sequence[chrom])) + "-rc\t"+str(int(sequence[chrom])-up+1)+"\t"+str(int(sequence[chrom])-down+1)+"\t"+str(int(sequence[chrom])-down+501)+"\t"+str(down)+"\t"+str(up)
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -i -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()
            elif down - 500 <= 0 and up + 500 > int(sequence[chrom]) :
                seq_pos = chrom + ":1-" + str(sequence[chrom])
                seq_id = chrom + "-1-" + str(sequence[chrom]) + "-rc\t"+str(int(sequence[chrom])-up+1)+"\t"+str(int(sequence[chrom])-down+1)+"\t"+str(sequence[chrom])+"\t"+str(down)+"\t"+str(up)
                exc_sequence[seq_id] = str(subprocess.check_output('samtools faidx %s %s -i -n 100000 | grep -v ">"'%(fasta_file,seq_pos),shell=True)).strip().upper()

fout1 = open("gene_identity.log", "w")
fout2 = open(sys.argv[1]+".raw_intact_gene.fasta", "w")
fout3 = open(sys.argv[1]+".no_intact.fasta", "w")
#fout4 = open(sys.argv[1]+".partial_sequence_no_stop_codon.fasta", "w")
#fout5 = open(sys.argv[1]+".partial_no_start_codon.fasta", "w")
#fout6 = open(sys.argv[1]+".partial_no_start_end_codon.fasta", "w")
#fout7 = open(sys.argv[1]+".aligned_sequence_too_short.fasta", "w")
fout8 = open(sys.argv[1]+"_500bp_extend", "w")

for key,value in exc_sequence.items() :
    fout8.write(">"+key+"\n"+value+"\n")
    chrom_1 = key.split("\t")[0].split("-")[0]
    start_1 = key.split("\t")[1]
    end_1 = key.split("\t")[2]
    length_1 = key.split("\t")[3]
    raw_start1 = key.split("\t")[4]
    raw_end1 = key.split("\t")[5]
    new_site1 = ">gi="+chrom_1+"\t"+start_1+"\t"+end_1+"\t"+length_1+"="+raw_start1+"\t"+raw_end1
    if len(value) < 810 :
        fout3.write(new_site1+"\n"+value+"\n")
    else :
        list = ATG_dect(value)
        strand = key.split("\t")[0].split("-")[3]
        chr1 = key.split("\t")[0].split("-")[0]
        start1 = key.split("\t")[0].split("-")[1]
        end1 = key.split("\t")[0].split("-")[2]
        if list[0] == "Type0" :
            no_stop_codon = {}
            full_sequence = {}
            too_short = {}
            for i in list[1] :
                for key1,value1 in i.items() :
                    list1 = end_dect(value1)
                    if len(list1) == 2 :
                        if list1[0] == "Type2":
                            no_stop_codon[key1] = list1[1]
                        else :
                            pos2 = key1.split("-")[0]+"-"+str(int(key1.split("-")[0])+len(list1[1])-1)
                            full_sequence[pos2] = list1[1]
                    else :
                        pos3 = key1.split("-")[0]+"-"+str(int(key1.split("-")[0])+len(list1[2])-1)
                        too_short[pos3] = list1[2]
            if bool(full_sequence) :
                new_dict1 = {}
                maximum_key = max(full_sequence, key=lambda k: sum(len(v) for v in full_sequence[k]))
                longest_seq = full_sequence[maximum_key]
#               new_dict1[maximum_key] = longest_seq
                up1=maximum_key.split("-")[0]
                down1 = maximum_key.split("-")[1]
                if strand == "pc" :
                    new_up =  int(start1)+int(up1)-1
                    new_end = int(start1)+int(down1)-1
                    fout1.write(key + " has full length sequence\n")
                    fout2.write(">"+chr1+":"+str(new_up)+"-"+str(new_end)+" "+key+" "+maximum_key+" postive\n"+longest_seq+"\n")
                if strand == "rc" :
                    new_up =  int(end1)-int(down1)+1
                    new_end = int(end1)-int(up1)+1
                    fout1.write(key + " has full length sequence\n")
                    fout2.write(">"+chr1+":"+str(new_up)+"-"+str(new_end)+" "+key+" "+maximum_key+" reverse\n"+longest_seq+"\n")
            elif bool(too_short) :
                fout1.write(key + " may be partial sequence too short\n")
                fout3.write(new_site1+"\n"+value+"\n")
            elif bool(no_stop_codon) :
                fout1.write(key + " may be partial sequence no stop codon\n")
                fout3.write(new_site1+"\n"+value+"\n")
        elif list[0] == "Type1" :
            list2 = stop_dect(value)
            if list2 == "No any sequence" :
                fout1.write(key + " No any sequence\n")
                fout3.write(new_site1+"\n"+value+"\n")
            else :
                fout1.write(key + " may be partial sequence no start codon\n")
                fout3.write(new_site1+"\n"+value+"\n")
        elif list[0] == "Type3" :
            fout1.write(key + " no start and stop codon\n")
            fout3.write(new_site1+"\n"+value+"\n")
