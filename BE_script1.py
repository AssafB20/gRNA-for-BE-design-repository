#!/usr/bin/env python
# coding: utf-8

# My postdoc research involves finding patients with retinal degenerative mutations that are amebale for base editing correction.
# Later, we would like to develop "disease-in-a-dish" model that recapitulates their degeneration features using iPSCs, and to try to reverse the disease phenotype by correction of the mutation using base editing techniques. Potentially, this kind of treatment will join the arsenal of genetic treatments like gene editing (Crispr-KO,i,a) and gene augmentation, that are either have already FDA approval (luxturna - gene augmentation using AAV for genetic retinal degenerative disease) or are on clinical trial phase (EDITAS's Crispr-KO of dominant genetic retinal degenerative disease).
# 
# The first step in my research is to find the right candidates for the treatment. We work on a data base of patients that were found to have been diagnosed with retinal degeneration clinically, and had Whole Exome Sequencing that found the presumably causative pathologic point mutation.
# 
# Base editing mechanism has several limitations:
# 1. There are (at this point) mainly Adenine and Cytidine base editors. Adenine changes A/T to G/C. Cytidine C/G to T/A. That means that one can treat G/C to A/T or T/A to C/G mutations. {-/- the second (-) means the complement base}
# 
# 2. Each base editor has different activity window. That means that we can have unwanted editing of bases other than the point mutation - "a bystander"
# 
# 3. Base editors are part of complex that incorporates cas9, which guides together with the guide RNA the complex to the target mutation. cas9 is responsible to bind the PAM site downstream to the mutation. Each cas9 can bind specific PAM sequences.
# 
# The way I used to check whether a patient fits for base editing is with a tool from a website (rgenome.net/be-designer/).
# It asks you to input a sequence (raw or fasta) and choose the type of cas9 and type of base editor.
# Unfortunately, there is no way to choose ALL the cas9 types in one step and there are 13 kinds of it, so for every sequence I had to run the website algorithm 13 times..And I had about 30 patient sequences.. 
# The output of this program gives you the gRNA if it finds one for the respected input, activity window, and bystanders..
# The gRNA is usually 20 base long sequence DNA sequence upstream to the pam site. This sequence guides the BE complex to the complement strand of the mutation strand that we want to edit. The mutation should be in the activity window - for Adenine base editor is at positions 14-16 of the gRNA counting from the 3' end (= the pam site 5').
# 
# For my final project I decided to design an algorithem that will receive a list of sequences where the pathogenic point mutation is the only one marked in capital letter. Another input is the type of base editor. The output will be a list of dictionaries where {seq(n): [cas9, gRNA, activity window, number of bystanders]} for each sequence in the input list in a "results.txt" text file, and in the return of the function for future coding manipulations.
# 
# This way we can potentially save a lot of time in the next base editing target searches. The only time consuming in this process is tagging the pathogenic mutation with capital letter, but this is already how the person from the sequencing core plugging the reading to the data base.
# 
# I did all the coding from scratch because I couldn't find any website or any paper that let you do exactly what I was trying to do, and also because there is no available code from the website I was using at rgenome.net/be-designer/ .
# 
# I chose to use python because we were learning in class how to manipulate strings mainly in this language. I have found regex tools to be very efficient in manipulating long strings. Also, I ended up with few loops in my script, so this is another reason why choosing python for this task.
# 
# The most "labor intensive" part of the code was to write donwn the dictionary "pam", which I did manually..I only realized mid-way that I could have write a function to do it..but it was too late.
# 
# I compared my results file to the work I have already done before through the website, and I have found them to be identical.
# This new method will not only save time, but also solves the human error factor from doing multiple repetitions of copy pasting from the website to a file.
# 
# I believe that after I add a feature for finding possible "miss-matches" throughout the whole genome (which I need to learn to do maybe with BioPython?), my script can be valuable for many people in this field who design their base editing experiments.

# In[1]:


pam = {'spCas9':['agg','cgg','ggg','tgg'], 
       'spCas9-VQR':['agaa','agac','agag','agat','cgaa','cgac','cgag','cgat','ggaa','ggac','ggag','ggat','tgaa','tgac','tgag','tgat'],
       'spCas9-EQR':['agag','cgag','ggag','tgag'],'spCas9-VRER':['agcg','cgcg','ggcg','tgcg'],
       'saCas9':['aagaat','aagagt','aaggat','aagggt','acgaat','acgagt','acggat','acgggt','aggaat','aggagt','agggat','aggggt','atgaat','atgagt','atggat','atgggt',
                 'cagaat','cagagt','caggat','cagggt','ccgaat','ccgagt','ccggat','ccgggt','cggaat','cggagt','cgggat','cggggt','ctgaat','ctgagt','ctggat','ctgggt',
                 'gagaat','gagagt','gaggat','gagggt','gcgaat','gcgagt','gcggat','gcgggt','gggaat','gggagt','ggggat','gggggt','gtgaat','gtgagt','gtggat','gtgggt',
                 'tagaat','tagagt','taggat','tagggt','tcgaat','tcgagt','tcggat','tcgggt','tggaat','tggagt','tgggat','tggggt','ttgaat','ttgagt','ttggat','ttgggt'],
       'saCas9-KKH':
                ['aaaaat','aaaagt','aaagat','aaaggt','acaaat','acaagt','acagat','acaggt','agaaat','agaagt','agagat','agaggt','ataaat','ataagt','atagat','ataggt',
                 'caaaat','caaagt','caagat','caaggt','ccaaat','ccaagt','ccagat','ccaggt','cgaaat','cgaagt','cgagat','cgaggt','ctaaat','ctaagt','ctagat','ctaggt',
                 'gaaaat','gaaagt','gaagat','gaaggt','gcaaat','gcaagt','gcagat','gcaggt','ggaaat','ggaagt','ggagat','ggaggt','gtaaat','gtaagt','gtagat','gtaggt',
                 'taaaat','taaagt','taagat','taaggt','tcaaat','tcaagt','tcagat','tcaggt','tgaaat','tgaagt','tgagat','tgaggt','ttaaat','ttaagt','ttagat','ttaggt',
                 
                 'aacaat','aacagt','aacgat','aacggt','accaat','accagt','accgat','accggt','agcaat','agcagt','agcgat','agcggt','atcaat','atcagt','atcgat','atcggt',
                 'cacaat','cacagt','cacgat','cacggt','cccaat','cccagt','cccgat','cccggt','cgcaat','cgcagt','cgcgat','cgcggt','ctcaat','ctcagt','ctcgat','ctcggt',
                 'gacaat','gacagt','gacgat','gacggt','gccaat','gccagt','gccgat','gccggt','ggcaat','ggcagt','ggcgat','ggcggt','gtcaat','gtcagt','gtcgat','gtcggt',
                 'tacaat','tacagt','tacgat','tacggt','tccaat','tccagt','tccgat','tccggt','tgcaat','tgcagt','tgcgat','tgcggt','ttcaat','ttcagt','ttcgat','ttcggt',
           
                 'aagaat','aagagt','aaggat','aagggt','acgaat','acgagt','acggat','acgggt','aggaat','aggagt','agggat','aggggt','atgaat','atgagt','atggat','atgggt',
                 'cagaat','cagagt','caggat','cagggt','ccgaat','ccgagt','ccggat','ccgggt','cggaat','cggagt','cgggat','cggggt','ctgaat','ctgagt','ctggat','ctgggt',
                 'gagaat','gagagt','gaggat','gagggt','gcgaat','gcgagt','gcggat','gcgggt','gggaat','gggagt','ggggat','gggggt','gtgaat','gtgagt','gtggat','gtgggt',
                 'tagaat','tagagt','taggat','tagggt','tcgaat','tcgagt','tcggat','tcgggt','tggaat','tggagt','tgggat','tggggt','ttgaat','ttgagt','ttggat','ttgggt',

                 'aataat','aatagt','aatgat','aatggt','actaat','actagt','actgat','actggt','agtaat','agtagt','agtgat','agtggt','attaat','attagt','attgat','attggt',
                 'cataat','catagt','catgat','catggt','cctaat','cctagt','cctgat','cctggt','cgtaat','cgtagt','cgtgat','cgtggt','cttaat','cttagt','cttgat','cttggt',
                 'gataat','gatagt','gatgat','gatggt','gctaat','gctagt','gctgat','gctggt','ggtaat','ggtagt','ggtgat','ggtggt','gttaat','gttagt','gttgat','gttggt',
                 'tataat','tatagt','tatgat','tatggt','tctaat','tctagt','tctgat','tctggt','tgtaat','tgtagt','tgtgat','tgtggt','tttaat','tttagt','tttgat','tttggt'
                ],
       'xCas9':['agt','cgt','ggt','tgt'],
       'Cas9-NG (NGK>NGM)':['ag','cg','gg','tg']}


# 
# SpCas9 from Streptococcus pyogenes: 5'-NGG-3'
# 
# SpCas9-VQR from Streptococcus pyogenes: 5'-NGAN-3'
# 
# SpCas9-EQR from Streptococcus pyogenes: 5'-NGAG-3'
# 
# SpCas9-VRER from Streptococcus pyogenes: 5'-NGCG-3'
# 
# SaCas9 from Staphylococcus aureus: 5'-NNGRRT-'3
# 
# SaCas9-KKH from Staphylococcus aureus: 5'-NNNRRT-'3
# 
# xCas9 3.7 (TLIKDIV SpCas9) from Streptococcus pyogenes: 5'-NGT-3'
# 
# N - any base
# 
# R - A or G
# 
# K - G or T
# 
# M - A or C

# In[2]:


def comp_strand(strand):  #this function gets a sequence and returns its complement strand with keeping the mutation in capital.
 old_chars = "TGacgt"
 replace_chars = "ACtgca"
 tab = str.maketrans(old_chars,replace_chars)
 return strand.translate(tab)[::-1]   # returns the complement strand 5' -> 3'


# In[3]:


def targetfinder_file(filename, be): #filename - a text file with a list of sequences with mutation in capital letter, be - base-editor
   
    import re
    y={}
    x=0
    mut_loc=0
    #f= open('results.txt','w')
    temp = open(filename,'r').read().split('\n\n')       #the seqs in the file are separated by blank lines..so split the seqs
    for seq in temp:
        seq = seq.replace(" ","")             # for each seq, ignore the spaces within
        if be == 'a':                       # input 'a' was for adenine base editor
            mut_loc = seq.find('A')            # find the position of 'A' in the seq - which is the ONLY capital letter in it
            win = [14,15,16]                 # the activity window of adenine base editor
            if mut_loc == -1:               # if it didn't find 'A' in the seq, I guess the user gave us a strand with a 'T' mutation
                seq = comp_strand(seq)      # call a function that will give the complement strand, but keep the mutation in capital letter 
                mut_loc = seq.find('A')     # and now search for the 'A' position
                
        if be == 'c':                      # same goes for c base editor
            mut_loc = seq.find('C')
            win = [13,14,15,16,17,18,19]
            if mut_loc == -1:
                seq = comp_strand(seq)
                mut_loc = seq.find('C')
        
        for key, value in pam.items():                # pam is a dic with key-cas9 and values-their pam sites
                for p in value:                       
                    for x in re.finditer(p, seq[mut_loc:]):           # for every pam sites that is located AFTER the mutation position...because base editors' pams are downstream to the mutation
                        x = x.start() + mut_loc                       # x is the position of the first base of the pam site
                        if mut_loc in range(x-win[-1],x-win[0]+1):          # if the mutation located in the activity window
                            bystand = seq.count(be,x-win[-1],x-win[0]+1)    # search for another 'a' (or 'c') in this window
                            editwin = seq[x-win[-1]:x-win[0]+1]             # the window sequence itself
                            y.setdefault(seq,[]).append([key,'5-'+seq[x-20:x+len(p)]+'-3',editwin,bystand]) # append to a dic where key is the seq, values are the type of cas9, gRNA, editing window, how many if any bystanders
                        
    #for k,v in y.items():
     #   f.write(str(k) + ' >>> ' + str(v) + '\n\n')          #write the results (list of dic items) to a text file
                          
    open(filename,'r').close()
    #f.close()
    return y           #besides writing it to a text file, I thought it was important to return the list for future coding


# In[4]:


def targetfinder_string(seq_list, be): #filename - a text file with a list of sequences with mutation in capital letter, be - base-editor
    #be = input('for adenine base editing send a, for cytidine base editing send c: ')
    import re
    y={}
    x=0
    mut_loc=0
    #f= open('results.txt','w')
    seq_list = seq_list.split('\n\n')       #the seqs in the file are separated by blank lines..so split the seqs
    for seq in seq_list:
        seq = seq.replace(" ","")             # for each seq, ignore the spaces within
        if be == 'a':                       # input 'a' was for adenine base editor
            mut_loc = seq.find('A')            # find the position of 'A' in the seq - which is the ONLY capital letter in it
            win = [14,15,16]                 # the activity window of adenine base editor
            if mut_loc == -1:               # if it didn't find 'A' in the seq, I guess the user gave us a strand with a 'T' mutation
                seq = comp_strand(seq)      # call a function that will give the complement strand, but keep the mutation in capital letter 
                mut_loc = seq.find('A')     # and now search for the 'A' position
                
        if be == 'c':                      # same goes for c base editor
            mut_loc = seq.find('C')
            win = [14,15,16]
            if mut_loc == -1:
                seq = comp_strand(seq)
                mut_loc = seq.find('C')
        
        for key, value in pam.items():                # pam is a dic with key-cas9 and values-their pam sites
                for p in value:                       
                    for x in re.finditer(p, seq[mut_loc:]):           # for every pam sites that is located AFTER the mutation position...because base editors' pams are downstream to the mutation
                        x = x.start() + mut_loc                       # x is the position of the first base of the pam site
                        if mut_loc in range(x-win[-1],x-win[0]+1):          # if the mutation located in the activity window
                            bystand = seq.count(be,x-win[-1],x-win[0]+1)    # search for another 'a' (or 'c') in this window
                            editwin = seq[x-win[-1]:x-win[0]+1]             # the window sequence itself
                            y.setdefault(seq,[]).append([key,'5-'+seq[x-20:x+len(p)]+'-3',editwin,bystand]) # append to a dic where key is the seq, values are the type of cas9, gRNA, editing window, how many if any bystanders
                        
    #for k,v in y.items():
     #   f.write(str(k) + ' >>> ' + str(v) + '\n\n')          #write the results (list of dic items) to a text file
                          
    #open(filename,'r').close()
    #f.close()
    return y           #besides writing it to a text file, I thought it was important to return the list for future coding

