My postdoc research involves finding patients with retinal degenerative mutations that are amenable for base editing correction. Later, we would like to develop "disease-in-a-dish" model that recapitulates their degeneration features using iPSCs, and to try to reverse the disease phenotype by correction of the mutation using base editing techniques. Potentially, this kind of treatment will join the arsenal of genetic treatments like gene editing (Crispr-KO,i,a) and gene augmentation, that are either have already FDA approval (luxturna - gene augmentation using AAV for genetic retinal degenerative disease) or are on clinical trial phase.

The first step in my research is to find the right candidates for the treatment. We work on a data base of patients that were found to have been diagnosed with retinal degeneration clinically, and had Whole Exome Sequencing that found the presumably causative pathologic point mutation.

Base editing mechanism has several limitations:

There are (at this point) mainly Adenine and Cytidine base editors. Adenine changes A/T to G/C. Cytidine C/G to T/A. That means that one can treat G/C to A/T or T/A to C/G mutations. {-/- the second (-) means the complement base}

Each base editor has different activity window. That means that we can have unwanted editing of bases other than the point mutation - "a bystander"

Base editors are part of complex that incorporates cas9, which guides together with the guide RNA the complex to the target mutation. cas9 is responsible to bind the PAM site downstream to the mutation. Each cas9 can bind specific PAM sequences.

The way I used to check whether a patient fits for base editing is with a tool from a website (rgenome.net/be-designer/). It asks you to input a sequence (raw or fasta) and choose the type of cas9 and type of base editor. Unfortunately, there is no way to choose ALL the cas9 types in one step and there are 13 kinds of it, so for every sequence I had to run the website algorithm 13 times..And I had about 30 patient sequences.. 
The output of this program gives you the gRNA if it finds one for the respected input, activity window, and bystanders.. The gRNA is usually 20 base long sequence DNA sequence upstream to the pam site. This sequence guides the BE complex to the complement strand of the mutation strand that we want to edit. The mutation should be in the activity window - for Adenine base editor is at positions 14-16 of the gRNA counting from the 3' end (= the pam site 5').

For my final project I decided to design an algorithem that will receive a list of sequences where the pathogenic point mutation is the only one marked in capital letter. Another input is the type of base editor. The output will be a list of dictionaries where {seq(n): [cas9, gRNA, activity window, number of bystanders]} for each sequence in the input list in a "results.txt" text file, and in the return of the function for future coding manipulations.
