#!/usr/bin/env python

######################
# E. Roberson        #
# Created 2015-11-20 #
######################

###########
# Imports #
###########
import argparse
import fastqp
import multiprocessing as mp
import numpy
import pandas as pd
import re
import sys
import time
from collections import Counter
import random
import operator

####################
# Version and name #
####################
SCRIPT_PATH = sys.argv[0]
SCRIPT_NAME = SCRIPT_PATH.split('/')[-1].split('\\')[-1]
VERSION = '1.1.0'

########
# fxns #
########
def add_dictionaries( d1, d2 ):
        for key in d2:
                if key in d1:
                        d1[key] = d1[key] + d2[key]
                else:
                        d1[key] = d2[key]
        return d1

#####################
# class Definitions #
#####################
class Peptide:
        def __init__( self, amino_acids, dna, sample, amino_acid_list=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y'], in_pep_count=None, in_dna_dict=None, in_sample_dict=None, in_df=None ):
                self.peptide_seq = amino_acids

                if in_pep_count is None:
                        self.peptide_count = 1
                else:
                        self.peptide_count = in_pep_count

                if in_dna_dict is None:
                        self.dna_seqs = {dna:1}
                else:
                        self.dna_seqs = in_dna_dict

                if in_sample_dict is None:
                        self.sample_count = {sample:1}
                else:
                        self.sample_count = in_sample_dict

                if in_df is None:
                        self.df = pd.DataFrame( pd.np.empty(( len( self.peptide_seq ), len( amino_acid_list ) ))*0, columns = amino_acid_list )
                        self.df_set = False
                else:
                        self.df = in_df
                        self.df_set = True

        def update( self, dna, sample ):
                self.peptide_count += 1
                self.dna_seqs[dna] = self.dna_seqs.get( dna, 0 ) + 1
                self.sample_count[sample] = self.sample_count.get( sample, 0 ) + 1

        def get_sample_count( self, sample ):
                return self.sample_count.get( sample, 0 )

        def set_up_df( self ):
                if self.df_set == False:
                        for idx in range( len( self.peptide_seq ) ):
                                self.df.ix[idx, self.peptide_seq[idx]] = self.peptide_count
                        self.df_set = True
	def get_df( self ):
                if self.df_set == False:
                        self.set_up_df()
                return self.df

        def set_motif( self ):
                self.peptide_seq = ''.join( self.df.idxmax( axis=1 ) )

        def __add__( self, other ):
                        tmp_peptide = self.peptide_seq
                        tmp_count = self.peptide_count + other.peptide_count
                        tmp_seq_dict = add_dictionaries( self.dna_seqs, other.dna_seqs )
                        tmp_samp_dict = add_dictionaries( self.sample_count, other.sample_count )
                        tmp_pd = self.df + other.df

                        return Peptide( tmp_peptide, None, None, in_pep_count=tmp_count, in_dna_dict=tmp_seq_dict, in_sample_dict=tmp_samp_dict, in_df=tmp_pd )

        def __str__( self ):
                s = ''
                s += "Peptide: %s\n" % ( self.peptide_seq )
                s += "Peptide count: %s\n" % ( self.peptide_count )
                s += "\nDNAs\n"

                for key in self.dna_seqs:
                        s += "%s %s\n" % ( key, self.dna_seqs[key] )
                return s
########
# fxns #
########
def pairwise( seq1, seq2 ):
        if len( seq1 ) != len( seq2 ):
                return False

        num_mismatch = 0
        num_mismatch = sum(map(lambda x, y:x!=y, seq1,seq2))

        return num_mismatch

def process_fastq_file( sample_name, file_name, peptide_dna_length, five_prime_seq, amino_acid_translations ):
        ################
        # accumulators #
        ################
        local_peptide_dict = {}
        local_dropped_count = {sample_name:0}
        local_codon_table = {}

        ###############
        # regex setup #
        ###############
        five_prime_re = re.compile( five_prime_seq )
		
	#########################
        # read whole fastq file #
        #########################
        with fastqp.FastqReader( open( file_name, 'r' ) ) as FQOBJ:
                for record in FQOBJ:
                        match_it = five_prime_re.search( record.seq )

                        if not match_it:
                                continue

                        dna_seq_to_translate = record.seq[ ( match_it.end() ):( match_it.end() + peptide_dna_length ) ]

                        peptide = ''

                        kill_fragment = False

                        for i in range( 0, len( dna_seq_to_translate ), 3 ):
                                try:
                                        triplet = dna_seq_to_translate[i:i+3]
                                        peptide += amino_acid_translations[ triplet ]
                                        local_codon_table[ triplet ] = local_codon_table.get( triplet, 0 ) + 1
                                except:
                                        kill_fragment = True
                                        break

                        if kill_fragment == True:
                                local_dropped_count[sample_name] += 1
                        elif peptide in local_peptide_dict:
                                local_peptide_dict[peptide].update( dna_seq_to_translate, sample_name )
                        else:
                                local_peptide_dict[peptide] = Peptide( peptide, dna_seq_to_translate, sample_name )

        return { 'sample':sample_name, 'peptide_dictionary':local_peptide_dict, 'dropped':local_dropped_count, 'codon_dictionary':local_codon_table }

###############################
# Fxns for K-means Clustering #
###############################
def kmeans(sequences, k, max_iterations):
    centroid_new = random.sample(sequences, k)
    iterations = 0;
    while (iterations != max_iterations):
        centroid_old = centroid_new
        cluster_all = cluster_sequences(sequences, centroid_new)
        centroid_new = recenter(centroid_old, cluster_all)
        iterations += 1
    return cluster_all

def cluster_sequences(sequences, centroid_new):
    clusters = {}
    all_keys = sequences.keys()
    for alt_idx in range(len(all_keys)):
        alt_seq = all_keys[alt_idx]
        best_key = min([(centroid_new[i[0]], pairwise(centroid_new[i[0]], alt_seq)) for i in enumerate(centroid_new)], key=lambda t:t[1])[0]
        try:    
            clusters[best_key].append(alt_seq)
        except KeyError:
            clusters[best_key] = [alt_seq]
    return clusters

def recenter(centroid_old, clusters):
    centroid_new = []
    keys = sorted(clusters.keys())
    for k in keys:
        centroid_new.append(consensus(clusters[k]))
    return centroid_new

def consensus(cluster):
    consensus_seq = ""
    for j in range(len(cluster[0])):
        dummy = []
        for i in range(len(cluster)):
	    try:
                dummy.append(cluster[i][j])
            except:
		dummy.append(" ")
        top_peptide = Counter(dummy).most_common(1)
        consensus_seq += str(top_peptide[0][0])
    return consensus_seq

############
# run main #
############
if __name__ == "__main__":
        ##############################
        # constant translation table #
        ##############################
        translationDict = {'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'TGC':'C', 'TGT':'C', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'TTC':'F', 'TTT':'F', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'CAC':'H', 'CAT':'H', 'ATA':'I', 'ATC':'I', 'ATT':'I', 'AAA':'K', 'AAG':'K', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'TTA':'L', 'TTG':'L', 'ATG':'M', 'AAC':'N', 'AAT':'N', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAA':'Q', 'CAG':'Q', 'AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'TGG':'W', 'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X', 'TGA':'X'}

        ############
        # argparse #
        ############
        parser = argparse.ArgumentParser( prog=SCRIPT_NAME, epilog="%s v%s" % (SCRIPT_NAME, VERSION) )

        parser.add_argument( '--input', help="sample,FASTQ_path. Can be specified multiple times", action='append', required=False )
        parser.add_argument( '--input_file', help="path to CSV file in sample,FASTQ_path format", required=False )
        parser.add_argument( '--output', help="Output basename", required=False, default="phage_enrichments" )
        parser.add_argument( '--peptideMinCount', help="Minimum times a peptide must be seen overall to print", default=1, type=int )
        parser.add_argument( '--fivePrimeSeq', help="Sequence upstream of insert", default='TCCTTTAGTGGTACCTTTCTATTCTCACTCT' )
        parser.add_argument( '--peptideLength', help="Size of insert peptide", default=10, type=int )
        parser.add_argument( '--cores', help="Number of computer cores to use", default=1, type=int )
        parser.add_argument( '--collapseMaxMismatch', default=1, type=int )
        parser.add_argument( '--maxIterations', default=100, type=int )
        parser.add_argument( '--clusterNumber', default=100, type=int) 
        parser.add_argument( "--quiet", default=True, action='store_false', dest='verbose' )
        args = parser.parse_args()

        ################
        # fix up input #
        ################
        if args.input is None:
                args.input = []

        if args.input_file is not None:
                with open( args.input_file, 'r' ) as IN:
                        for line in IN:
                                line = line.rstrip()

                                if len( line ) > 0:
                                        args.input.append( line )

        if len( args.input ) == 0:
                sys.exit( "No input specified! Use --input and/or --input_file" )

        args.input = list( set( args.input ) ) # uniquifies the set
        args.input.sort()

        #####################
        # write out options #
        #####################
        options = "%s v%s\n\n" % ( SCRIPT_NAME, VERSION )
        options += "Options\n"
        options += "=======\n"
        options += "Input FASTQ: %s\n" % ( '; '.join( args.input ) )
        options += "Output File: %s\n" % ( args.output )
        options += "Cores: %s\n" % ( args.cores )
        options += "Peptide Minimum Count: %s\n" % ( args.peptideMinCount )
        options += "5' sequence to trim: %s\n" % ( args.fivePrimeSeq )
        options += "Peptide length to translate: %s\n" % ( args.peptideLength )
        options += "Max mismatch to collapse: %s\n" % ( args.collapseMaxMismatch )
        options += "Cluster Number: %s\n" % (args.clusterNumber )
        options += "Iterations: %s\n" % ( args.maxIterations )

        if args.verbose == True:
                print options

        peptideDnaLength = args.peptideLength * 3

        #######################
        # Process FASTQ files #
        #######################
        sample_list = []
        mp_inputs = []

        for input in args.input:
                sample, file = input.rstrip().split( ',' )
                mp_inputs.append( (sample, file) )

        ##########################################
        # run fastq processing on multiple cores #
        ##########################################
        proc_workers = mp.Pool( processes=args.cores )
        run_pool = [ proc_workers.apply_async( process_fastq_file, args=(x[0], x[1], peptideDnaLength, args.fivePrimeSeq, translationDict ) ) for x in mp_inputs ]
        mp_results = [ x.get() for x in run_pool ]

	#####################################################
        # merge data from sub procs into global report objs #
        #####################################################

        globalPeptideDict = {}
        droppedCount = {}
        global_codon_table = {}

        for i in range( len(mp_results) ):
                print "Merging " + mp_results[i]['sample']

                sample_list.append( mp_results[i]['sample'] )

                globalPeptideDict = add_dictionaries( globalPeptideDict, mp_results[i]['peptide_dictionary'] )
                droppedCount = add_dictionaries( droppedCount, mp_results[i]['dropped'] )
                global_codon_table = add_dictionaries( global_codon_table, { mp_results[i]['sample']:mp_results[i]['codon_dictionary'] } )

        sample_list = list( set( sample_list ) ) # uniqify sample_list
        
	########################
        # check the keys again #
        ########################
        for key in globalPeptideDict:
                if globalPeptideDict[key] is None or type( globalPeptideDict[key] ) == "NoneType":
                        print "Fracking problem is after merging"
                        print "Key: %s" % ( key )

                        for i in range( len( mp_results ) ):
                                loop_sample = mp_results[i]['sample']
                                loop_dict = mp_results[i]['peptide_dictionary']

                                print "Sample: %s\nOutput: %s" % ( loop_sample, loop_dict.get( key, None ) )

                        sys.exit(1)
                        
        ##################
        # Cluster Motifs #
        ##################
        
        kmeansCluster = kmeans(globalPeptideDict, args.clusterNumber,  args.maxIterations)
                        
	###################################################
        # now write sensible output from the input tables #
        ###################################################
        with open("unique_seq_counts.txt", "w") as f:
        	for key in kmeansCluster:
                	f.write(key + ":" + str( len(kmeansCluster[key])) +"\n")
        f.close()

	output = ['Motif']
	for currSample in sample_list:
		output.append(currSample)

	for key in kmeansCluster:
		dummy = []
		output.append(key);
		for currSample in sample_list:
			cSum = 0;
			for seq in kmeansCluster[key]:
				cSum = cSum + globalPeptideDict[seq].get_sample_count(currSample)
			output.append(cSum)
	with open("sample_counts_table.txt", "w") as f:
		a=0
		for item in output:
			if a < len(sample_list):
				f.write(str(item) + "\t")
				a = a+1
			else:
				f.write(str(item) +  "\n")
				a = 0
	f.close()

	with open("motif_dnaseq_abundance.txt", "w") as f:
		for motif in kmeansCluster.keys():
			k=0;
			for seq in kmeansCluster[motif]:
				if k==0:
					f.write(str(motif))
				else: 
					f.write("\t")
				f.write("\t" + str(seq))
				for sample in sample_list:
					f.write("\t"+str(globalPeptideDict[seq].get_sample_count(sample)))
				f.write("\n")
				k = k+1;
	f.close()
