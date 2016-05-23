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
def pairwise( seq1, seq2, max_mismatch ):
	if len( seq1 ) != len( seq2 ):
		return False
	
	num_mismatch = 0
	
	for idx in range( len( seq1 ) ):
		if seq1[idx] != seq2[idx]:
			num_mismatch += 1
			
			if num_mismatch > max_mismatch:
				return False
			
	return num_mismatch <= max_mismatch

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
			
			for index in range( 0, len( dna_seq_to_translate ), 3 ):
				try:
					triplet = dna_seq_to_translate[index:index+3]
					curr_amino_acid = amino_acid_translations[ triplet ]
					peptide += curr_amino_acid
					
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

	
	# for key in local_peptide_dict:
		# if local_peptide_dict[key] == None:
			# print "Problem!"
			# return None
	
	return { 'sample':sample_name, 'peptide_dictionary':local_peptide_dict, 'dropped':local_dropped_count, 'codon_dictionary':local_codon_table }

def merge_proteins( protein_dictionary, proteins_to_combine ):
	###############################
	# combine DNA dictionaries    #
	# combine sample dictionaries #
	# add peptide counts          #
	# add pandas DataFrames       #
	###############################
	new_dna_dict = {}
	new_sample_dict = {}
	new_peptide_counts = 0
	
	tmp_df = protein_dictionary[ proteins_to_combine[0] ].get_df()
	new_df = pd.DataFrame( numpy.zeros( tmp_df.shape ), columns = tmp_df.columns.values.tolist() )
	
	for key in proteins_to_combine:
		curr_peptide = protein_dictionary[key]
		
		new_dna_dict = add_dictionaries( new_dna_dict, curr_peptide.dna_seqs )
		new_sample_dict = add_dictionaries( new_sample_dict, curr_peptide.sample_count )
		
		new_peptide_counts += curr_peptide.peptide_count
		
		new_df += curr_peptide.get_df()
	
	motif_peptide = ''.join( new_df.idxmax( axis=1 ) )
	
	return Peptide( motif_peptide, None, None, in_pep_count = new_peptide_counts, in_dna_dict = new_dna_dict, in_sample_dict = new_sample_dict, in_df = new_df )
	
def protein_collapse ( inDictionary, max_mismatch ):
	print "Collapse!"
	
	all_keys = inDictionary.keys()
	all_keys.sort()
	
	collapsed = []
	collapse_dict = {}
	
	for ref_idx in range( len(all_keys)-1 ):
		if ref_idx in collapsed:
			continue
			
		ref_seq = all_keys[ref_idx]
	
		for alt_idx in range( ref_idx + 1, len( all_keys ) ):
			if alt_idx in collapsed:
				continue
				
			alt_seq = all_keys[alt_idx]
				
			if pairwise( ref_seq, alt_seq, max_mismatch ):
				collapsed.append( alt_idx )
				
				if ref_seq in collapse_dict:
					collapse_dict[ ref_seq ].append( alt_seq )
				else:
					collapse_dict[ ref_seq ] = [ alt_seq ]
	
	print "Now collapsing keys"
	for key in collapse_dict:
		seq_group = collapse_dict[key]
		seq_group.append( key )
		new_peptide = merge_proteins( inDictionary, seq_group )
		
		for val in seq_group:
			del inDictionary[val]
			
		inDictionary[ new_peptide.peptide_seq ] = new_peptide
		
		# for val in collapse_dict[key]:
			# try:
				# inDictionary[key] += inDictionary[val]
				# del inDictionary[val]
			# except Exception as e:
				# print "Problem!"
				# print "Key: " + key
				# print "Val: " + val
				# print "inDictionary key type, val type"
				# print type( inDictionary[key] )
				# print type( inDictionary[val] )
				# print str( inDictionary[key] )
				# print str( inDictionary[val] )
				# print "Error: %s\n" % ( e )
				# sys.stdout.flush()
				# sys.exit()
	
	return inDictionary
	
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
	parser.add_argument( '--noCollapse', dest="runCollapse", default=True, required=False, action='store_false' )
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
	options += "Collapse motifs: %s\n" % ( args.runCollapse )

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
	
	#############################
	# check the damn keys again #
	#############################
	for key in globalPeptideDict:
		if globalPeptideDict[key] is None or type( globalPeptideDict[key] ) == "NoneType":
			print "Fracking problem is after merging"
			print "Key: %s" % ( key )
			
			for i in range( len( mp_results ) ):
				loop_sample = mp_results[i]['sample']
				loop_dict = mp_results[i]['peptide_dictionary']
				
				print "Sample: %s\nOutput: %s" % ( loop_sample, loop_dict.get( key, None ) )
			
			sys.exit(1)
	
	###################
	# collapse motifs #
	###################
	if args.runCollapse:
		globalPeptideDict = protein_collapse( globalPeptideDict, args.collapseMaxMismatch )
		
		for key in globalPeptideDict:
			globalPeptideDict[key].set_motif()
	
	###################################################
	# now write sensible output from the input tables #
	###################################################
	with open( "%s_counts.csv" % ( args.output ), 'w' ) as COUNTS_FH, open( "%s_seqs.csv" % ( args.output ), 'w' ) as DNAS_FH, open( "%s_codon_table.csv" % ( args.output ), 'w' ) as CODON_FH:
		
		# write out the peptides and DNAs per peptide
		COUNTS_FH.write( "Peptide,%s\n" % ( ','.join( sample_list ) ) )
		
		for currPeptide in globalPeptideDict.keys():
			if globalPeptideDict[ currPeptide ].peptide_count >= args.peptideMinCount:
				# dna
				DNAS_FH.write( "%s,%s\n" % ( currPeptide, '|'.join( globalPeptideDict[ currPeptide ].dna_seqs.keys() ) ) )
				
				# write samples
				COUNTS_FH.write( "%s" % ( currPeptide ) )
				for currSample in sample_list:
					COUNTS_FH.write( ",%s" % ( globalPeptideDict[ currPeptide ].get_sample_count( currSample ) ) )
				COUNTS_FH.write( "\n" )
			else:
				for currSample in sample_list:
					if currSample in droppedCount.keys():
						droppedCount[currSample] += globalPeptideDict[ currPeptide ].get_sample_count( currSample )
					else:
						droppedCount[currSample] = globalPeptideDict[ currPeptide ].get_sample_count( currSample )

		COUNTS_FH.write( "__dropped" )
		
		for currSample in sample_list:
			if currSample in droppedCount:
				dropped_num = droppedCount[currSample]
			else:
				dropped_num = 0

			COUNTS_FH.write( ",%s" % ( dropped_num ) )
		COUNTS_FH.write( '\n' )
		
		# write codon usage table
		all_codons_obs = []
		
		for currSample in global_codon_table:
			for codon in global_codon_table[currSample]:
				all_codons_obs.append( codon )
		all_codons_obs = list( set( all_codons_obs ) )
		all_codons_obs.sort()
		
		CODON_FH.write( "Codon,%s\n" % ( ','.join( sample_list ) ) )
		
		for codon in all_codons_obs:
			CODON_FH.write( "%s" % ( codon ) )
			
			for sample in sample_list:
				CODON_FH.write( ",%s" % ( global_codon_table[sample].get( codon, 0 ) ) )
			CODON_FH.write( "\n" )
