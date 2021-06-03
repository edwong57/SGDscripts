#!/usr/bin/python
import sys, getopt, os
import re
from operator import itemgetter

## Parse Pelechano datafiles ###

def main(argv):

	sorted_keys = list()
	data_to_filter = dict()
	two_fold_data = dict()
	genome_data = dict()
	transcript_data = dict()
	feat_type = 'gene' ## specific feature type; use 'all' if no specific feature type
	filter_val = 3  ## number to filter by (count or score)
	filter_type = 'count'  # can filter by count('count') or score cutoff ('cutoff')
	overlap = '1' # amount of overlap a transcript has on an ORF
	
	try:
		opts, args = getopt.getopt(argv,"hi:o:l:f:v:t:s",["ifile=","ofile=", "lfile=","format=","value=","type=","seqfeat="])
	except getopt.GetoptError:
		print 'parse_files_for_transcripts.py -i <inputfile> -o <outputfile> -l <locusfile> -f <output_file_format> -v <filter value> -t <filter_type: count OR cutoff> -s <sequence_feature: gene, CDS, etc>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'parse_files_for_transcripts.py -i <inputfile> -o <outputfile> -l <locusfile> -f <output_file_format> -v <filter value> -t <filter_type: count OR cutoff> -s <sequence_feature: gene, CDS, etc>'
			sys.exit()
		elif opt in ("-i", "--ifile"): # gff3 format
			firstfile = arg
		elif opt in ("-o", "--ofile"):  # outfile name
			outfile = arg
		elif opt in ("-l","--lfile"):  # s_cer gff3 file
			locusfile = arg
		elif opt in ("-f","--format"):  #gff3 or tsv/wig file ## DOESN'T WORK
			file_format = arg
		elif opt in ("-v","--value"): ## filter value 
			filter_val = arg
		elif opt in ("-t","--type"): ## type of filter -- count or cutoff # Not tested
			filter_type = arg
		elif opt in ("-s","--seqfeat"): ## feature type # Not tested
			feat_type = arg
			
	# defining out files:
	unfiltered_file = "unfiltered_" + outfile
	non_match_file = "unmatched_" + outfile
	
		
	# make hashes
	
	(file_data) = _open_make_hash(firstfile)  # opens file to parse and makes a hash
	(genome_data) = _parse_sac_gff(locusfile) # opens and makes hash of gene annotations
	
#	print "file keys: " + ",".join(file_data.keys())
#	print "s_cer gff keys: "+ ",".join(genome_data.keys())
	
	## iterate over all the 'gene' feature types and find the transcripts that cover the entire thing
	#, non_matches)
	
	(unfiltered_matches, filtered_transcripts, unmatched_transcripts) = _get_transcripts(file_data, genome_data, filter_type, filter_val, feat_type, overlap)
	
	## right now, just make a gff file
	_print_gff(filtered_transcripts, outfile)
	_print_gff(unfiltered_matches, unfiltered_file)
	_print_gff(unmatched_transcripts, non_match_file)
		
	## if -f flag is used, then just reformat data ###
# 	if 'format' in vars():
# #		reformatted_data = _reformat_data(data_to_filter, mapfile)
# 		
# 				### creating reformatted file
# 	
# 		final_file = open(outfile, 'w')
# 	
# 		## print header row
# 		# final_file.write ("regulator feature name\tregulator gene name\ttarget feature name\ttarget gene name\tvalue\tstrain\n")
# 		

def _print_gff(dict_to_print, outfile):
	print "making GFF file"
	#print headers #
	newfile = open(outfile, 'w')
	newfile.write("## GFF file for transcripts\n")
	
	chr_order = ['chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI']
	
	for chr in chr_order:
		transcript_data = dict_to_print[chr]
		
		for key in sorted(transcript_data.keys()):
			print 'feature: ' + key + '# of transcripts: ' + str(len(transcript_data[key]))
			
			if len(transcript_data[key]) < 1:
				continue
			
			for track in transcript_data[key]:
				print 'adding ' + key + ' to file'
#				print "|".join(track.values())
				notes = track.get('notes', ".")
				newfile.write("\t".join([chr,'rtracklayer_'+key,'sequence_feature',track['start'], track['stop'], track['score'],track['strand'],'.',notes]))
				newfile.write("\n")
	#	sys.exit

# def _make_bed_file(dict_to_print, outfile):
#

def _get_transcripts(file_data, genome_data, filter_type, filter_val, feat_type, overlap):
	transcripts = list()
	chrom_plus = list()
	chrom_minus = list()
	
	newdata = dict()
	unfiltered_data = dict()
	unmatched_trans = dict()
	remove_list = dict()

	for chromosome in genome_data.keys():
		if chromosome == 'chrmt':  # skip mito chromosome
			print 'skipping mito'
			continue
		
		matching_file_data = file_data[chromosome]
	#	unmatched_trans = matching_file_data	
#		print 'chromosome:' + chromosome + ':' + matching_file_data[0]['chr']
		#sys.exit()
		
		## split by strand ##
		for transcript in matching_file_data:
			if transcript['strand'] == '+':
				chrom_plus.append(transcript)
			else:
				chrom_minus.append(transcript)
		
		print '# transcripts to search: ' + str(len(matching_file_data))
		print '# plus strand: ' + str(len(chrom_plus))
		print '# minus strand: ' + str(len(chrom_minus))
		
		if feat_type != 'all':  ## use specific feature type if specified
		
			for element in genome_data[chromosome]:
#				print feat_type + ':' + element['feat_type']
				if element['feat_type'] != feat_type:
#					print 'feature types don\'t match'
					continue
				else: #matching feature types
					data_to_search = chrom_plus #default is plus strand
				
					if element['strand'] == '-':
						data_to_search = chrom_minus
				
					transcripts = _find_overlapping_transcripts(element, data_to_search, overlap)  # one feature, finding overlaps
					
# 					print 'number of transcripts: ' + str(len(transcripts))
# 					if len(transcripts) > 0:
# 						print '|'.join([d['key'] for d in transcripts])
# 						sys.exit()
#				print element['feat_type'] + "->" + chromosome + " " + str(element['start']) + " to " + str(element['stop']) + ":" + element['notes'] + " has " + str(len(transcripts)) + ' number of matches\n'
# 					notes_array = element['notes'].split(';')

					feat_id = element['feat_name']

# 					print "feat:" + feat_id
					if len(transcripts) == 0: # skip if no matches
						continue
						
					## add to unfiltered list and Add 'key' to remove from list
					if chromosome in unfiltered_data.keys():
							unfiltered_data[chromosome][feat_id] = transcripts
							remove_list[chromosome].extend([d['id'].rstrip() for d in transcripts])
					else:
							unfiltered_data[chromosome] = dict()
							unfiltered_data[chromosome][feat_id] = transcripts
							remove_list[chromosome] = [d['id'].rstrip() for d in transcripts]
					
					## filter transcripts by score or number 
					if filter_type == 'count':  # filter by count
						max_index = int(filter_val) - 1
						if chromosome in newdata.keys():
							newdata[chromosome][feat_id] = transcripts[0:filter_val]
						else:
							newdata[chromosome] = dict()
							newdata[chromosome][feat_id] = transcripts[0:filter_val]
					else:  # filter by score
						for each in transcripts:
							if each['score'] >= filter_val: # if score is > or = to cut off, add to array
								newdata[chromosome][feat_id].append(each)
					
# 					print "transcripts: " 
# 					for item in newdata[chromosome][feat_id]:
# 						print "keys: " + "|".join(item.keys())
# 						print 'vals: ' + ','.join(item.values())
#					sys.exit()
					
					
## take matched transcripts out of original set; use 'id' to find duplicates between original list and matched list
####################################################################################################################
# remove matching transcripts from original chromosome list #

# a = [x for x in a if x['link'] not in b]
## unmatched transcripts for chrV: 19983 ##
	print "|".join(file_data.keys()) 
	print "keys for remove list: " + ":".join(remove_list.keys())
	
	all_ids_to_remove = list()
	for one in remove_list.keys():
		all_ids_to_remove.extend(remove_list[one])

	for each_chr in remove_list.keys():
		print "removing matching transcripts from " + each_chr
		print "# of matches (length of remove list): " + str(len(remove_list[each_chr]))
		print "original # of transcripts (length of original file with transcripts on that chr): " + str(len(file_data[each_chr]))
		
		unmatched_trans[each_chr] = dict()
		unmatched_trans[each_chr][each_chr + "_non_ORF"] = [x for x in file_data[each_chr] if x['id'] not in all_ids_to_remove]
# 	
	return (unfiltered_data, newdata, unmatched_trans)
	
def _calc_overlap(trans_start, trans_stop, feat_start, feat_stop):
	feat_length = feat_stop - feat_start
	trans_length = trans_stop - trans_start
	
	if trans_start > feat_start and trans_stop > feat_stop:  # transcript starts inside the feature and ends outside the feature
		return (feat_stop - trans_start)/feat_length
	elif trans_start <= feat_start and trans_stop < feat_stop: # transcript starts outside of feature and ends inside feature
		return (trans_stop - feat_start)/feat_length
	elif trans_start > feat_start and trans_stop < feat_stop: # transcript starts and ends within the feature
		return (trans_length/trans_stop)
	else:  # transcript completely overlaps the feature -- trans_start <= feat_start AND trans_stop >= feat_stop
		return 1
	
def _find_overlapping_transcripts(feat_element, trans_data, overlap):
	
	feat_start = int(feat_element['start'])
	feat_stop = int(feat_element['stop'])
	strand = feat_element['strand']
	notes = feat_element['notes']
	feat_type = feat_element['feat_type']
	feat_name = feat_element['feat_name'] #notes_array[0].replace("ID=","")
	
	## calculate the range for start and stops to cover the overlap of the ORF
	feat_length = int(feat_stop - feat_start)
	feat_overlap = int(float(feat_length) * float(overlap))
	min_feat_start = int(feat_start - feat_overlap)
	min_feat_stop = int(feat_start + feat_overlap) # start plus the percent of overlap
	max_feat_start = int(feat_stop - feat_overlap) # distance from stop that will overlap req. amount
	
#	print 'finding transcripts for ' + feat_start + " to " + feat_end + ", " + strand + ':' + notes + '\n'
	
	
	filtered_data = list()
	sort_data = list()
	match_list = list()
	print "start # of transcipts: "+ str(len(trans_data))
	slice = 0
	
	for each in trans_data:
#		print "start transcript chr " + each['chr'] + ' and strand: ' + each['strand'] + " feat: "+ feat_element['chr'] + ',' + feat_element['strand']

		if each['chr'] != feat_element['chr'] or each['strand'] != feat_element['strand']:
			continue
		
		trans_start = int(each['start'])
		trans_stop = int(each ['stop'])
		# transcript length
		trans_length = trans_stop - trans_start
		
# 		print "f keys:" + ",".join(feat_element.keys())
# 		print "f values:" + "|".join(feat_element.values())
# 		print "t keys:" + ",".join(each.keys())
# 		print "t values: "+ "|".join(each.values())
# 		
		# skip if transcript length is less than overlap requirement of the ORF length
		if trans_length < feat_overlap:
			continue
		
		# if it is greater than or equal to overlap length, then see if it overlaps the ORF the correct amount then add to the 
		# 1. transcript_start is less than or eq to ORF start and trans_stop is greater than or eq to overlap req of the ORF (5' overlap AND full coverage)
		# 2. transcript start is greater than ORF start, but less than max_feat_start and transcript stop is greater than ORF stop (3' overlap)
		# 3. transcript start is greater than ORF start and less than max_feat_start AND transcript stop is less than ORF stop 
		# (get the diff between trans and ORF start and then add that to the min_feat_stop to find new min_feat_stop)
		
# 		print "transcript: " + str(trans_start) + " to " + str(trans_stop) + "; length: "+ str(trans_length)
# 		print "gene: " + str(feat_start) + " to " + str(feat_stop) + "; " + str(feat_overlap) + "bp overlap"
# 		
		
		if (trans_start <= feat_start and min_feat_stop <= trans_stop) or (feat_start <= trans_start <= max_feat_start and trans_stop > feat_stop) or (feat_start <= trans_start <= max_feat_start and (min_feat_stop + (trans_start - feat_start)) <= trans_stop):
		
		# calculate the % overlap #
			each['per_overlap'] = str(_calc_overlap(trans_start, trans_stop, feat_start, feat_stop))
# 			print 'match for ' + feat_name + ": feat coord - " + str(feat_start) + "," + str(feat_stop) + " and match: " + str(trans_start) +"-"+ str(trans_stop)
# 			print 'chr:' + each['chr'] + ' and strand ' + each['strand'] + ' match feature ' + feat_element['chr'] + ' and strand ' + feat_element['strand']
			each['notes'] = notes # add notes
			each['feat'] = feat_name
			print each['id'] + ' overlaps ' + each['feat']
			filtered_data.append(each) # add to filtered data array
			# add to list of list elements to remove from list at the end
	#		match_list.append(slice)
	#		sys.exit()
		else: 
			continue
#			print 'no match'
#			sys.exit()
	#	slice = slice + 1 # increment slice
#	print '# transcripts before sort: ' + str(len(filtered_data))
	# listsorted = sorted(XWordDict, key=lambda x: int(operator.itemgetter("pos")(x)))
	#sort_data = sorted(filtered_data, int(key=itemgetter('score')), reverse=True)
	#sort_data = sorted(filtered_data, key = lambda x: int(operator.itemgetter("score")(x)), reverse=True)
	#listsorted = sorted(XWordDict, key=lambda x: int(x['pos']))	
	
	sort_data = sorted(filtered_data, key=lambda x: (int(x['score']), int(x['per_overlap'])), reverse=True)  # sort by score and then by % overlap
	
# 	if len(sort_data) > 0:
# 		for obj in sort_data:
# 			print "keys: " + "|".join(obj.keys())
# 			print "vals: " + ",".join(obj.values())
# 		sys.exit()
		
	return sort_data
		
def _parse_sac_gff(file):
	file_obj = open(file, 'r')
#	file_data = open(file).readlines()
	col_count = 0

	ordered_keys = list()
	data_dict = dict()
	feat_data = list()
	headers = dict()
	
	feat_count = 1
	
	for line in file_obj:
	# skip line if it is a comment
		if (re.match("^#+", line)):
#			print "comment line: " + line
			continue

		line_array = list()
#		print "LINE:"+ line + "\n"
		
		#line_array = line.split(",") # split by commas
		
		length = len(line_array)  # check how many columns there are
		
		if len(line_array) <= 1:  # if it is tab-delimited, then split by \t
			line_array = line.split("\t")
			
			if len(line_array) <= 1: ## skip any line that can't be parsed
				continue
				
		key = line_array[0]
		
		notes_array = line_array[8].rstrip().split(';')
		feat_name = notes_array[0].replace("ID=","")

#		print "data: " + line

		# if key exists already, then add to data array
		if key in data_dict.keys(): 
			#feat_count = feat_count + 1 ## increase feat_count number
			data_dict[key].append({'feat_name': feat_name, 'chr': line_array[0],'start' : line_array[3], 'stop': line_array[4], 'strand' : line_array[6], 'feat_type': line_array[2],'notes': line_array[8].rstrip()})
		
		else:
		#	print "new chromosome: " + line_array[0]
		## else make a new key
			data_dict[line_array[0]] = [{'feat_name': feat_name, 'chr':line_array[0],'start' : line_array[3], 'stop': line_array[4], 'strand' : line_array[6], 'feat_type': line_array[2], 'notes':line_array[8].rstrip()}]
			
		# add to feat_count
# 	for chr in data_dict.keys():
# 		print "chr:" + chr
# 		for row in data_dict[chr]:
# 			if row['feat_type'] == 'gene':
		
			
	return (data_dict)

def _open_make_hash(file):
	file_obj = open(file, 'r')
#	file_data = open(file).readlines()
	col_count = 0

	ordered_keys = list()
	data_dict = dict()
	feat_data = list()
	headers = dict()
	
	for line in file_obj:
	# skip line if it is a comment
		if (re.match("^#+", line)):
#			print "comment line: " + line
			continue

		line_array = list()
#		print "LINE:"+ line + "\n"
		
		line_array = line.split(",") # split by commas
		
		length = len(line_array)  # check how many columns there are
		
		if len(line_array) <= 1:  # if it is tab-delimited, then split by \t
#			print "data: " + line
			line_array = line.split("\t")
			
		key = line_array[0]
		
		# if key exists already, then add to data array
		if key in data_dict.keys(): 
#			feat_count = feat_count + 1 ## increase feat_count number
			data_dict[key].append({'id': "".join(line_array[0:5]), 'chr': line_array[0], 'start' : line_array[3], 'stop': line_array[4], 'strand' : line_array[6], 'score': line_array[5].rstrip()})
		
		else:
#			print "new chromosome: " + line_array[0]
		## else make a new key
			data_dict[line_array[0]] = [{'id': "".join(line_array[0:5]),'chr': line_array[0], 'start' : line_array[3], 'stop': line_array[4], 'strand' : line_array[6], 'score': line_array[5].rstrip()}]
			
		# add to feat_count
#	for key in data_dict.keys():
#		print key + ": num transcripts = " + str(len(data_dict[key]))#
#	sys.exit()
	return (data_dict)


if __name__ == "__main__":
   main(sys.argv[1:])