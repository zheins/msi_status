import sys
import os
import csv

DETAILED_HEADERS = ['Hugo', 'Amino_Acid', 'Site', 'Total_Events', 'Total_Insertions', 'Total_Deletions', 'Total_Truncating_Events', 'Missense_Count', 'Nonsense_Count', 'Inframe_Count', 'Frameshift_Count', 'Inframe_INS_Count', 'Inframe_DEL_Count', 'Frameshift_INS_Count', 'Frameshift_DEL_Count', 'Truncating_Frameshift_INS_Count', 'Truncating_Frameshift_DEL_Count', 'Truncating_Frameshift_Count', 'Tumor_Type_List', 'Tumor_Type_Count', 'Primary_Tumor_list', 'Primary_Tumor_Count', 'Metastatic_Tumor_List', 'Metastatic_Tumor_Count', 'OQL_All', 'OQL_Primary', 'OQL_Metastatic', 'Cohort_All', 'Cohort_Primary', 'Cohort_Metastatic']
INT_TYPE_HEADERS = ['Total_Events', 'Total_Insertions', 'Total_Deletions', 'Total_Truncating_Events', 'Missense_Count', 'Nonsense_Count', 'Inframe_Count', 'Frameshift_Count', 'Inframe_INS_Count', 'Inframe_DEL_Count', 'Frameshift_INS_Count', 'Frameshift_DEL_Count', 'Truncating_Frameshift_INS_Count', 'Truncating_Frameshift_DEL_Count', 'Truncating_Frameshift_Count', 'Tumor_Type_Count', 'Primary_Tumor_Count', 'Metastatic_Tumor_Count']
LIST_TYPE_HEADERS = ['OQL_All', 'OQL_Primary', 'OQL_Metastatic', 'Cohort_All', 'Cohort_Primary', 'Cohort_Metastatic']
TUMOR_LIST_HEADERS = ['Tumor_Type_List', 'Primary_Tumor_list', 'Metastatic_Tumor_List']
def merge_studies(fname_study1, fname_study2):
	common_sites = get_common_sites(fname_study1, fname_study2)

	study1_file = open(fname_study1, 'rU')
	study1_reader = csv.DictReader(study1_file, dialect = 'excel-tab')
	merged_file_output = ['\t'.join(DETAILED_HEADERS)]
	
	merged_lines = {}
	#merge lines with common sites
	for line1 in study1_reader:
		print line1
		if line1['Site'] in common_sites:
			data_to_merge = {}
			for key in DETAILED_HEADERS:
				if key in INT_TYPE_HEADERS:
					data_to_merge[key] = int(line1[key])
				elif key in LIST_TYPE_HEADERS:
					data_to_merge[key] = line1[key].split(' ')
				else:
					data_to_merge[key] = line1[key]
			merged_lines[line1['Site']] = data_to_merge
		else:
			merged_line = [line1.get(key) for key in DETAILED_HEADERS]
			merged_file_output.append('\t'.join(merged_line))	
			print merged_line
	
	study2_file = open(fname_study2, 'rU')	
	study2_reader = csv.DictReader(study2_file, dialect = 'excel-tab')

	for line2 in study2_reader:
		print line2
		if line2['Site'] in common_sites:
			data_to_merge = merged_lines.get(line2['Site'])
			for key in DETAILED_HEADERS:
				if key in INT_TYPE_HEADERS:
					data_to_merge[key] = str(data_to_merge.get(key) + int(line2[key]))
				elif key in LIST_TYPE_HEADERS:
					new_list = data_to_merge.get(key)
					new_list.extend(line2[key].split(' '))
					data_to_merge[key] = ' '.join(list(set(new_list)))
				elif key in TUMOR_LIST_HEADERS:
					data_to_merge[key] = '; '.join([data_to_merge.get(key),line2[key]])
			merged_lines[line2['Site']] = data_to_merge
		else:
			merged_line = [line2.get(key) for key in DETAILED_HEADERS]
			merged_file_output.append('\t'.join(merged_line))	
			print merged_line


	for site,data in merged_lines.iteritems():
		merged_line = []
		for key in DETAILED_HEADERS:
			merged_line.append(data.get(key))
		merged_file_output.append('\t'.join(merged_line))

	return merged_file_output


def get_common_sites(fname_study1, fname_study2):
	study1_file = open(fname_study1, 'rU')
	study2_file = open(fname_study2, 'rU')	

	study1_reader = csv.DictReader(study1_file, dialect = 'excel-tab')
	study2_reader = csv.DictReader(study2_file, dialect = 'excel-tab')
	
	study1_sites = [line['Site'] for line in study1_reader]
	study2_sites = [line['Site'] for line in study2_reader]

	common_sites = [site for site in study2_sites if site in study1_sites]
	return common_sites	

def main ():
	try:
	    fname_study1 = sys.argv[1]
	    fname_study2 = sys.argv[2]
	except IndexError:
		print 'Not enough input arguments.'
		sys.exit(2)

	if not os.path.exists(fname_study1) or not os.path.exists(fname_study2):
	# 	study1_file = open(fname_study1, 'rU')
	# 	study2_file = open(fname_study2, 'rU')
	# else:
		print 'bad filenames'
		sys.exit(2)	

	cancer1 = os.path.basename(fname_study1).split('-')[0]
	cancer2 = os.path.basename(fname_study2).split('-')[0]
	d = os.path.dirname(fname_study1)
	fname_merged = d + '-'.join(['/MERGED',cancer1,cancer2,'Detailed_Summary']) + '.txt'



	merged_output = merge_studies(fname_study1, fname_study2)

	fh = open(fname_merged,'w')
	fh.write('\n'.join(merged_output))
	fh.close()

	print 'Merged study written to: ', fname_merged



if __name__ == '__main__':
	main()