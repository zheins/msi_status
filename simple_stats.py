import sys
import os
import csv

CANCER_TYPES = ['Colorectal Cancer', 'Anal Cancer','Endometrial Cancer']
POLE_HOT_SPOTS = ['P286R','P286H','P286S','V411L']
POLE_CURATED_MUTS = ['S297F','S297Y','F367S','P436R','M444K','S459F','A456P','L424V','L424I']
#MSI_STATUS_MAP = {'MSI':'Yes','High':'Yes','Stable':'No','MSS':'No','':'Other',' ':'Other','Pending':'Other','Inconclusive':'Other','NA':'Other'}

def get_hotspots():
	if os.path.exists('hotspots_data.txt'):
		filedata = open('hotspots_data.txt').read().split('\n')[1:]
	else:
		print 'hotspots_data.txt not in current working directory'
		sys.exit(2)

	hotspots = {}

	for line in filedata:
		data = line.split('\t')
		gene = data[0]
		aa_changes = data[1].split(';')

		hotspots[gene] = aa_changes
	return hotspots


def tumor_breakdown(cancer_type, sample_data, pt2sample):
	cancer_breakdown_output = ['\t'.join(['PATIENT_ID','SAMPLE_ID','CANCER_TYPE_DETAILED','NO_MUTATIONS','RATIO (INDELS/TOTAL_EVENTS)','POLE_AA_CHANGES','NUM_POLE_MUTATIONS','NUM_POLE_HOTSPOTS','NUM_CURATED_POLE_MUTS','NUM_OTHER_HOTSPOTS','HOT_SPOTS_PRESENT','MSI_STATUS','GENE_PANEL','OTHER_AA_CHANGES'])]
	cancer_breakdown_data = []

	#retrieve hotspots
	hotspots = get_hotspots()

	#for each patient,sampleList grab variant data for input cancer type
	for ptid,sampleList in pt2sample.items():
		for sid in sampleList:
			data =  [data for sample,data in sample_data.iteritems() if sample == sid]

			for d in data:
				if d.get('CANCER_TYPE') == cancer_type:
					variant_data = d.get('VARIANT',[])
					cancer_type_detailed = d.get('CANCER_TYPE_DETAILED')					
					gene_panel = d.get('GENE_PANEL')
					msi_status = d.get('MSI_STATUS','')
					
					num_mutations = str(len(variant_data))
					
					total_events = 0
					indel_count = 0
					pole_count = 0
					pole_hotspots = 0
					hotspot_count = 0
					curated_count = 0
					pole_annotations = []
					hotspot_annotations = []
					other_annotations = []
					
					#grab variant data and annotation information
					for vd in variant_data:
						gene = vd.get('Hugo_Symbol')
						total_events += 1
						
						#count indel events
						if vd.get('Variant_Type') == 'INS' or vd.get('Variant_Type') == 'DEL':
							indel_count += 1
						
						#increment appropriate counts and append amino acid changes to appropriate annotation list
						if gene == 'POLE':
							pole_count += 1
							poleAnnot = vd.get('Amino_Acid_Change').replace('p.','')
							pole_annotations.append(poleAnnot)
							
							if poleAnnot in POLE_HOT_SPOTS:
								pole_hotspots += 1
							elif poleAnnot in POLE_CURATED_MUTS:
								curated_count += 1
						
						elif gene in hotspots.keys():
							aaChange = vd.get('Amino_Acid_Change').replace('p.','')
						
							if aaChange in hotspots.get(gene):
								hotspot_count += 1
								hsAnnot = gene+':'+aaChange
								hotspot_annotations.append(hsAnnot)

						else:
							aaChange = vd.get('Amino_Acid_Change').replace('p.','')
							protAnnot = gene+':'+aaChange
							other_annotations.append(protAnnot)
					
					#format counts and annotations for output
					pole_count = str(pole_count)
					pole_hotspots = str(pole_hotspots)
					curated_count = str(curated_count)
					hotspot_count = str(hotspot_count)
					pole_annotations = '; '.join(pole_annotations)
					other_annotations = '; '.join(sorted(other_annotations))
					hotspot_annotations = '; '.join(sorted(hotspot_annotations))

					#calculate ratio of indels to total events
					if total_events == 0:
						ratio = str(0.00)
					else:
						ratio = str(round(1.0*indel_count/total_events,2))
																		
					#normalize msi status
					#norm_msi_status = MSI_STATUS_MAP.get(msi_status,'Invalid MSI value')

					#data to append for output file data
					line = [ptid,sid,cancer_type_detailed,num_mutations,ratio,pole_annotations,pole_count,pole_hotspots,curated_count,hotspot_count,hotspot_annotations,msi_status,gene_panel,other_annotations]	
					cancer_breakdown_data.append('\t'.join(line))

	#sort output by patient ID and get rid of any possibly duplicated data
	cancer_breakdown_output.extend(sorted(list(set(cancer_breakdown_data))))
	return '\n'.join(cancer_breakdown_output)

def main():
	clin_filename = sys.argv[1]
	clin_detailed_filename = sys.argv[2]
	maf_filename = sys.argv[3]

	if os.path.exists(clin_filename) and os.path.exists(maf_filename) and os.path.exists(clin_detailed_filename):
		clin_file = open(clin_filename, 'rU')
		clin_detailed_file = open(clin_detailed_filename, 'rU')
		maf_file = open(maf_filename, 'rU')
	else:
		print 'bad filenames'
		sys.exit(2)

	sample_data = {}
	pt2sample = {}
	cancer_type = ' '.join(list(map(lambda v: v.capitalize(),' '.join(clin_detailed_filename.replace('.txt','').split('_')[-2:]).split())))
	
	# process general clin file
	clin_reader = csv.DictReader(clin_file, dialect = 'excel-tab')
	for line in clin_reader:
		if line['CANCER_TYPE'] in CANCER_TYPES:
			sample_data[line['SAMPLE_ID']] = {'CANCER_TYPE':line['CANCER_TYPE'], 'CANCER_TYPE_DETAILED':line['CANCER_TYPE_DETAILED'], 'GENE_PANEL':line['GENE_PANEL']}
			ptid = line['PATIENT_ID']
			sid_list = pt2sample.get(ptid,[])
			sid_list.append(line['SAMPLE_ID'])
            
			pt2sample[ptid] = list(set(sid_list))

	# process supplemental clin file
	clin_detailed_reader = csv.DictReader(clin_detailed_file, dialect = 'excel-tab')
	for line in clin_detailed_reader:
		if line['SAMPLE_ID'] in sample_data:
			if line.get('MSI_STATUS','') != '':
				sample_data[line['SAMPLE_ID']]['MSI_STATUS'] = line['MSI_STATUS']


	# process maf
	maf_reader = csv.DictReader(maf_file, dialect = 'excel-tab')
	for line in maf_reader:
		# construct a variant dictionary with fields of interest from maf, then add it to the sample data sample dictionary variant list
		variant = {}
		sid = line['Tumor_Sample_Barcode']
		if sid in sample_data:
			variant['Variant_Type'] = line['Variant_Type']
			variant['Hugo_Symbol'] = line['Hugo_Symbol']
			variant['Amino_Acid_Change'] = line['HGVSp_Short']
			if sample_data[sid].get('VARIANT', '')  == '':
				sample_data[sid]['VARIANT'] = [variant]
			else:
				sample_data[sid]['VARIANT'].append(variant)


	cancer_breakdown_output = tumor_breakdown(cancer_type, sample_data, pt2sample)
	filename = 'cancer_breakdown_'+cancer_type.replace(' ','_')+'.txt'
	fh = open(filename,'w')
	fh.write(cancer_breakdown_output)
	fh.close()				

	# print len(sample_data)
	# for sample, data in sample_data.iteritems():
	# 	print sample
	# 	print data

if __name__ == '__main__':
	main()