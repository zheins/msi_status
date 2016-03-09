import sys
import os
import csv

CANCER_TYPES = ['Colorectal Cancer', 'Anal Cancer']
POLE_HOT_SPOTS = ['P286R','P286H','P286S','V411L']
def tumorBreakdown(ccrType,sample_data,pt2sample):
	cancer_breakdown_output = ['\t'.join(['PATIENT_ID','SAMPLE_ID','NO_MUTAIONS','RATIO (INDELS/TOTAL_EVENTS)','NO_POLE_MUTATIONS','POLE_ANNOTATIONS','NO_POLE_HOTSPOTS','MSI_STATUS','GENE_PANEL','OTHER_PRO_ANNOTATIONS'])]
	cancer_breakdown_data = []

	#for each patient,sampleList grab variant data for input cancer type
	for ptid,sampleList in pt2sample.items():

		for sid in sampleList:
			data =  [data for sample,data in sample_data.iteritems() if sample == sid]

			for d in data:
				if d.get('CANCER_TYPE') == ccrType:
					varData = d.get('VARIANT',[])
					totalEvents = 0
					indelCount = 0
					poleCount = 0
					poleHotSpots = 0
					poleAnnotations = []
					protAnnotations = []
					
					#grab variant data and annotation information
					for vd in varData:
						totalEvents += 1
						if vd.get('Variant_Type') == 'INS' or vd.get('Variant_Type') == 'DEL':
							indelCount += 1
						
						if vd.get('Hugo_Symbol') == 'POLE':
							poleCount += 1
							poleAnnot = vd.get('Amino_Acid_Change').replace('p.','')
							poleAnnotations.append(poleAnnot)
							if poleAnnot in POLE_HOT_SPOTS:
								poleHotSpots += 1
						else:
							aaChange = vd.get('Amino_Acid_Change').replace('p.','')
							protAnnot = vd.get('Hugo_Symbol')+':'+aaChange
							protAnnotations.append(protAnnot)

					msiStatus = d.get('MSI_STATUS','')
					genePanel = d.get('GENE_PANEL')
					poleAnnotations = '; '.join(poleAnnotations)
					protAnnotations = '; '.join(sorted(protAnnotations))


					if totalEvents == 0:
						ratio = 0.00

					else:
						ratio = round(1.0*indelCount/totalEvents,2)


					line = [ptid,sid,str(len(varData)),str(ratio),str(poleCount),poleAnnotations,str(poleHotSpots),msiStatus,genePanel,protAnnotations]
					
					cancer_breakdown_data.append('\t'.join(line))

	cancer_breakdown_output.extend(list(set(cancer_breakdown_data)))
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
	# process general clin file
	clin_reader = csv.DictReader(clin_file, dialect = 'excel-tab')
	for line in clin_reader:
		if line['CANCER_TYPE'] in CANCER_TYPES:
			sample_data[line['SAMPLE_ID']] = {'CANCER_TYPE':line['CANCER_TYPE'], 'CANCER_TYPE_DETAILED':line['CANCER_TYPE_DETAILED'], 'GENE_PANEL':line['GENE_PANEL']}
			ptId = line['PATIENT_ID']
			sidList = pt2sample.get(ptId,[])
			sidList.append(line['SAMPLE_ID'])
            
			pt2sample[ptId] = list(set(sidList))

	# process supplemental clin file
	clin_detailed_reader = csv.DictReader(clin_detailed_file, dialect = 'excel-tab')
	for line in clin_detailed_reader:
		if line['SAMPLE_ID'] in sample_data:
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
			variant['Amino_Acid_Change'] = line['Amino_Acid_Change']
			if sample_data[sid].get('VARIANT', '')  == '':
				sample_data[sid]['VARIANT'] = [variant]
			else:
				sample_data[sid]['VARIANT'].append(variant)

	for ccrType in CANCER_TYPES:
		cancer_breakdown_output = tumorBreakdown(ccrType,sample_data,pt2sample)

		fname = 'cancer_breakdown_'+ccrType+'.txt'
		fh = open(fname,'w')
		fh.write(cancer_breakdown_output)
		fh.close()				

	# print len(sample_data)
	# for sample, data in sample_data.iteritems():
	# 	print sample
	# 	print data

if __name__ == '__main__':
	main()