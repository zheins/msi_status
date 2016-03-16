import sys
import os
import csv
import re

INFRAME = ['In_Frame_Ins', 'In_Frame_Del']
FRAMESHIFT = ['frameshift_deletion', 'Frame_Shift_Ins', 'Frame_Shift_Del']

patient_data = {}
sample_to_patient = {}

def process_clinical_file(clin_file):
    clin_reader = csv.DictReader(clin_file, dialect = 'excel-tab')
    for line in clin_reader:
        patient = patient_data.get(line['PATIENT_ID'],{})
        sample_to_patient[line['SAMPLE_ID']] = line['PATIENT_ID']
        try:
            patient['CANCER_TYPE'].append(line['CANCER_TYPE'])
            patient['CANCER_TYPE_DETAILED'].append(line['CANCER_TYPE_DETAILED'])
            patient['SAMPLES'].append(line['SAMPLE_ID'] + ':' + line['SAMPLE_TYPE']) 
        except KeyError:
            patient['CANCER_TYPE'] = [line['CANCER_TYPE']]
            patient['CANCER_TYPE_DETAILED'] = [line['CANCER_TYPE_DETAILED']]
            patient['SAMPLE_TYPE'] = [line['SAMPLE_TYPE']]
            patient['SAMPLES'] = [line['SAMPLE_ID'] + ':' + line['SAMPLE_TYPE']]
        patient_data[line['PATIENT_ID']] = patient

def process_maf_file(maf_file):
    maf_reader = csv.DictReader(maf_file, dialect = 'excel-tab')
    for line in maf_reader:
        # construct a variant dictionary with fields of interest from maf, then add it to the sample data sample dictionary variant list
        variant = {}
        try:
            pid = sample_to_patient[line['Tumor_Sample_Barcode']]
        except KeyError:
            continue
        variant_type = line['Variant_Classification']
        if variant_type in INFRAME:
            variant_type = 'inframe'
        elif variant_type in FRAMESHIFT:
            variant_type = 'frameshift'
            
        variant['Variant_Type'] = variant_type 
        variant['Hugo_Symbol'] = line['Hugo_Symbol']
        variant['SAMPLE_ID'] = line['Tumor_Sample_Barcode']

        aa_change = line['HGVSp_Short']
        if len(aa_change.strip()) > 0 and aa_change != 'p.=':
            m = re.search('p.([A-Z*]*)([0-9]*)(\S*)', aa_change)
            ref = m.group(1)
            loc = m.group(2)
            alt = m.group(3)

            variant['ref'] = ref
            variant['loc'] = loc
            variant['alt'] = alt
       
        print aa_change
        print variant.get('loc',999)
        print variant.get('Hugo_Symbol')
        print variant.get('SAMPLE_ID')
        print line['Start_Position']
        if int(variant.get('loc',999)) < 350 and variant_type != 'inframe' and variant_type != 'frameshift':
            continue

        if patient_data[pid].get('VARIANTS', '')  == '':
            patient_data[pid]['VARIANTS'] = [variant]
        else:
            patient_data[pid]['VARIANTS'].append(variant)

def main():
    clin_filename = sys.argv[1]
    maf_filename = sys.argv[2]

    if os.path.exists(clin_filename) and os.path.exists(maf_filename):
        clin_file = open(clin_filename, 'rU')
        maf_file = open(maf_filename, 'rU')
    else:
        print 'bad filenames'
        sys.exit(2)

	
    process_clinical_file(clin_file)
    process_maf_file(maf_file)

    for patient, data in patient_data.iteritems():
        print patient, data
        print
        


if __name__ == '__main__':
	main()
