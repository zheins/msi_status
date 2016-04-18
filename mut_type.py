import sys
import os
import csv
import re
import time
import datetime

VARIANT_CLASSIFICATIONS = ['Nonsense_Mutation', 'Missense_Mutation', 'In_Frame_Ins', 'In_Frame_Del', 'frameshift_deletion', 'Frame_Shift_Ins', 'Frame_Shift_Del']
DETAILED_HEADERS = ['Hugo', 'Amino_Acid', 'Site', 'Total_Events', 'Total_Insertions', 'Total_Deletions', 'Total_Truncating_Events', 'Missense_Count', 'Nonsense_Count', 'Inframe_Count', 'Frameshift_Count', 'Inframe_INS_Count', 'Inframe_DEL_Count', 'Frameshift_INS_Count', 'Frameshift_DEL_Count', 'Truncating_Frameshift_INS_Count', 'Truncating_Frameshift_DEL_Count', 'Truncating_Frameshift_Count', 'Tumor_Type_List', 'Tumor_Type_Count', 'Primary_Tumor_list', 'Primary_Tumor_Count', 'Metastatic_Tumor_List', 'Metastatic_Tumor_Count', 'OQL_All', 'OQL_Primary', 'OQL_Metastatic', 'Cohort_All', 'Cohort_Primary', 'Cohort_Metastatic']
GENES = ['TP53']
ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H%M')

patient_data = {}
sample_to_patient = {}
amino_acid_data = {}
truncating_events = {}
cancer_type_list = set()


def process_clinical_file(clin_file):
    clin_reader = csv.DictReader(clin_file, dialect = 'excel-tab')
    for line in clin_reader:
        patient = patient_data.get(line['PATIENT_ID'],{})
        sample_to_patient[line['SAMPLE_ID']] = line['PATIENT_ID']
        samples = {}
        try:
            patient['CANCER_TYPE'].append(line['CANCER_TYPE'])
            patient['CANCER_TYPE_DETAILED'].append(line['CANCER_TYPE_DETAILED'])
            patient['SAMPLE_TYPE'].append(line['SAMPLE_TYPE'])
        except KeyError:
            patient['CANCER_TYPE'] = [line['CANCER_TYPE']]
            patient['CANCER_TYPE_DETAILED'] = [line['CANCER_TYPE_DETAILED']]
            patient['SAMPLE_TYPE'] = [line['SAMPLE_TYPE']]
        patient_data[line['PATIENT_ID']] = patient


        samples = {'SAMPLE_ID':line['SAMPLE_ID'], 'SAMPLE_TYPE':line['SAMPLE_TYPE'],'CANCER_TYPE':line['CANCER_TYPE']}
        if patient_data[line['PATIENT_ID']].get('SAMPLES','') == '':
            patient_data[line['PATIENT_ID']]['SAMPLES'] = [samples]
        else:
            patient_data[line['PATIENT_ID']]['SAMPLES'].append(samples)

        cancer_type_list.add(line['CANCER_TYPE'])


def process_maf_file(maf_file):
    maf_reader = csv.DictReader(maf_file, dialect = 'excel-tab')
    for line in maf_reader:
        # construct a variant dictionary with fields of interest from maf, then add it to the sample data sample dictionary variant list
        variant = {}
        try:
            pid = sample_to_patient[line['Tumor_Sample_Barcode']]
        except KeyError:
            continue
        
        if line['Hugo_Symbol'] not in GENES:
            continue
        variant['Hugo_Symbol'] = line['Hugo_Symbol']

        variant_class = line['Variant_Classification']
        if variant_class not in VARIANT_CLASSIFICATIONS:
            continue

        if variant_class == 'frameshift_deletion':
            variant_class = 'Frame_Shift_Del'
            
        variant['Variant_Class'] = variant_class 
        variant['SAMPLE_ID'] = line['Tumor_Sample_Barcode']
        variant['Variant_Type'] = line['Variant_Type']

        aa_change = line['HGVSp_Short']
        if len(aa_change.strip()) > 0 and aa_change != 'p.=':
            m = re.search('p.([A-Z*]*)([0-9]*)(\S*)', aa_change)
            ref = m.group(1)
            loc = m.group(2)
            alt = m.group(3)

            variant['ref'] = ref
            variant['loc'] = loc
            variant['alt'] = alt
       
        if int(variant.get('loc',999)) > 350 and alt.find('*') < 0:
            continue

        if patient_data[pid].get('VARIANTS', '')  == '':
            patient_data[pid]['VARIANTS'] = [variant]
        else:
            patient_data[pid]['VARIANTS'].append(variant)

def process_amino_acids():
    for patient, data in patient_data.iteritems():
        try:
            variant_data = data['VARIANTS']
        except KeyError:
            continue

        samples = data.get('SAMPLES')
        for var_data in variant_data:
            event_data = {}
            sid = var_data.get('SAMPLE_ID')

            key = (var_data.get('Hugo_Symbol'),var_data.get('ref'),var_data.get('loc'))
            
            var_class = var_data.get('Variant_Class')

            event_data['Variant_Type'] = var_data.get('Variant_Type')
            event_data['AA_Change'] = var_data.get('alt')

            sample_data = [v for i,v in enumerate(samples) if v.get('SAMPLE_ID') == sid][0]
            event_data['Sample_Type'] = sample_data.get('SAMPLE_TYPE')
            event_data['Cancer_Type'] = sample_data.get('CANCER_TYPE')
            event_data['SAMPLE_ID'] = sid 

            aa = amino_acid_data.get(key,{})
            try:
                aa[var_class].append(event_data)
            except KeyError:
                aa[var_class] = [event_data]

            amino_acid_data[key] = aa


def process_truncating_events():
    for aa,data in amino_acid_data.iteritems():
        for k,var_data in data.iteritems():
            for vd in var_data:
                if vd.get('AA_Change').find('*') >= 0:
                    vd['Amino_Acid'] = aa 

                    try:
                        truncating_events[k].append(vd)
                    except KeyError:
                        truncating_events[k] = [vd]


def format_tumor_counts(tumor_type_dict):
    formatted_tumor_type_str = []
    for k,v in tumor_type_dict.iteritems():
        s = k+' ('+str(v)+')' 
        formatted_tumor_type_str.append(s)
    return formatted_tumor_type_str


def filter_aa_data_set(aa_data_set, filter_type, filter):
    filtered_aa_data_set = []
    if aa_data_set != []:
        if filter_type == 'cancer':
            return [data for i,data in enumerate(aa_data_set) if data.get('Cancer_Type') == filter]            
        elif filter_type == 'sample':
            return [data for i,data in enumerate(aa_data_set) if data.get('Sample_Type') == filter]
        elif filter_type == '':
            return aa_data_set
    return filtered_aa_data_set


def process_aa_data(filter_type, filter):
    detailed_output = ['\t'.join(DETAILED_HEADERS)]

    for aa,aa_data_set in amino_acid_data.iteritems():
        det_output_data = ''
        
        #in frame data
        if_ins_data = filter_aa_data_set(aa_data_set.get('In_Frame_Ins',[]), filter_type, filter) 
        if_del_data = filter_aa_data_set(aa_data_set.get('In_Frame_Del',[]), filter_type, filter) 

        #frame shift data
        fs_ins_data = filter_aa_data_set(aa_data_set.get('Frame_Shift_Ins',[]), filter_type, filter) 
        fs_del_data = filter_aa_data_set(aa_data_set.get('Frame_Shift_Del',[]), filter_type, filter) 

        #missense and nonsense data
        missense_data = filter_aa_data_set(aa_data_set.get('Missense_Mutation',[]), filter_type, filter) 
        nonsense_data = filter_aa_data_set(aa_data_set.get('Nonsense_Mutation',[]), filter_type, filter) 

        #variables for counts/output
        truncating_frameshift_ins_count = 0
        truncating_frameshift_del_count = 0
        truncating_frameshift_count = 0

        missense_count = len(missense_data)
        nonsense_count = len(nonsense_data)
        
        inframe_ins_count = len(if_ins_data)
        inframe_del_count = len(if_del_data)
        inframe_count = inframe_ins_count + inframe_del_count

        frameshift_ins_count = len(fs_ins_data)
        frameshift_del_count = len(fs_del_data)
        frameshift_count = frameshift_ins_count + frameshift_del_count


        total_events = missense_count + nonsense_count + inframe_count + frameshift_count
        if total_events == 0:
            continue
        total_insertions = inframe_ins_count + frameshift_ins_count
        total_deletions = inframe_del_count + frameshift_del_count        

        #truncating data
        trunc_fs_ins_data = filter_aa_data_set(truncating_events.get('Frame_Shift_Ins',[]), filter_type, filter) 
        trunc_fs_del_data = filter_aa_data_set(truncating_events.get('Frame_Shift_Del',[]), filter_type, filter) 
        trunc_nonsense_data = filter_aa_data_set(truncating_events.get('Nonsense_Mutation',[]), filter_type, filter) 
        
        truncating_frameshift_ins_count = len([data for data in trunc_fs_ins_data if data.get('Amino_Acid') == aa])
        truncating_frameshift_del_count = len([data for data in trunc_fs_del_data if data.get('Amino_Acid') == aa]) 
        truncating_frameshift_count = truncating_frameshift_ins_count + truncating_frameshift_del_count

        truncating_nonsense_count = len([data for data in trunc_nonsense_data if data.get('Amino_Acid') ==  aa])
        total_truncating_events = truncating_frameshift_count + truncating_nonsense_count

        #tumor data
        tumor_type_list = {}
        primary_tumor_list = {}
        metastatic_tumor_list = {}

        #tumor data
        primary_tumor_count = 0
        metastatic_tumor_count = 0
        tumor_type_count = 0

        #oql lists
        oql_base = 'MUT = ' + aa[1] + str(aa[2])
        oql_all  = set()
        oql_primary = set()
        oql_metastatic = set()

        #cohort lists
        cohort_all = []
        cohort_primary = []
        cohort_metastatic = []

        for data_set in [if_ins_data,if_del_data,fs_ins_data,fs_del_data,missense_data,nonsense_data]:
            primary_count = 0
            met_count = 0
            for vd in data_set:
                tumor_type_list[vd.get('Cancer_Type')] = tumor_type_list.get(vd.get('Cancer_Type'),0) + 1
                oql_all.add(oql_base+vd.get('AA_Change'))
                cohort_all.append(vd.get('SAMPLE_ID'))

                if vd.get('Sample_Type') == 'Primary':
                    primary_tumor_list[vd.get('Cancer_Type')] = primary_tumor_list.get(vd.get('Cancer_Type'),0) + 1
                    primary_count += 1
                    oql_primary.add(oql_base+vd.get('AA_Change'))
                    cohort_primary.append(vd.get('SAMPLE_ID'))

                if vd.get('Sample_Type') == 'Metastasis':
                    metastatic_tumor_list[vd.get('Cancer_Type')] = metastatic_tumor_list.get(vd.get('Cancer_Type'),0) + 1
                    met_count += 1
                    oql_metastatic.add(oql_base+vd.get('AA_Change'))
                    cohort_metastatic.append(vd.get('SAMPLE_ID'))

            primary_tumor_count += primary_count 
            metastatic_tumor_count += met_count 
            tumor_type_count += primary_count + met_count

        #format tumor type lists
        tumor_type_list = '; '.join(format_tumor_counts(tumor_type_list))
        primary_tumor_list = '; '.join(format_tumor_counts(primary_tumor_list))
        metastatic_tumor_list = '; '.join(format_tumor_counts(metastatic_tumor_list))

        #format oql queries
        oql_all = ' '.join(list(oql_all))
        oql_primary = ' '.join(list(oql_primary))
        oql_metastatic = ' '.join(list(oql_metastatic))

        #format cohort lists
        cohort_all = ' '.join(cohort_all)
        cohort_primary = ' '.join(cohort_primary)
        cohort_metastatic = ' '.join(cohort_metastatic)

        #format data for output files
        gene_aa_site = '\t'.join([aa[0],aa[1],str(aa[2])])
        totals = '\t'.join(map(str,[total_events,total_insertions,total_deletions,total_truncating_events,missense_count,nonsense_count,inframe_count,frameshift_count]))
        in_frame_det_counts = '\t'.join(map(str,[inframe_ins_count,inframe_del_count]))
        frameshift_det_counts = '\t'.join(map(str,[frameshift_ins_count,frameshift_del_count]))
        trunc_fs_counts = '\t'.join(map(str,[truncating_frameshift_ins_count,truncating_frameshift_del_count,truncating_frameshift_count]))
        tumor_types = '\t'.join([tumor_type_list,str(tumor_type_count)])
        primary_tumors = '\t'.join([primary_tumor_list,str(primary_tumor_count)])
        metastatic_tumors = '\t'.join([metastatic_tumor_list,str(metastatic_tumor_count)])
        oql_queries = '\t'.join([oql_all,oql_primary,oql_metastatic])
        cohort_queries = '\t'.join([cohort_all,cohort_primary,cohort_metastatic])

        #put formatted data into detailed outputs
        det_output_data += '\t'.join([gene_aa_site,totals,in_frame_det_counts,frameshift_det_counts,trunc_fs_counts,tumor_types,primary_tumors,metastatic_tumors,oql_queries,cohort_queries])  
        detailed_output.append(det_output_data) 

    return detailed_output

def main():
    clin_filename = sys.argv[1]
    maf_filename = sys.argv[2]


    if os.path.exists(clin_filename) and os.path.exists(maf_filename):
        clin_file = open(clin_filename, 'rU')
        maf_file = open(maf_filename, 'rU')
    else:
        print 'bad filenames'
        sys.exit(2)

    print 'Processing clinical file'
    process_clinical_file(clin_file)
    print 'Processing MAF file'
    process_maf_file(maf_file)
    print 'Processing amino acid data'
    process_amino_acids()
    print 'Processing truncating events'
    process_truncating_events()
    print

    #generate detailed summary for all data
    print 'Generating comprehensive summary results.'
    detailed_output = process_aa_data('','')    
    fname_detailed = ts + '-Detailed_Comprehensive_Summary.txt'
    fh = open(fname_detailed,'w')
    fh.write('\n'.join(detailed_output))
    fh.close()
    print 'Detailed comprehensive summary written to: ', fname_detailed
    print

    #generate detailed summaries by cancer types
    print 'Data by cancer types written to dir:',ts+'-CancerBreakdown/'
    num_files = 0
    for cancer_type in cancer_type_list:
        detailed_output = process_aa_data('cancer', cancer_type)
        if len(detailed_output) == 1:
            continue
        num_files += 1
        cancer_name = cancer_type.replace(' ','_')
        fname_cancer_detailed = ts+'-CancerBreakdown/'+cancer_name+'-Detailed_Summary.txt'
        
        d = os.path.dirname(fname_cancer_detailed)
        if not os.path.exists(d):
            os.makedirs(d)

        fh = open(fname_cancer_detailed,'w')
        fh.write('\n'.join(detailed_output))
        fh.close()
    print 'Succesfully generated files for ',str(num_files)+'/'+str(len(cancer_type_list)), 'cancer types.'
    print

    #generate detailed summaries by Primary/Metastasis sample types
    print 'Data by sample types written to dir:',ts+'-SampleBreakdown/'
    for sample_type in ['Primary','Metastasis']:
        detailed_output = process_aa_data('sample', sample_type)

        sample_name = sample_type.replace(' ','_')
        fname_sample_detailed = ts+'-SampleBreakdown/'+sample_name+'-Detailed_Summary.txt'
        
        d = os.path.dirname(fname_sample_detailed)
        if not os.path.exists(d):
            os.makedirs(d)

        fh = open(fname_sample_detailed,'w')
        fh.write('\n'.join(detailed_output))
        fh.close()
    print 'Succesfully generated files for Primary and Metastatic samples'
    print



if __name__ == '__main__':
	main()
