#v8
#none-align to noalign,delete none-annotation
#change 'expn' to 'expr'
#add some list.txt for rmarkdown part
#recalculate the log2


import os
import sys
import pandas as pd
import subprocess
import shutil
import re
import numpy as np
import argparse
import datetime




def directory_structure(ori_samplefolder, outputpath):
    os.chdir(ori_samplefolder)
    # Change the work path to original samplefolder (the path saved all *.fastq files)
    # capture the fastq file with initial name
    global path
    path = []
    for file in os.listdir(ori_samplefolder):
        if file.endswith(".fastq"):
            samplename = os.path.splitext(file)[0]
            sampleext = os.path.splitext(file)[1]
            samplename = samplename.replace('_', '')
            new_samplefolder = os.path.join(outputpath,samplename)
            if os.path.exists(new_samplefolder):
                path.append(new_samplefolder)
                print("Work folder: \n"+new_samplefolder+ "\n is exsisting..")
            else:
                os.makedirs(new_samplefolder)
                shutil.copy(file,new_samplefolder+'/'+samplename+sampleext)
                path.append(new_samplefolder)
                print("Work folder: \n"+ new_samplefolder+'\n will be created and fastq files will be copy to..')
        elif file.endswith(".gz"):
            subprocess.call(["gunzip "+file], shell=True)
            samplename = os.path.splitext(file)[0]
            samplename_copy_use = os.path.split(file)[0]
            sampleext = os.path.splitext(file)[1]
            samplename = samplename.replace('_', '')
            new_samplefolder = os.path.join(outputpath,samplename)
            if os.path.exists(new_samplefolder):
                path.append(new_samplefolder)
                print("Work folder: \n"+new_samplefolder+ "\n is exsisting..")
            else:
                os.makedirs(new_samplefolder)
                shutil.copy(samplename_copy_use,new_samplefolder+'/'+samplename+'.fastq')
                path.append(new_samplefolder)
                print("Work folder: \n"+ new_samplefolder+'\n will be created and fastq files will be copy to..')
    print("# The work folder will be: \n" + outputpath)


def parameters(spe,outputpath):
    species = None
    refname = None
    refpath = None
    if spe == 'hsa':
        species = ' -o hsa'
        refname = 'hg38'
        refpath = '/n/data1/genomes/indexes/hg38/hg38.fa'
    elif spe == 'mmu':
        species = ' -o mmu'
        refname = 'mm9'
        refpath = '/n/data1/genomes/indexes/mm9/mm9.fa'
    elif spe == 'dme':
        species = ' -o dme'
        refname = 'dm6'
        refpath = '/n/data1/genomes/indexes/dm6/dm6.fa'
    # if ref.endswith(".fa"):
    #     fullrefname = os.path.split(ref)[1]
    #     refname = os.path.splitext(fullrefname)[0]
    #     print(refname)
    else:
        print("invalid species code, only support hsa, mmu and dme")

    if os.path.isfile(outputpath+'/'+'ref.txt'):
        os.remove(outputpath+'/'+'ref.txt')
        f = open(outputpath + '/' + 'ref.txt', 'w+')
        f.write(refpath+'\n')

    else:
        f = open(outputpath +'/'+'ref.txt','w+')
        f.write(refpath+'\n')
        f.close()
    return(species,refname,refpath)


def parameters_custom(spe,ref):
    species = None
    refname = None
    refpath = None
    if spe == 'hsa':
        species = ' -o hsa'
        refname = os.path.splitext(ref)[0]
        refpath = ref
    elif spe == 'mmu':
        species = ' -o mmu'
        refname = os.path.splitext(ref)[0]
        refpath = ref
    elif spe == 'dme':
        species = os.path.splitext(ref)[0]
        refname = 'dm6'
        refpath = ref
    # if ref.endswith(".fa"):
    #     fullrefname = os.path.split(ref)[1]
    #     refname = os.path.splitext(fullrefname)[0]
    #     print(refname)
    else:
        print("invalid species code, only support hsa, mmu and dme")
    return(species,refname,refpath)


def samplelist(ori_samplefolder,outputpath):
    if os.path.isfile(outputpath+'/'+'samplelist.txt'):
        os.remove(outputpath+'/'+'samplelist.txt')
        f = open(outputpath + '/' + 'samplelist.txt', 'w+')
        f.write(ori_samplefolder+'\n')
        for file in os.listdir(ori_samplefolder):
            if file.endswith(".fastq"):
                samplename = os.path.splitext(file)[0]
                f.write(samplename + '.fastq\n')
        f.close()
    else:
        f = open(outputpath +'/'+'samplelist.txt','w+')
        f.write(ori_samplefolder+'\n')
        for file in os.listdir(ori_samplefolder):
            if file.endswith(".fastq"):
                samplename = os.path.splitext(file)[0]
                f.write(samplename + '.fastq\n')
        f.close()


def alignments(referencebase_bwa):
    for i in range(len(path)):
        os.chdir(path[i])
        for file in os.listdir(os.getcwd()):
            if file.endswith(".fastq"):
                samplename, extension = os.path.splitext(file)
                print(
                "# The alignment will be running as :\n /n/apps/CentOS7/bin/bwa aln -t 3 " + referencebase_bwa + " " + samplename + ".fastq >" + samplename + ".sai")
                subprocess.call(["/n/apps/CentOS7/bin/bwa aln -t 3 " + referencebase_bwa + " "+samplename + ".fastq >" + samplename + ".sai"+'\n'],
                                shell=True)
                print('# building the .sam file :\n'
                "/n/apps/CentOS7/bin/bwa samse -n 10 " + referencebase_bwa + " "+samplename + ".sai " + samplename + ".fastq >" + samplename + ".sam"+'\n')
                subprocess.call([
                                    "/n/apps/CentOS7/bin/bwa samse -n 10 " + referencebase_bwa + " " + samplename + ".sai " + samplename + ".fastq >" + samplename + ".sam"],
                                shell=True)


def annotation(spec,referencebase,outputpath):
    # #path of tools
    runtools_annotation = 'perl /n/ngs/tools/tcga/v0.2.7/code/annotation/annotate.pl'
    runtools_annotation_alignment_stats = 'perl /n/ngs/tools/tcga/v0.2.7/code/library_stats/alignment_stats.pl'
    runtools_annotation_graph = 'perl /n/ngs/tools/tcga/v0.2.7/code/library_stats/graph_libs.pl'
    runtools_annotation_tcga = 'perl /n/ngs/tools/tcga/v0.2.7/code/custom_output/tcga/tcga.pl'
    runtools_expr = 'perl /n/ngs/tools/tcga/v0.2.7/code/library_stats/expression_matrix.pl'
    mirnabase = ' -m mirna_21a'

    fullcommand_annotation = runtools_annotation + mirnabase + " -u "+referencebase + spec + " -p " + outputpath
    fullcommand_align_stats = runtools_annotation_alignment_stats + " -p " + outputpath
    fullcommand_graph = runtools_annotation_graph + " -p " + outputpath
    fullcommand_tcga = runtools_annotation_tcga + mirnabase + spec + " -g "+ referencebase + " -p " + outputpath
    fullcommand_expr = runtools_expr + mirnabase +spec+ " -p " + outputpath


    print("# The annotation script will run as: ")
    print(fullcommand_annotation+'\n')
    subprocess.call([fullcommand_annotation], shell=True)

    print("# The statistics of annotation script will run as:")
    print(fullcommand_align_stats+'\n')
    subprocess.call([fullcommand_align_stats], shell=True)

    print("# The visualization of stats script will run as:")
    print(fullcommand_graph+'\n')

    print("# The tcga results output script will run as:")
    print(fullcommand_tcga+'\n')
    subprocess.call([fullcommand_tcga], shell=True)

    print("# The expression matrix script will run as: ")
    print(fullcommand_tcga+'\n')
    subprocess.call([fullcommand_expr], shell=True)


def rename(outputpath):
    df = pd.read_csv(outputpath + '/' + 'alignment_stats.csv')

    df.columns = [u'Library', u'Index', u'Total Reads greater 15bp', u'Adapter dimers',
                  u'Adapter dimers', u'Adapter at 1-14bp', u'Adapter at 15-25bp',
                  u'Adapter at 26-35bp', u'Adapter after 35bp',
                  u'Aligned Reads Post-Filter', u'percent Aligned Reads', u'Unaligned Reads',
                  u'percent Unaligned Reads', u'Filtered Reads without XA',
                  u'Softclipped Reads', u'Chastity Failed Reads Post-Filter',
                  u'detected miRNA', u'detected miRNA Covered by 10 Reads',
                  u'Total miRNA reads', u'Crossmapped miRNA reads', u'mature miRNA reads', u'star miRNA reads',
                  u'precursor miRNA reads', u'miRNA loop reads', u'unannotated miRNA reads', u'snoRNA reads',
                  u'tRNA reads', u'rRNA reads', u'snRNA reads', u'scRNA reads', u'srpRNA reads',
                  u'Other RepeatMasker RNAs reads', u'RNA (No CDS) reads', u'3UTR reads', u'5UTR reads',
                  u'Coding Exon reads', u'Intron reads', u'LINE reads', u'SINE reads', u'LTR reads', u'Satellite reads',
                  u'RepeatMasker DNA reads', u'RepeatMasker Low complexity reads',
                  u'RepeatMasker Simple repeat reads', u'RepeatMasker Other reads',
                  u'RepeatMasker Unknown reads', u'Unknown reads', u'Total miRNA',
                  u'Crossmapped miRNA', u'mature miRNA', u'star miRNA',
                  u'precursor miRNA', u'miRNA loop', u'unannotated miRNA',
                  u'snoRNA', u'tRNA', u'rRNA', u'snRNA', u'scRNA', u'srpRNA',
                  u'Other RepeatMasker RNAs', u'RNA(No CDS)', u'3UTR',
                  u'5UTR', u'Coding Exon', u'Intron', u'LINE', u'SINE',
                  u'LTR', u'Satellite', u'RepeatMasker DNA',
                  u'RepeatMasker Low complexity', u'RepeatMasker Simple repeat',
                  u'RepeatMasker Other', u'RepeatMasker Unknown', u'Unknown']
    df1 = df.iloc[:, 0:47]
    df2 = (df.iloc[:, 47:76]).replace('%', '', regex=True).astype('float') / 100
    df3 = pd.concat([df1, df2], axis=1)
    df3.to_csv(outputpath + '/' + "alignment_stats.csv", index=False, header=True)

    if os.path.isfile(outputpath + '/expn_matrix.txt') == True:
        os.rename(outputpath + '/expn_matrix.txt', outputpath + '/miRNA_expr_raw_matrix.txt')

    if os.path.isfile(outputpath + '/expn_matrix_norm.txt') == True:
        rpm = pd.read_table(outputpath + '/expn_matrix_norm.txt', index_col=0, header=0)
        rpm = rpm.round(6)
        rpm1 = rpm + 1
        rpm = rpm.astype(str)
        rpm = rpm.replace('0.000000', '0')


        log = np.log2(rpm1)
        log = pd.DataFrame(log)
        log = log.round(6)
        log = log.astype(str)
        log = log.replace('0.000000', '0')
        log.to_csv(outputpath + '/' + 'miRNA_expr_RPM_plus_1_log2_matrix.txt', sep='\t')

        rpm.to_csv(outputpath + '/' + 'miRNA_expr_RPM_matrix.txt', sep='\t')
        os.remove(outputpath + '/expn_matrix_norm.txt')
    if os.path.isfile(outputpath + '/expn_matrix_norm_log.txt') == True:
        os.remove(outputpath + '/expn_matrix_norm_log.txt')


def table_merge(outputpath):

    result_table_list_1 = ['3_UTR.txt', '5_UTR.txt', 'Intron.txt', 'LINE.txt', 'LTR.txt', 'rmsk_DNA.txt', 'rmsk_RNA.txt', 'rmsk_Simple_repeat.txt', 'rmsk_Unknown.txt',
                           'rRNA.txt', 'scRNA.txt', 'SINE.txt', 'snoRNA.txt', 'snRNA.txt', 'srpRNA.txt', 'tRNA.txt', 'Satellite.txt']

    merged_result = None

    #create the empty txt table file if that table not exist under sample_features
    for i in range(len(path)):
        sample_name = os.path.split(path[i])[1]
        workpath = path[i] + '/' + sample_name + '_features'
        os.chdir(workpath)
        for filename in result_table_list_1:
            if os.path.exists(workpath + '/' +filename):
                continue
            else:
                f = open(filename,'w+')
                f.close()

    for filename in result_table_list_1:
        result = None
        for i in range(len(path)):
            sample_name = os.path.split(path[i])[1]
            file_content = pd.read_csv(path[i] + '/' + sample_name + '_features' + '/' + filename, names=[sample_name], index_col=0)
            if result is None:
                result = file_content
            else:
                result = pd.concat([result, file_content], axis = 1)
        result['type'] = filename.split('.')[0]
        if merged_result is None:
            merged_result = result
        else:
            merged_result = pd.concat([merged_result, result], axis=0)

    # set the missing value
    merged_result = merged_result.fillna('NA')

    #switch the type column with the second column
    type = merged_result.pop('type')
    merged_result.insert(0,'type',type)

    merged_result.to_csv(outputpath+'/'+"Feature_summary_table.txt", sep='\t')
    return merged_result


def filter_miRNA_table():

    for i in range(len(path)):
        sample_name = os.path.split(path[i])[1]
        workpath = path[i] + '/' + sample_name + '_features'
        os.chdir(workpath)

        name = []
        category = []
        mir_id = []
        value = []
        content = []

        file = open(workpath+"/"+ "miRNA.txt")
        reader = file.readlines()
        for line in reader:
            if (re.split('[,\s]',line)[1] == 'unannotated') or (re.split('[,\s]',line)[1]) == "precursor" or (re.split('[,\s]',line)[1] == "stemloop"):
                name.append(re.split('[,\s]', line)[0])
                category.append(re.split('[,\s]', line)[1])
                mir_id.append(re.split('[,\s]', line)[1])
                value.append(re.split('[,\s]', line)[2])
            else:
                name.append(re.split('[,\s]', line)[0])
                category.append(re.split('[,\s]', line)[1])
                mir_id.append(re.split('[,\s]', line)[2])
                value.append(re.split('[,\s]', line)[3])
        for i in range(len(value)):
            content.append(name[i] + " " + category[i] + " " + mir_id[i] + "," + value[i])

        output = open(workpath +"/" + "miRNA_filtered.txt", "w")
        for i in range(len(content)):
            output.write(content[i])
            output.write('\n')


# def table_merge_miRNA():
#
#     result_table_list_1 = ['miRNA_filtered.txt']
#
#     merged_result = None
#     print(path)
#     #create the empty txt table file if that table not exist under sample_features
#     for i in range(len(path)):
#         sample_name = os.path.split(path[i])[1]
#         print(sample_name)
#         workpath = path[i] + '/' + sample_name + '_features'
#         os.chdir(workpath)
#         for filename in result_table_list_1:
#             if os.path.exists(workpath + '/' +filename):
#                 continue
#             else:
#                 f = open(filename,'w+')
#                 f.close()
#
#     for filename in result_table_list_1:
#         result = None
#         for i in range(len(path)):
#             sample_name = os.path.split(path[i])[1]
#             file_content = pd.read_csv(path[i] + '/' + sample_name + '_features' + '/' + filename, names=[sample_name], index_col=0)
#             print(file_content)



    #         if result is None:
    #             result = file_content
    #         else:
    #             result = pd.concat([result, file_content], axis = 1)
    #     result['type'] = filename.split('.')[0]
    #     if merged_result is None:
    #         merged_result = result
    #     else:
    #         merged_result = pd.concat([merged_result, result], axis=0)
    #
    # # set the missing value
    # merged_result = merged_result.fillna('NA')
    #
    # #switch the type column with the second column
    # type = merged_result.pop('type')
    # merged_result.insert(0,'type',type)
    # return merged_result
    #


def table_merge_miRNA(outputpath):

    result_table_list_1 = ['miRNA_filtered.txt']

    result = None

    #create the empty txt table file if that table not exist under sample_features
    for i in range(len(path)):
        sample_name = os.path.split(path[i])[1]
        workpath = path[i] + '/' + sample_name + '_features'
        os.chdir(workpath)
        for filename in result_table_list_1:
            if os.path.exists(workpath + '/' +filename):
                continue
            else:
                f = open(filename,'w+')
                f.close()

    table_head = ['name', 'type', 'mir_id']

    for filename in result_table_list_1:
        result = None
        for i in range(len(path)):

            sample_name = os.path.split(path[i])[1]
            table_head.append(sample_name)


            file_content = pd.read_csv(path[i] + '/' + sample_name + '_features' + '/' + filename, header= None, index_col=0)


            if result is None:

                result = file_content
            else:
                result = pd.concat([result, file_content], axis = 1)

    result = result.fillna(0)       #########can not set the fill content as NA, since it leads misplacement in next step, to FIX!!!!!##########
    # set the missing value
    result.to_csv(outputpath+'/'+"miRNA_expr_summary_table.txt",sep='\t',header=False)
    merged_result = pd.read_csv(outputpath+'/'+"miRNA_expr_summary_table.txt", sep='\s+',names=table_head)
    merged_result.to_csv(outputpath+'/'+"miRNA_expr_summary_table.txt",sep='\t',index=False)
    return result
    #
    #


def merge_mirna_taglen(outputpath):
    result = None

    for i in range(len(path)):

        sample_name = os.path.split(path[i])[1]
        # print(sample_name)
        workpath = path[i] + '/' + sample_name + '_features'
        os.chdir(workpath)
        file_content = pd.read_csv('filtered_taglengths.csv',header=0)
        file_content.insert(0,'sample',sample_name)
        file_content.to_csv('filtered_taglengths.txt', sep='\t',index=False, header=True)

        if result is None:
            result = file_content
        else:
            result = pd.concat([result, file_content], axis = 0)
    df = pd.DataFrame(result)
    col_list = list(df)
    col_list.remove('sample')
    col_list.remove('taglen')
    df['total'] = df[col_list].sum(axis=1)
    df.to_csv(outputpath+'/'+"merged_taglen_table.csv", index=False,header=True)


def remove_files():
    print("# The redundant file will be transformed and removed as : ")
    for i in range(len(path)):
        os.chdir(path[i])
        for file in os.listdir(os.getcwd()):
            if file.endswith(".fastq"):
                os.remove(file)
            if file.endswith(".sam"):
                samplename = os.path.splitext(file)[0]
                print("\n samtools view -Su "+ file + ">" + samplename +".bam ")
                subprocess.call(["samtools view -Su "+ file + ">" + samplename +".bam "],shell=True)

                print("\n samtools sort " + samplename +".bam -o " + samplename + ".sorted.bam ")
                subprocess.call(["samtools sort " + samplename +".bam -o " + samplename + ".sorted.bam "], shell=True)

                print("\n samtools index "+samplename+".sorted.bam")
                subprocess.call(["samtools index "+samplename+".sorted.bam"], shell=True)
                os.remove(file)
                os.remove(samplename+'.bam')
                os.remove(samplename+'.sai')


def main():
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-i", "--input", action="store",dest='fastq_dir',required=True, help="Please enter the path of sample .fastq folder")
    requiredNamed.add_argument("-o", "--ouput", action="store",dest='out_dir',required=True, help="Please enter the path of directory to save result")
    requiredNamed.add_argument("-s", "--species", action="store",dest='specie',required=True, help="Please enter the species type of your sample: (hsa, mmu or dme)")

    parser.add_argument("--noalign",action="store_true",default=False, dest='alignment', help="if you already alignment .sam files under your sample work folder, and you may skip alignment step please enter this argument")
    parser.add_argument("-r", "--reference", action="store", default=None, dest='ref', help="indicate the custom reference instead of default e.g. /n/data1/genomes/indexes/hg38/hg38.fa")
    args = parser.parse_args()

    f = open(args.fastq_dir+ '/' + "logging.txt", 'w+')
    logging_path= (args.fastq_dir+ '/' + "logging.txt")
    sys.stdout = f
    print('# The used arguments are:')
    a = 'python '
    for i in range(len(sys.argv)):
        a = a + ' ' + sys.argv[i]

    b = a.replace("'", '')
    b = b.replace("[", '')
    b = b.replace(",", '')
    b = b.replace("]", '')

    print(b)

    directory_structure(args.fastq_dir, args.out_dir)

    # parameter settings
    if args.ref is None:
        tcga_arg = parameters(args.specie,args.out_dir)
        print('# The using reference genome index is:')
        print(tcga_arg[2] + '\n')
    else:
        tcga_arg = parameters_custom(args.specie,args.ref)
        print('# The using reference genome index is:')
        print(tcga_arg[2] + '\n')

    samplelist(args.fastq_dir, args.out_dir)


    if args.alignment == False:
        alignments(tcga_arg[2])

    annotation(tcga_arg[0], tcga_arg[1], args.out_dir)


    rename(args.out_dir)

    table_merge(args.out_dir)

    filter_miRNA_table()

    table_merge_miRNA(args.out_dir)

    merge_mirna_taglen(args.out_dir)

    shutil.copy("/n/core/Bioinformatics/research/cxu/mirna_pipeline/python_pipeline/rmarkdown_mod.Rmd", args.out_dir)
    remove_files()

    today = datetime.date.today()
    print('\n')
    print('Date:')
    print(today)


    f.close()
    if os.path.isfile(args.out_dir + '/logging.txt') == True:
        os.remove(args.out_dir + '/logging.txt')
    shutil.move(logging_path,args.out_dir)
if __name__ == '__main__':
    main()
