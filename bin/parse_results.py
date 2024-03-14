#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import logging
import sys
from pathlib import Path
import re
from collections import OrderedDict

logger = logging.getLogger()

def parse_busco(out_short_summary):
    dict_busco = {}
    with open(out_short_summary, 'r') as infbusco:
        for idxline, line in enumerate(infbusco.readlines()):
            # line = line.strip("\n")
            if line.startswith("#"):
                match_db_obj = re.match(r'.+\(.+number of genomes\:\s(\d+),\snumber of BUSCOs\:\s(\d+)\)', line)
                #if match_db_obj:
                #    print(match_db_obj.group(1), match_db_obj.group(2))
            else:
                match_summary_obj = re.match(r'\s+C\:(\d+.\d+)\%\[S\:(\d+.\d)\%,D\:(\d+.\d+)\%\],F\:(\d+.\d+)\%,M:(\d+.\d+)%,n:(\d+)', line)
                if match_summary_obj:
                    dict_busco['pctcomplete'] = match_summary_obj.group(1)
                    dict_busco['pctsingle'] = match_summary_obj.group(2)
                    dict_busco['pctduplicated'] = match_summary_obj.group(3)
                    dict_busco['pctfragmented'] = match_summary_obj.group(4)
                    dict_busco['pctmissing'] = match_summary_obj.group(5)
                    dict_busco['total_busco_searched_genes'] = match_summary_obj.group(6)
                #    print(idxline)
                elif idxline == 9:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['ncomplete'] = n
                elif idxline == 10:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nsingle'] = n
                elif idxline == 11:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nduplicated'] = n
                elif idxline == 12:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nfragmented'] = n
                elif idxline == 13:
                    match_n = re.match(r'\s+(\d+)\s+.+', line)
                    n = match_n.group(1)
                    dict_busco['nmissing'] = n
    return(dict_busco)

def parse_compleasm(out_short_summary):
    dict_busco = {}
    with open(out_short_summary, 'r') as infbusco:
        for idxline, line in enumerate(infbusco.readlines()):
            # line = line.strip("\n")
            if line.startswith("#"):
                match_db_obj = re.match(r'## lineage:\s(\S+)$', line)
                    # S + D
                    
                    # S (Single Copy Complete Genes): The BUSCO genes that can be entirely aligned in the assembly, with only one copy present.
                    
                    # D (Duplicated Complete Genes): The BUSCO genes that can be completely aligned in the assembly, with more than one copy present.
                    
                    # F (Fragmented Genes, subclass 1): The BUSCO genes which only a portion of the gene is present in the assembly, and the rest of the gene cannot be aligned.
                    # PLUS
                    # I (Fragmented Genes, subclass 2): The BUSCO genes in which a section of the gene aligns to one position in the assembly, while the remaining part aligns to another position.
                    
                    # M (Missing Genes): The BUSCO genes with no alignment present in the assembly.
                    
            elif line.startswith("S"):
                match_n = re.match(r'S:(\S+)%,\s(\d+)$', line)
                dict_busco['pctsingle'] = match_n.group(1)
                dict_busco['nsingle'] = match_n.group(2)
            elif line.startswith("D"):
                match_n = re.match(r'D:(\S+)%,\s(\d+)$', line)
                dict_busco['pctduplicated'] = match_n.group(1)
                dict_busco['nduplicated'] = match_n.group(2)
            elif line.startswith("F"):
                match_n = re.match(r'F:(\S+)%,\s(\d+)$', line)
                dict_busco['pctfragmented'] = match_n.group(1)
                dict_busco['nfragmented'] = match_n.group(2)
            elif line.startswith("I"):
                match_n = re.match(r'I:(\S+)%,\s(\d+)$', line)
                dict_busco['pctincomplete'] = match_n.group(1)
                dict_busco['nincomplete'] = match_n.group(2)
            elif line.startswith("M"):
                match_n = re.match(r'M:(\S+)%,\s(\d+)$', line)
                dict_busco['pctmissing'] = match_n.group(1)
                dict_busco['nmissing'] = match_n.group(2)
            else:
                continue
        
    # COMPLETE = S + D
    dict_busco['pctcomplete'] = dict_busco['pctsingle'] + dict_busco['pctduplicated']
    dict_busco['ncomplete'] = dict_busco['nsingle'] + dict_busco['nduplicated']

    # FRAGMENTED GENES IN BUSCO = 
    # F (Fragmented Genes, subclass 1): The BUSCO genes which only a portion of the gene is present in the assembly, and the rest of the gene cannot be aligned.
    # PLUS
    # I (Fragmented Genes, subclass 2): The BUSCO genes in which a section of the gene aligns to one position in the assembly, while the remaining part aligns to another position.
    # dict_busco['pctfragmented'] = dict_busco['pctfragmented'] + dict_busco['pctincomplete']
    # dict_busco['nfragmented'] = dict_busco['nfragmented'] + dict_busco['nincomplete']

    return(dict_busco)

def parse_results_to_table(genomes_ids, ale_res, reapr_res, busco_re_summary, quast_res, file_out, busco, merfin_qv_res, merfin_comp_res):
                # [ ...
                # [ sample_id, [ ale_res ], [reapr_res], [busco_re_summary], [quast_res]],
                # [ sample_id_B, [ ale_res_B ], [reapr_res_B], [busco_re_summary_B], [quast_res_B]]
                #  ... ]
                # ale_res=expand('evaluate_assembly/{sample}/ALEScore_{sample}.finished',sample=samples['sample']),
                # reapr_res=expand('evaluate_assembly/{sample}/REAPR_{sample}.finished',sample=samples['sample']),
                # busco_res=expand('evaluate_assembly/{sample}/busco/BUSCO_{sample}.finished', sample=samples['sample']),
                # quast_res='evaluate_assembly/quast_results/QUAST.OK'
        output = "evaluate_assembly/results.html"
        items = []
        items_classes = []
        my_list = [item.strip() for item in genomes_ids[1:-1].split(',')]
        my_list_ale = [item.strip() for item in ale_res[1:-1].split(',')]
        my_list_reapr = [item.strip() for item in reapr_res[1:-1].split(',')]
        my_list_busco = [item.strip() for item in busco_re_summary[1:-1].split(',')]
        my_list_quast = [item.strip() for item in quast_res[1:-1].split(',')]
        my_merfin_qv_file = [item.strip() for item in merfin_qv_res[1:-1].split(',')]
        my_merfin_completeness_file = [item.strip() for item in merfin_comp_res[1:-1].split(',')]
        for idx, pseudo_samp in enumerate(my_list):
                # samp = pseudo_samp.split("/")[-1]
                samp = pseudo_samp
                dict_sample = OrderedDict()
                dict_classes = OrderedDict()
                dict_sample['name'] = samp
#                        print(pseudo_samp, samp)
                # busco_summary_file = "evaluate_assembly/{}/run_{}/short_summary_{}.txt".format(pseudo_samp, samp, samp)
                if busco:
                    busco_summary_file = my_list_busco[idx]
                    if Path(busco_summary_file).exists():
                            dict_samp = OrderedDict(parse_busco(busco_summary_file))
                    else:
                            dict_samp['ncomplete'] = 'X'
                            dict_samp['pctcomplete'] = 'X'
                            dict_samp['nduplicated'] = 'X'
                            dict_samp['pctduplicated'] = 'X'
                            dict_samp['nfragmented'] = 'X'
                            dict_samp['pctfragmented'] = 'X'
                else:
                    compleasm_summary_file = my_list_busco[idx]
                    if Path(compleasm_summary_file).exists():
                            dict_samp = OrderedDict(parse_compleasm(compleasm_summary_file))
                # genome_path = samples.loc[samp, "assembly"]
                # genome_name  = genome_path.split('/')[-1]
                # ale_file = "evaluate_assembly/{}/ALEoutput.txt".format(pseudo_samp, genome_name)
                ale_file = my_list_ale[idx]
                if Path(ale_file).exists():
                    with open(ale_file) as infreapr:
                            for line in infreapr.readlines():
                                    match_score = re.match(r'#\sALE_score:\s(-\d+.\d+)', line)
                                    if match_score:
                                            ale_score = float(match_score.group(1))
                dict_sample['ale'] = float(ale_score)


                # MATCH MERFIN RESULTS - COMPLETENESS
                merfin_file_completeness = my_merfin_completeness_file[idx]
                if Path(merfin_file_completeness).exists():
                    with open(merfin_file_completeness) as infmc:
                            for line in infmc.readlines():
                                    match_comp = re.match(r'^COMPLETENESS:\s+(\S+)$', line)
                                    if match_comp:
                                            comp_score = float(match_comp.group(1))
                dict_sample['merfin_completeness'] = float(comp_score)

                # MATCH MERFIN RESULTS - QV*
                merfin_file_qv = my_merfin_qv_file[idx]
                if Path(merfin_file_qv).exists():
                    with open(merfin_file_qv) as infmc:
                            for line in infmc.readlines():
                                    match_comp = re.match(r'^Merfin\sQV\*:\s(\S+)$', line)
                                    if match_comp:
                                            qv_score = float(match_comp.group(1))
                dict_sample['merfin_qv_ast'] = float(qv_score)

#                        ale_file = "evaluate_assembly/{}/ALEoutput.txt".format(pseudo_samp, genome_name)
#                        if path.exists(ale_file):
#                                with open(ale_file, 'r') as inale:
#                                        nline = 0
#                                        for line in inale.readlines():
#                                                if nline == 0:
#                                                        line_list = line.split(" ")
#                                                        dict_sample['ale'] = float(line_list[-1].strip('\n'))
                # BEGIN - PARSING REAPR RESULTS
                reapr_file = my_list_reapr[idx]
                if Path(reapr_file).exists():
                    with open(reapr_file) as infreapr:
                            for line in infreapr.readlines():
                                    match_errors = re.match(r'^(\d+)\serrors.$', line)
                                    if match_errors:
                                            reapr_errors = match_errors.group(1)
                                    match_fcd_errors = re.match(r'FCD errors within a contig:\s(\d+)', line)
                                    if match_fcd_errors:
                                            fcd_errors = match_fcd_errors.group(1)
                                    match_low_frag_cov = re.match(r'Low fragment coverage within a contig:\s(\d+)', line)
                                    if match_low_frag_cov:
                                            low_frag_errors = match_low_frag_cov.group(1)
                    dict_sample['reapr_total_errors'] = reapr_errors
                    dict_sample['reapr_fcd'] = fcd_errors
                    dict_sample['reapr_low'] = low_frag_errors
                # END - PARSING REAPR RESULTS
#                        reapr_file = "evaluate_assembly/{}/reapr_results/05.summary.report.txt".format(pseudo_samp)
#                        if path.exists(reapr_file):
#                                with open(reapr_file, 'r') as inreapr:
#                                        for line in inreapr.readlines():
#                                                if line.endswith('errors:\n'):
#                                                        line_list = line.split(" ")
                #                                 dict_sample['reapr'] = line_list[0]
                # parsing quast results
                # df_quast = pd.read_table('evaluate_assembly/quast_results/report.tsv', sep='\t')
                df_quast = pd.read_table(my_list_quast[idx], sep='\t')
                array_sample = df_quast.iloc[:,1].values
                print(array_sample)
                dict_sample['genomesize'] = '{:,}'.format(int(array_sample[6]))
                dict_sample['contigs']    = '{:,}'.format(int(array_sample[12]))
                dict_sample['n50']        = '{:,}'.format(int(array_sample[16]))
                dict_sample['largest']    = '{:,}'.format(int(array_sample[13]))

                dict_sample['genomesize'] = int(array_sample[6])
                dict_sample['contigs']    = int(array_sample[12])
                dict_sample['n50']        = int(array_sample[16])
                dict_sample['largest']    = int(array_sample[13])
                # concatenating results
                # dict_sample['genomesize_class'] = "tg-lboi"
                dict_appended = {**dict_sample, **dict_samp}
                for k, v in dict_appended.items():
                        key_class = '{}_class'.format(k)
                        dict_classes[key_class] = "tg-lboi"
                items.append(dict_appended)
                items_classes.append(dict_classes)
        print(items)
        ale_scores = []
        for it in items:
            ale_scores.append(it['ale'])
        ale_scores = np.array(ale_scores)
        for it in items:
            # min-max normalization
            
            it['ale_norm'] = '{:.2f}'.format(((it['ale']-np.nanmin(ale_scores))/(np.nanmax(ale_scores)-np.nanmin(ale_scores))))
            # print(it['ale'])
        #myList = [list(col) for col in zip(*[d.values() for d in items])]
        #myList_argmax = np.argmax(myList, axis=1)
        #print("argmax",myList_argmax)

        #for idx, arg in enumerate(myList_argmax):
#               key_checked = list(dict_appended.keys())[idx]
#                       key_class = "{}_class".format(key_checked)
#                       items_classes[arg][key_class] = "mark"
                #print("K",items[arg][key_checked])
        # res_items = []
        # df = pd.DataFrame(items)
        # df_unnameA = df.drop('name', axis=1)
        # df_unnameA = df_unnameA.drop('ale_norm', axis=1)
        # print(df_unnameA)
        # df_unname = df_unnameA.apply(pd.to_numeric)
        # print("HERE")
        # df_unname['pct_nonduplicated'] = 100. - df_unname['pctduplicated']
        # df_unname['pct_integral'] = 100. - df_unname['pctfragmented']
        # df_unname['pct_found'] = 100. - df_unname['pctmissing']
        # df_sub = df_unname[['genomesize', 'contigs', 'n50', 'largest', 'pctcomplete', 'pct_nonduplicated', 'pct_integral', 'pct_found']]
        # df_sub['name'] = df['name']
        # df_sub_sorted = df_sub.sort_values('name')
        # df_sub_sorted.columns = ['Genome Size (bp)', 'Number of Contigs', 'N50', 'Largest Contig (bp)', 'BUSCO Complete Genes (%)', 'BUSCO Single-Copy Genes (%)', 'BUSCO Non-fragmented Genes (%)', 'BUSCO Found Genes (%)', 'Assembly']
        # df_sub_sorted = df_sub_sorted[['Assembly', 'Genome Size (bp)', 'Number of Contigs', 'N50', 'Largest Contig (bp)', 'BUSCO Complete Genes (%)', 'BUSCO Single-Copy Genes (%)', 'BUSCO Non-fragmented Genes (%)', 'BUSCO Found Genes (%)']]
        # k = df_sub_sorted.style.hide_index().background_gradient('viridis', axis=0, subset=['Genome Size (bp)', 'Number of Contigs', 'N50', 'Largest Contig (bp)', 'BUSCO Complete Genes (%)', 'BUSCO Single-Copy Genes (%)', 'BUSCO Non-fragmented Genes (%)', 'BUSCO Found Genes (%)'])
        # with open('assets/results_heat.html', 'w') as fheat:
        # 	fheat.write(k.render())

# 		for t in list(zip(items, items_classes)):
# 				nd = {**t[0], **t[1]}
# 				res_items.append(nd)
# 		loader = jinja2.FileSystemLoader('template.html')
# 		env = jinja2.Environment(loader=loader)
# 		template = env.get_template('')
# 		output_jinja2 = template.render(items=res_items)
# #               print(output_jinja2)
# 		with open(output[0], 'w') as outfile:
# 				outfile.write(output_jinja2)
        # c = dict([(k,[a[k],b[k]]) for k in items])
        c = pd.DataFrame(items)
        if busco:
            # {'name': 'Sample_AssemblerA', 'ale': -100937528.909207, 'reapr_total_errors': '402', 'reapr_fcd': '196', 'reapr_low': '206', 'genomesize': 4090859, 
            # 'contigs': 2, 'n50': 4065161, 'largest': 4065161, 'pctcomplete': '6.2', 'pctsingle': '5.7', 'pctduplicated': '0.5', 'pctfragmented': '3.0', 'pctmissing': '90.8', 'total_busco_searched_genes': '758', 'ncomplete': '47', 'nsingle': '43', 'nduplicated': '4', 'nfragmented': '23', 'nmissing': '688'}, {'name': 'Sample_AssemblerB', 'ale': -100937528.909207, 'reapr_total_errors': '402', 'reapr_fcd': '196', 'reapr_low': '206', 'genomesize': 4090859, 'contigs': 2, 'n50': 4065161, 'largest': 4065161, 'pctcomplete': '6.2', 'pctsingle': '5.7', 'pctduplicated': '0.5', 'pctfragmented': '3.0', 'pctmissing': '90.8', 'total_busco_searched_genes': '758', 'ncomplete': '47', 'nsingle': '43', 'nduplicated': '4', 'nfragmented': '23', 'nmissing': '688'}
            dict_names = {'name': 'Assembly', 'ale': 'ALE score (neglog)', 'reapr_total_errors': 'REAPR erros', 'reapr_fcd': 'REAPR fcd', 'reapr_low': 'REAPR low',
             'genomesize': 'Assembly length', 'contigs': 'contigs', 'n50': 'N50', 'largest': 'Largest contig',
             'pctcomplete': 'BUSCO complete (%)', 'pctsingle': 'BUSCO single (%)', 'pctduplicated': 'BUSCO duplicated (%)', 'pctfragmented':'BUSCO fragmented (%)',
              'pctmissing': 'BUSCO missing (%)',
               'ncomplete': 'BUSCO complete', 'nsingle': 'BUSCO single', 'nduplicated': 'BUSCO duplicated',
                'nfragmented': 'BUSCO fragmented', 'nmissing': 'BUSCO missing','total_busco_searched_genes': 'BUSCO SEARCHED GENES', 'ale_norm':'ALE normalized',
                'merfin_completeness': 'Merfin Completness','merfin_qv_ast': 'Merfin QV*'}

            c.rename(columns=dict_names, inplace=True)


            # c.columns = ['Assembly', 'ALE score (neglog)', 'REAPR erros', 'REAPR fcd', 'REAPR low',
            #  'Assembly length', 'contigs', 'N50', 'Largest contig', 'BUSCO complete (%)', 'BUSCO single (%)', 'BUSCO duplicated (%)', 'BUSCO fragmented (%)', 'BUSCO missing (%)', 'BUSCO complete', 'BUSCO single', 'BUSCO duplicated', 'BUSCO fragmented', 'BUSCO missing', 'ALE normalized']
        else:
            dict_names = {'name': 'Assembly', 'ale': 'ALE score (neglog)', 'reapr_total_errors': 'REAPR erros', 'reapr_fcd': 'REAPR fcd', 'reapr_low': 'REAPR low',
             'genomesize': 'Assembly length', 'contigs': 'contigs', 'n50': 'N50', 'largest': 'Largest contig',
             'pctcomplete': 'COMPLEASM complete (%)', 'pctsingle': 'COMPLEASM single (%)', 'pctduplicated': 'COMPLEASM duplicated (%)', 'pctfragmented':'COMPLEASM fragmented Class I (%)', 'pctincomplete':'COMPLEASM fragmented Class II (%)',
              'pctmissing': 'COMPLEASM missing (%)',
               'ncomplete': 'COMPLEASM complete', 'nsingle': 'COMPLEASM single', 'nduplicated': 'COMPLEASM duplicated', 'nfragmented': 'COMPLEASM fragmented Class I', 
               'nincomplete': 'COMPLEASM fragmented Class II', 'nmissing': 'COMPLEASM missing','total_busco_searched_genes': 'COMPLEASM SEARCHED GENES',
                'ale_norm':'ALE normalized', 'merfin_completeness': 'Merfin Completness','merfin_qv_ast': 'Merfin QV*'}
            c.rename(columns=dict_names, inplace=True)
        # c.to_excel("evaluate_assembly/results.xlsx", index=False)
        c.to_csv(file_out, index=False, sep="\t")
        print("Success ! The results summary table has been written ! \n The results can be view in:\n \t- Excel format in file evaluate_assembly/results.xlsx \n \t- HTML format in file evaluate_assembly/results.html \n \t- HTML heatmap in file evaluate_assembly/results_head.html \n \t- CSV format in file evaluate_assembly/results.csv")


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Transform record fasta header if it contains things like trailing whitespace or characters |':- that could break the REAPR execution.",
        epilog="Example: python check_header_fasta.py sample.fasta sample_valid.fasta",
    )
    parser.add_argument(
        "-g",
        "--genomes_ids",
        # nargs='+',
        metavar="GENOMES_IDS",
        # type=Path,
        help="Fasta input.",
    )
    parser.add_argument(
        "-a",
        "--ale_res",
        # nargs='+',
        metavar="ALE_RES",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "-r",
        "--reapr_res",
        # nargs='+',
        metavar="REAPR_RES",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "-b",
        "--busco_re_summary",
        # nargs='+',
        metavar="BUSCO_RE_SUMMARY",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "-q",
        "--quast_res",
        # nargs='+',
        metavar="QUAST_RES",
        # type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "--merfin_qv_res",
        # nargs='+',
        metavar="MERFIN_QV_RES",
        # type=Path,
        help="Merfin hist log outputs.",
    )

    parser.add_argument(
        "--merfin_comp_res",
        # nargs='+',
        metavar="MERFIN_COMP_RES",
        # type=Path,
        help="Merfin completness output.",
    )

    parser.add_argument(
        "-f",
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Fasta output.",
    )

    parser.add_argument(
        "--busco",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        metavar="BUSCO",
        help="Active script to parse BUSCO results instead of COMPLEASM.",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    # if not args.file_out.is_file():
    #     logger.error(f"The given input file {args.file_in} was not found!")
    #     sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    # check_fasta(args.file_in, args.file_out)
    print(args)
    parse_results_to_table(args.genomes_ids, args.ale_res, args.reapr_res, args.busco_re_summary, args.quast_res, args.file_out, args.busco, args.merfin_qv_res, args.merfin_comp_res)


if __name__ == "__main__":
    sys.exit(main())
        