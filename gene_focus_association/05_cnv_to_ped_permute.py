import sys

target_CNV_type = sys.argv[1]
WORKING_DIR = "/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc"
cnv_file = "/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/input/cfile_final.cnv"
phe_file = "/stanley/huang_lab/home/ychen/proj-CNV/misc/EAS_all_cohorts_final.phe"
genelist_file = "/stanley/huang_lab/home/ychen/proj-CNV/misc/ensemble_exon_genelist.hg19.txt"

dup_threshold_count = 8
del_threshold_count = 8
all_threshold_count = 8


if __name__ == "__main__":
    transcript_by_chr = {}
    ids_affected_gene = {}
    for i in range(1, 24):
        transcript_by_chr[str(i)] = [] # [ [BP1, BP2, transcript_id, gene_id] , [BP1, BP2, transcript_id, gene_id] ... ]
        # ids_affected_gene= {} # {chr -> {gene -> {id -> CNV_number}}}

    # read genelist file
    gene2info = {}
    line_count = 0
    for line in open(genelist_file, "r"):
        parsed_line = line.strip("\n").split("\t")
        if line_count == 0:
            col2idx = {}
            for i in range(len(parsed_line)):
                col2idx[parsed_line[i]] = i
        else:
            if parsed_line[col2idx['cdsStartStat']] != 'cmpl' or parsed_line[col2idx['cdsEndStat']] != 'cmpl':
                continue
            CHR = parsed_line[col2idx['chrom']][3: ]
            if CHR == 'X':
                CHR = '23'
            if CHR in [chr for chr, _ in transcript_by_chr.items()]:
                transcript_by_chr[CHR].append([int(parsed_line[col2idx['cdsStart']]), int(parsed_line[col2idx['cdsEnd']]), parsed_line[col2idx['name']], parsed_line[col2idx['name2']]])
                ids_affected_gene[parsed_line[col2idx['name2']]] = {}
                gene2info[parsed_line[col2idx['name2']]] = [CHR, parsed_line[col2idx['name2']], "0", parsed_line[col2idx['cdsStart']], "1", "2"]
        line_count += 1
    
    # # Read .fam file
    barcodes = []
    barcode2FID = {}
    line_count = 0
    for line in open(phe_file, "rt"):
        parsed_line = line.strip("\n").split()
        if line_count == 0:
            col2idx = {}
            for i in range(len(parsed_line)):
                col2idx[parsed_line[i]] = i
        else:
            barcodes.append(parsed_line[col2idx['IID']])
            barcode2FID[parsed_line[col2idx['IID']]] = parsed_line[col2idx['FID']]
        line_count += 1
    
    print(str(len(barcodes)) + " individuals in phenotype file.")
    
    gene_cnvIDs = {}
    # Read .cnv file
    line_count = 0
    for line in open(cnv_file, "r"):
        parsed_line = line.strip("\n").split()
        if line_count == 0:
            col2idx = {}
            for i in range(len(parsed_line)):
                col2idx[parsed_line[i]] = i
        else:
            barcode = parsed_line[col2idx['IID']]
            BP1 = int(parsed_line[col2idx['BP1']])
            BP2 = int(parsed_line[col2idx['BP2']])
            CHR = parsed_line[col2idx['CHR']]
            CNV_number = int(parsed_line[col2idx['TYPE']])
            if (CNV_number > 2 and target_CNV_type == 'dup') or (CNV_number < 2 and target_CNV_type == 'del') or target_CNV_type == "all":
                for [gene_BP1, gene_BP2, transcript, gene] in transcript_by_chr[CHR]:
                    if (BP1 > gene_BP1 and BP1 < gene_BP2) or (BP2 > gene_BP1 and BP2 < gene_BP2) or (BP1 < gene_BP1 and BP2 > gene_BP2):
                        ids_affected_gene[gene][barcode] = CHR + "-" + str(BP1) + "-" + str(BP2)
        line_count += 1
    
    assoc_genes = []
    carriers2genes = {}
    for gene, barcode2CNV in ids_affected_gene.items():
        if (len(barcode2CNV.items()) >= dup_threshold_count and target_CNV_type == "dup") or (len(barcode2CNV.items()) >= del_threshold_count and target_CNV_type == "del") or (len(barcode2CNV.items()) >= all_threshold_count and target_CNV_type == "all"):
            assoc_genes.append(gene)
            carriers_str = ";".join(sorted([barcode + "-" + CNV for barcode, CNV in barcode2CNV.items()]))
            if not carriers2genes.get(carriers_str):
                carriers2genes[carriers_str] = []
            carriers2genes[carriers_str].append(gene)
            
    
    print("{} genes kept.".format(len(assoc_genes)))
    print("{} independent tests.".format(len(carriers2genes.keys())))

    dupIDs = [
        'Ma_xajd*205403310065_R06C02',
        'Release_Qin_biox6*200647740031_R01C01',
        'Release_Qin_biox6*200647740031_R02C01',
        'Release_Qin_biox6*200647740031_R03C01',
        'Release_Qin_biox6*200647740031_R04C01',
        'Release_Qin_biox6*200647740031_R05C01',
        'Release_Qin_biox6*200647740031_R06C01',
        'Release_Qin_biox6*200647740031_R07C01'
    ]
    
    output = open(WORKING_DIR + "/input/cnv_bed/cnv_permute." + target_CNV_type + ".ped", "w+")
    for barcode in barcodes:
        if barcode in dupIDs or "Tai1_multiplex" in barcode or "tai2trios" in barcode or "Staging_Tsuang_tai3" in barcode:
            continue
        row = [barcode2FID[barcode], barcode, "0", "0", "0", "0"]
        for gene in assoc_genes:
            if ids_affected_gene[gene].get(barcode):
                row += ["1", "2"]    # has CNV
            else:
                row += ["1", "1"]    # no CNV
        output.write("\t".join(map(str, row)) + "\n")      
        
    output = open(WORKING_DIR + "/input/cnv_bed/cnv_permute." + target_CNV_type + ".map", "w+")  
    for gene in assoc_genes:    
        output.write(" ".join(gene2info[gene]) + "\n")  
