import os
import glob

def extract_data_from_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        copyn_data = {}
        for idx, line in enumerate(lines):
            if "CopyN Case/Control" in line:
                i = idx + 1
                while i < len(lines) and lines[i].strip():
                    parts = lines[i].strip().split()
                    #print (parts[2])
                    copyn = parts[0]
                    if copyn == "1":
                     copyn = "del"
                    if copyn == "3":
                     copyn = "dup"
                    case = parts[1]
                    control = parts[3]
                    copyn_data[copyn] = (case, control)
                    i += 1
                break
    return copyn_data

def main():
    # 根据您给出的路径模式，查找所有匹配的文件
    file_pattern = '/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/intermediate/gene/combined_ENSG*_allCNV_Allfreq_20kb/*.log'
    files = glob.glob(file_pattern)
    
    output_file = '/stanley/huang_lab/home/ychen/proj-CNV/07CNV_gene_association/01EAS_gene_assoc/output/EAS_freq_output.csv'
    with open(output_file, 'w') as out:
        out.write('File Name,CopyN,Case,Control\n')  # 写入表头
        for file_path in files:
            file_name = os.path.basename(file_path)
            data = extract_data_from_file(file_path)
            for copyn, (case, control) in data.items():
                out.write(f'{file_name},{copyn},{case},{control}\n')

if __name__ == '__main__':
    main()
