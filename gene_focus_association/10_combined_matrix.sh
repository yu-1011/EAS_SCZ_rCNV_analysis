#!/bin/bash
#$ -cwd
#$ -pe smp 4 -binding linear:4
#$ -l h_rt=12:00:00,h_vmem=15G,
#$ -N gene_assoc_gene_cnv_matrix

directory=$1  # e.g., /path/to/your/files
output_file="merged_matrix_$(date +%Y%m%d).txt"

cd $directory

# Use the first two columns of the first file as the starting content of the output file
awk '{print $1, $2}' "$(find $directory -maxdepth 1 -type f | head -n 1)" > ${directory}/"NSEG_"${output_file}
# Iterate through each file in the directory and append the NSEG column to the outp
for file in *.indiv; do
    awk 'NR==1{print FILENAME}NR>1{print $4}' $file | sed 's!${directory}/!!g' | sed 's!.cnv.indiv!!g'|paste ${directory}/"NSEG_"$output_file - > ${directory}/"NSEG_"temp.txt && mv ${directory}/"NSEG_"temp.txt ${directory}/"NSEG_"${output_file}
done

#Use awk to calculate the sum of each column and print only those columns with a sum greater than 8
cat ${directory}/"NSEG_"${output_file} | awk '
{
    for (i=1; i<=NF; i++) sum[i] += $i
}
END {
    for (i=1; i<=NF; i++) {
        if (sum[i] > 8) cols[i] = 1
    }
}
{
    for (i=1; i<=NF; i++) {
        if (cols[i]) printf "%s ", $i
    }
    printf "\n"
}
' > ${directory}/"NSEG_"${output_file}"_filtered.txt"

cat ${directory}/"NSEG_"${output_file}"_filtered.txt" | awk 'NR==1{print}' | awk '
{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str = a[1,j]
        for(i=2; i<=NR; i++){
            str = str" "a[i,j]
        }
        print str
    }
}' > ${directory}/gene_filtered.txt

# Iterate through each file in the directory and append the KB column to the outp
awk '{print $1, $2}' "$(find $directory -maxdepth 1 -type f | head -n 1)" > ${directory}/"KB_"${output_file}
for file in *.indiv; do
    awk 'NR==1{print FILENAME}NR>1{print $5}' $file | sed 's!${directory}/!!g' | sed 's!.cnv.indiv!!g'| paste ${directory}/"KB_"$output_file - > ${directory}/"KB_"temp.txt && mv ${directory}/"KB_"temp.txt ${directory}/"KB_"${output_file}
done

cat ${directory}/gene_filtered.txt|awk '
NR==FNR { genes[$1]=1; next }
FNR==1 {
    for (i=1; i<=NF; i++) {
        if ($i in genes) {
            cols[i]=1;
            printf "%s\t", $i;
        }
    }
    printf "\n";
    next;
}
{
    for (i=1; i<=NF; i++) {
        if (i in cols) {
            printf "%s\t", $i;
        }
    }
    printf "\n";
}
' - ${directory}/"KB_"${output_file} > ${directory}/"KB_"${output_file}"_filtered.txt"

# Iterate through each file in the directory and append the KBAVG column to the outp
awk '{print $1, $2}' "$(find $directory -maxdepth 1 -type f | head -n 1)" > ${directory}/"KBAVG_"${output_file}
for file in *.indiv; do
    awk 'NR==1{print FILENAME}NR>1{print $6}' $file | sed 's!${directory}/!!g' | sed 's!.cnv.indiv!!g'| paste ${directory}/"KBAVG_"$output_file - > ${directory}/"KBAVG_"temp.txt && mv ${directory}/"KBAVG_"temp.txt ${directory}/"KBAVG_"${output_file}
done

cat ${directory}/gene_filtered.txt|awk '
NR==FNR { genes[$1]=1; next }
FNR==1 {
    for (i=1; i<=NF; i++) {
        if ($i in genes) {
            cols[i]=1;
            printf "%s\t", $i;
        }
    }
    printf "\n";
    next;
}
{
    for (i=1; i<=NF; i++) {
        if (i in cols) {
            printf "%s\t", $i;
        }
    }
    printf "\n";
}
' - ${directory}/"KBAVG_"${output_file} > ${directory}/"KBAVG_"${output_file}"_filtered.txt"

# Iterate through each file in the directory and append the COUNT column to the outp
awk '{print $1, $2}' "$(find $directory -maxdepth 1 -type f | head -n 1)" > ${directory}/"COUNT_"${output_file}
for file in *.indiv; do
    awk 'NR==1{print FILENAME}NR>1{print $7}' $file | sed 's!${directory}/!!g' | sed 's!.cnv.indiv!!g'| paste ${directory}/"COUNT_"$output_file - > ${directory}/"COUNT_"temp.txt && mv ${directory}/"COUNT_"temp.txt ${directory}/"COUNT_"${output_file}
done

cat ${directory}/gene_filtered.txt|awk '
NR==FNR { genes[$1]=1; next }
FNR==1 {
    for (i=1; i<=NF; i++) {
        if ($i in genes) {
            cols[i]=1;
            printf "%s\t", $i;
        }
    }
    printf "\n";
    next;
}
{
    for (i=1; i<=NF; i++) {
        if (i in cols) {
            printf "%s\t", $i;
        }
    }
    printf "\n";
}
' - ${directory}/"COUNT_"${output_file} > ${directory}/"COUNT_"${output_file}"_filtered.txt"

