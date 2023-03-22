for file in star_code/*_mapped.tsv
do echo starcode -i $file --print-clusters -s -d 1 -o ${file%%.*}_sc_out.tsv
starcode -i $file --print-clusters -s -d 1 -o ${file%%.*}_sc_out.tsv
done