#!/bin/bash

gwas_info_file="../C100001554_And_Friends_3Metabolites.csv"

# Output file
output="chr19_snp_table.tsv"

# Print header first, then append pipeline output
{
    echo -e "snp\tntraits\tminmaf\tminss\tmaxss"

    awk -F',' '
    NR==1 {
        for(i=1; i<=NF; i++) {
            if($i=="raw_data_path") file_col=i;
            if($i=="snp") snp_col=i;
            if($i=="af") maf_col=i;
            if($i=="chrom") chrom_col=i;
            if($i=="sample_size") ss_col=i;
        }
        next;
    }
    NR>1 {
        print $file_col, $snp_col, $maf_col, $chrom_col, $ss_col;
    }
    ' "$gwas_info_file" \
    | while read trait_file snp_col_name maf_col_name chrom_col_name ss_col_name; do

        [[ "$trait_file" == *.gz ]] && cat_cmd="zcat" || cat_cmd="cat"

        header=$($cat_cmd "$trait_file" | head -1)
        IFS=$'\t' read -r -a cols <<< "$header"

        snp_col_idx=""
        maf_col_idx=""
        chrom_col_idx=""
        ss_col_idx=""
        for i in "${!cols[@]}"; do
            [[ "${cols[$i]}" == "$snp_col_name" ]] && snp_col_idx=$((i+1))
            [[ "${cols[$i]}" == "$maf_col_name" ]] && maf_col_idx=$((i+1))
            [[ "${cols[$i]}" == "$chrom_col_name" ]] && chrom_col_idx=$((i+1))
            [[ "${cols[$i]}" == "$ss_col_name" ]] && ss_col_idx=$((i+1))
        done

        trait=$(basename "$trait_file" | sed 's/\..*//')

        $cat_cmd "$trait_file" | awk -F'\t' -v trait="$trait" \
            -v snp_col="$snp_col_idx" -v maf_col="$maf_col_idx" \
            -v chrom_col="$chrom_col_idx" -v ss_col="$ss_col_idx" \
            'NR > 1 && $chrom_col == "19" { print $snp_col "\t"trait"\t"$maf_col"\t"$ss_col }'

    done \
    | awk -F'\t' '
    {
      snp=$1
      trait=$2
      maf=$3
      ss=$4
      traits[snp][trait]=1
      if( (snp in minmaf) == 0 || maf < minmaf[snp]) minmaf[snp] = maf
      if( (snp in minss) == 0 || ss < minss[snp]) minss[snp] = ss
      if( (snp in maxss) == 0 || ss > maxss[snp]) maxss[snp] = ss
    }
    END {
      for(snp in traits) {
        ntraits = 0
        for(trait in traits[snp]) ntraits++
        printf "%s\t%s\t%.5f\t%s\t%s\n", snp, ntraits, minmaf[snp], minss[snp], maxss[snp]
      }
    }
    '

} > "$output"
