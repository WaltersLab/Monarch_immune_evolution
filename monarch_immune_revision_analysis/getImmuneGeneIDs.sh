grep "Immune" Genestats_Filtered_2021-05-15.csv | cut -d',' -f 1 | sed  's/\"//g' > Genes_Filtered.txt

