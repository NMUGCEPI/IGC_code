##Ligand-receptor interactions among different cell populations were analyzed using CellPhoneDB.
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --counts-data=gene_name --output-path=out
cellphonedb plot dot_plot
cellphonedb plot heatmap_plot cellphonedb_meta.txt