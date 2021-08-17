#Two ways for effective genome size calculation (from deeptools)

#=========================================================================================
#1. Calculated by non-N bases (These values only appropriate if multimapping reads are included)
#=========================================================================================
#faCount Saccharomyces_cerevisiae.EF4.67.dna.toplevel_chrM.fa -summary -dinuc > Saccharomyces_cerevisiae.EF4.67.dna.toplevel_chrM.fa.summary
#faCount Zea_mays.B73_RefGen_v4.dna.toplevel_chrM.fa -summary -dinuc > Zea_mays.B73_RefGen_v4.dna.toplevel_chrM.fa.summary
#faCount Homo_sapiens.GRCh38.dna.primary_assembly_chrM.fa -summary -dinuc  > Homo_sapiens.GRCh38.dna.primary_assembly_chrM.fa.summary
#faCount Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa -summary -dinuc > Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa.summary
#faCount Caenorhabditis_elegans.WBcel235.dna.toplevel_chrM.fa -summary -dinuc > Caenorhabditis_elegans.WBcel235.dna.toplevel_chrM.fa.summary
#faCount Danio_rerio.GRCz11.dna.primary_assembly_chrM.fa  -summary -dinuc  > Danio_rerio.GRCz11.dna.primary_assembly_chrM.fa.summary
#faCount Arabidopsis_thaliana.TAIR10.dna.toplevel_chrM.fa  -summary -dinuc  >Arabidopsis_thaliana.TAIR10.dna.toplevel_chrM.fa.summary
#faCount Drosophila_melanogaster.BDGP6.22.dna.toplevel_chrM.fa -summary -dinuc > Drosophila_melanogaster.BDGP6.22.dna.toplevel_chrM.fa.summary
#faCount Plasmodium_falciparum.EPr1.dna.toplevel_chrM.fa -summary -dinuc > Plasmodium_falciparum.EPr1.dna.toplevel_chrM.fa.summary

sacCer3=12157105
zm3=2135083061
pfa2=23292622
danRer11=1373471384
hg38=3099750718
tair10=119667750
mm10=2730871774
dm6=143726002
ce11=100286401

#2. The number of regions in the genome that are uniquely mappable (If there’s any MAPQ filter applied).

#Read_length	hg38	mm10	dm6	danRer11	ce11
#50	2701495761	2308125349	125464728	1195445591	95159452
#75	2747877777	2407883318	127324632	1251132686	96945445
#100	2805636331	2467481108	129789873	1280189044	98259998
#150	2862010578	2494787188	129941135	1312207169	98721253
#200	2887553303	2520869189	132509163	1321355241	98672758