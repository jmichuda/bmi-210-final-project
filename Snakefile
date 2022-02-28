from src.constants import MAF_FILES, ONCOKB_API_KEY


# to do: 
# add raw, unzip and annotated to suffix instead of prefix

rule download_data:
	output:
		all_maf_files = "source_data/gdc_download.tar.gz"
	run:
		import gdown
		url = "https://drive.google.com/uc?id=1RYKzSFCo7cCUo-FQye-wcoF0rNMqGmgq"
		gdown.download(url, output.all_maf_files, quiet=False)


rule unzip:
	input: 
		all_maf_files = rules.download_data.output.all_maf_files
	output:
		maf_files = expand("maf_files/{maf_file}.maf.gz", maf_file = MAF_FILES)
	shell:
		"mkdir -p maf_files/raw/ && "
		"tar -xf {input.all_maf_files} -C maf_files/"

rule unzip_maf_files:
	input: 
		gz_maf = "maf_files/{maf_file}.maf.gz"
	output:
		maf = "maf_files/{maf_file}.maf"
	shell:
		"mkdir -p maf_files/unzip/ && "
		"gunzip -c {input.gz_maf} > {output.maf}"

rule annotate_maf_files:
	input:
		maf = rules.unzip_maf_files.output.maf
	output:
		maf = "maf_files/{maf_file}_annotated.maf"
	shell:
		"python -m src.MafAnnotator -i {input.maf} -o {output.maf} -b {ONCOKB_API_KEY}"


rule cat_maf_files:
	input:
		all_unzip = expand(rules.annotate_maf_files.output.maf, maf_file = MAF_FILES)
	output:
		maf = "maf_files/concatenated_maf.parquet"
	run:
		import pandas as pd
		mafs = []
		for i in input.all_unzip:
			mafs.append(pd.read_csv(i,sep = "\t",index_col=False,low_memory=False))
		pd.concat(mafs).to_parquet(output.maf)


rule subset_level_1:
	input:
		maf = rules.cat_maf_files.output.maf
	output:
		variants = "ontology/level_1.csv"
	run:
		import pandas as pd
		df = pd.read_parquet(input.maf)
		df = df.loc[df['LEVEL_1'].notna()]
		df = df[["Hugo_Symbol","HGVSp_Short","LEVEL_1"]].drop_duplicates()
		df.to_csv(output.variants,index=False)
		

rule make_owl:
	input:
		level_1 = rules.subset_level_1.output.variants
	output:
		ontology = "ontology/oncokb.owl"
	shell:
		"python -m src.ontology {input.level_1} {output.ontology}"


rule download_tcga_copyalt:
	output:
		copyalt = "source_data/tcga_gene_copyalt.gz"
	shell:
		"wget -O {output.copyalt}  https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"


rule write_tcga_combined_variants_table:
	input:
		copyalt = rules.download_tcga_copyalt.output.copyalt
	output:
		tcga_var_tbl = "source_data/TCGA_AllVarTypes_by_Sample.tsv.gz"
	script:
		"src/Prep_TCGA_PanCan.R"


	
		
