from src.constants import MAF_FILES, ONCOKB_API_KEY,MAF_COLUMNS

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
		"poetry run python -m src.oncokb_annotator.MafAnnotator -i {input.maf} -o {output.maf} -b {ONCOKB_API_KEY}"


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


rule subset_maf_files:
	input:
		maf = rules.cat_maf_files.output.maf
	output:
		variants = "ontology/all_variants.csv"
	run:
		import pandas as pd
		df = pd.read_parquet(input.maf)
		df = df[MAF_COLUMNS]
		df = df=df.loc[df.iloc[:,-22:].any(axis=1)]
		df.to_csv(output.variants,index=False)
		


rule download_tcga_copyalt:
	output:
		copyalt = "source_data/tcga_gene_copyalt.gz"
	shell:
		"wget -O {output.copyalt}  https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"


rule write_tcga_combined_variants_table:
	input:
		copyalt = rules.download_tcga_copyalt.output.copyalt
	output:
		tcga_var_tbl = "source_data/TCGA_AllVarTypes_by_Sample.tsv",
		tcga_clinical = "source_data/TCGA_Clinical.tsv",
		tcga_cna = "source_data/TCGA_CNA.tsv",
		tcga_fusion = "source_data/TCGA_Fusions.tsv",
		tcga_maf = "source_data/TCGA_MAF.tsv"
	script:
		"src/Prep_TCGA_PanCan.R"

rule annotate_tcga_maf:
	input:
		tcga_maf = rules.write_tcga_combined_variants_table.output.tcga_maf,
		tcga_clinical = rules.write_tcga_combined_variants_table.output.tcga_clinical,
	output:
		tcga_maf = "source_data/tcga_annotated_maf.tsv"
	shell:
		"poetry run python -m src.oncokb_annotator.MafAnnotator -i {input.tcga_maf} -o {output.tcga_maf}  -c {input.tcga_clinical} -b {ONCOKB_API_KEY}"

rule annotate_fusion:
	input: 
		tcga_fusion = rules.write_tcga_combined_variants_table.output.tcga_fusion,
		tcga_clinical = rules.write_tcga_combined_variants_table.output.tcga_clinical
	output:
		tcga_fusion ="source_data/tcga_annotated_fusion.tsv"
	shell:
		"poetry run python -m src.oncokb_annotator.FusionAnnotator -i {input.tcga_fusion} -o {output.tcga_fusion} -c {input.tcga_clinical} -b {ONCOKB_API_KEY}"


rule annotate_cnv:
	input: 
		tcga_cna = rules.write_tcga_combined_variants_table.output.tcga_cna,
		tcga_clinical = rules.write_tcga_combined_variants_table.output.tcga_clinical
	output:
		tcga_cna ="source_data/tcga_annotated_cna.tsv"
	shell:
		"poetry run python -m src.oncokb_annotator.CnaAnnotator -i {input.tcga_cna} -o {output.tcga_cna} -c {input.tcga_clinical} -b {ONCOKB_API_KEY}"


rule annotate_clinical:
	input: 
		tcga_clinical = rules.write_tcga_combined_variants_table.output.tcga_clinical,
		tcga_cna = rules.annotate_cnv.output.tcga_cna,
		tcga_fusion = rules.annotate_fusion.output.tcga_fusion,
		tcga_maf = rules.annotate_tcga_maf.output.tcga_maf,

	output:
		tcga_clinical ="source_data/tcga_annotated_clinical.tsv"
	shell:
		'poetry run python -m src.oncokb_annotator.ClinicalDataAnnotator -i {input.tcga_clinical} -o {output.tcga_clinical} -a "{input.tcga_maf},{input.tcga_cna},{input.tcga_fusion}"'



rule generate_civic:
	output:
		civic = "source_data/civic_evidence.csv"
	shell:
		"poetry run python -m src.generate_civic {output.civic}"

rule make_oncokb_owl:
	input:
		annotate_maf = rules.annotate_tcga_maf.output.tcga_maf,
		fusion  = rules.annotate_fusion.output.tcga_fusion,
		cna  = rules.annotate_cnv.output.tcga_cna,
		clinical = rules.annotate_clinical.output.tcga_clinical,
	output:
		owl = "ontology/oncokb.owl"
	shell:
		"poetry run python -m src.ontology {input.annotate_maf} {input.fusion} {input.cna} {output.owl}"

rule make_oncokb_civic_owl:
	input:
		civic = rules.generate_civic.output.civic,
		owl = rules.make_oncokb_owl.output.owl
	output:
		owl = "ontology/oncokb_civic.owl"
	shell:
		"poetry run python -m src.add_civic_onto {input.owl} {input.civic} {output.owl}"




rule run_inference:
	input:
		tcga = rules.write_tcga_combined_variants_table.output.tcga_var_tbl,
		onto = rules.make_oncokb_civic_owl.output.owl
	output:
		therapies = "inference/tcga_samples.csv"
	threads: 96
	shell:
		"poetry run python -m src.run_inference {input.onto}  {input.tcga} {output.therapies} {threads}"




	
		
