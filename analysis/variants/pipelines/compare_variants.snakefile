import os
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

def getPath( relativePath ):
	return os.path.join( ROOT, relativePath )

vcfs = {
	#'10X': getPath( 'projects/wgs_10x/SVs/JK_HV31_3/outs/phased_variants.vcf.gz' ),
	'10X': getPath( 'users/gav/tmp/10X_phased_variants.vcf.gz' ),
	#'nanopore': '%s/projects/nanopore//snv-calling/longshot-chr22-build38.vcf' % ROOT
	'nanopore': getPath( 'users/gav/tmp/longshot-chr22-build38.vcf.gz' )
}

gnomAD = getPath( 'projects/reference/GRCh38/gnomad/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz' )

tools = {
	'bcftools': '/apps/well/bcftools/1.4.1/bin/bcftools',
	'tabix': '/apps/well/htslib/1.8-gcc5.4.0/bin/tabix',
	'qctool': 'qctool_v2.0.5',
	'bgenix': 'bgenix',
	'vcflib': '/apps/well/vcflib/20170323-gcc4.7.2/bin/'
}

rule all:
	input:
		vcfs = [ 'vcf/{platform}.normalised.vcf.gz'.format( platform = platform ) for platform in vcfs.keys() ],
		gnomAD = 'vcf/gnomAD.normalised.vcf.gz'

rule getSNPs:
	input:
		vcf = lambda x: vcfs[ x.platform ],
		mask = getPath( "shared/masks/20160622.strict_mask_GRCh38.whole_genome.bed" ),
		ref = getPath( "projects/reference/GRCh38/10X/refdata-GRCh38-2.1.0/fasta/genome.fa" )
	output:
		result = 'vcf/{platform}.normalised.vcf.gz',
		tbi = 'vcf/{platform}.normalised.vcf.gz.tbi',
		bgen = 'bgen/{platform}.normalised.bgen',
		bgi = 'bgen/{platform}.normalised.bgen.bgi'
	shell:
		"""
		{tools[bcftools]} view -r chr22:0- {input.vcf} | \
		{tools[vcflib]}/vcfbreakmulti | \
		{tools[vcflib]}/vcfallelicprimitives --tag-parsed SPLIT | \
		{tools[vcflib]}/vcfleftalign -r {input.ref} | \
		{tools[vcflib]}/vcfstreamsort | bgzip -c > {output.result}
		tabix -p vcf {output.result}
		{tools[qctool]} -g {output.result} -og {output.bgen} -bgen-bits 1 -bgen-compression zstd
		{tools[bgenix]} -g {output.bgen} -index
	"""

rule extractGnomAD:
	input:
		vcf = gnomAD,
		ref = "%s/projects/reference/GRCh38/10X/refdata-GRCh38-2.1.0/fasta/genome.fa" % ROOT
	output:
		vcf = 'vcf/gnomAD.normalised.vcf.gz',
		tbi = 'vcf/gnomAD.normalised.vcf.gz.tbi',
		freq = 'vcf/gnomeAD.commonPositionsAF_nfe.freq'
	shell:
		"""
		{tools[bcftools]} view -r chr22:0- {input.vcf} | 
		{tools[vcflib]}/vcfbreakmulti | \
		{tools[vcflib]}/vcfallelicprimitives --tag-parsed SPLIT | \
		{tools[vcflib]}/vcfleftalign -r {input.ref} | \
		{tools[vcflib]}/vcfstreamsort | bgzip -c > {output.vcf}
		{tools[tabix]} -p vcf {output.vcf}
		vcftools --gzvcf {output.vcf} --get-INFO AF_nfe --out {output.freq}
	"""

rule listAll:
	input:
		tenX = 'bgen/10X.split.pass.masked.snps.bgen.bgi',
		nanopore = 'bgen/nanopore.split.pass.masked.snps.bgen.bgi'
	output:
		result = 'snps/snps.sqlite'
	shell:
		"""
		rm -f {output.result}
		echo "ATTACH DATABASE '{input.tenX}' AS tenX; \
			ATTACH DATABASE '{input.nanopore}' AS nanopore ; \
			CREATE TABLE Variant( chromosome TEXT, position INT, alleleA TEXT, alleleB TEXT ); \
			CREATE TEMP TABLE tmp( chromosome TEXT, position INT, alleleA TEXT, alleleB TEXT ) ; \
			INSERT INTO tmp SELECT chromosome, position, allele1, allele2 FROM tenX.Variant;  \
			INSERT INTO tmp SELECT chromosome, position, allele1, allele2 FROM nanopore.Variant; \
			INSERT INTO Variant SELECT DISTINCT chromosome, position, allele1, allele2, CASE WHEN TV IS NULL THEN 0 ELSE 1 END, CASE WHEN NV IS NULL THEN 0 ELSE 1 END \
				FROM tmp \
				LEFT OUTER JOIN tenX.Variant TV ON TV.chromosome == chromosome AND TV.position == position AND TV.allele1 == allele1 AND TV.allele2 == allele2, \
				LEFT OUTER JOIN nanopore.Variant NV ON NV.chromosome == chromosome AND NV.position == position AND NV.allele1 == allele1 AND NV.allele2 == allele2 \
			;" | sqlite3 {output.result}
	"""
