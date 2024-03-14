
include { BOWTIE_BUILD } from '../../modules/nf-core/bowtie/build/main'
include { BOWTIE_ALIGN } from '../../modules/nf-core/bowtie/align/main'
include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'




workflow READ_MAPPING {

    take:
    reads        // queue channel; [ sample_id, [ file(read1), file(read2) ] ]
    reference    // file( "path/to/reference" )

    main:
    // Quality Check input reads
    // READ_QC ( reads )

    // Align reads to reference
    Channel.empty()
        .set { aligned_reads_ch }

	reference.map{ meta, assemblies -> {
		def meta_new = meta.clone();
		meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
		[ meta_new, assemblies['pri_asm'] ]
	 } }.set{ reference_ch }
	

	reference_ch.combine(reads).map{ meta_asm, asm, meta_reads, reads -> {
		def meta_cor_reads = meta_reads.clone();
		meta_cor_reads.put("id", meta_asm.get("id"));
		[ meta_cor_reads, reads, meta_asm, asm ]
	}
	}.set{ inputs_reads_asm_ch }
	inputs_reads_asm_ch.view{ "ASSEMBLIES+READS COMBINED CHANNEL: $it"}


    if( params.aligner == 'bowtie' ){
		inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ asm ] }.set{ reference_bowtie_ch }
		BOWTIE_BUILD ( reference_bowtie_ch )

		inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_cor_reads, reads ] }.set{ reads_ch_cor_meta }
		BOWTIE_ALIGN ( reads_ch_cor_meta, BOWTIE_BUILD.out.index )
		
		aligned_reads_ch.mix( BOWTIE_ALIGN.out.bam )
            .set { aligned_reads_ch }
    } else if ( params.aligner == 'bwa' ) {
		inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_asm, asm ] }.set{ reference_bwa_ch }
        BWA_INDEX ( reference_bwa_ch )

		inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_cor_reads, reads ] }.set{ reads_ch_cor_meta }
        BWA_MEM ( reads_ch_cor_meta, BWA_INDEX.out.index, true )

		aligned_reads_ch.mix( BWA_MEM.out.bam )
            .set { aligned_reads_ch }
    }
    aligned_reads_ch.view{"ALIGNED READS: $it"}

	// Empty channel to store sorted bams
    Channel.empty()
        .set { bams_ch } // queue channel: [ sample_id, file(bam_file)]

	// Empty channel to store sorted bams index
    Channel.empty()
        .set { bams_index_ch } // queue channel: [ sample_id, file(bam_index_file)]

	left = inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [meta_asm['id'], meta_asm, asm] }
	right = aligned_reads_ch.map{ meta, bam -> [meta['id'], meta, bam] }

	joined_asm = left.join(right, failOnMismatch: true).view{ "### JOINED: $it"}

	joined_asm.multiMap{ asm_id, meta_asm, asm, meta_reads, bam ->
		bam: [meta_reads, bam]
		ref: [meta_asm, asm] }.set{ input_samtools_sort_ch }

	SAMTOOLS_SORT( input_samtools_sort_ch.bam,  input_samtools_sort_ch.ref )

	bams_ch.mix( SAMTOOLS_SORT.out.bam )
            .set { bams_ch }

	SAMTOOLS_INDEX( SAMTOOLS_SORT.out.bam )

	bams_index_ch.mix( SAMTOOLS_INDEX.out.bai )
            .set { bams_index_ch }

    emit:
	asm = input_samtools_sort_ch.ref	// queue channel: [ sample_id, file(asm_fasta) ]
    bam = bams_ch   // queue channel: [ sample_id, file(bam_file) ]
	bai = bams_index_ch // queue channel: [ sample_id, file(bai_file) ]

}