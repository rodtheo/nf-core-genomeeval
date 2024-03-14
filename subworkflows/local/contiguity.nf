include { QUAST } from '../../modules/nf-core/quast/main'

workflow CONTIGUITY_ASM {

    take:
    // reads        // queue channel; [ sample_id, [ file(read1), file(read2) ] ]
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
	

	// reference_ch.combine(reads).map{ meta_asm, asm, meta_reads, reads -> {
	// 	def meta_cor_reads = meta_reads.clone();
	// 	meta_cor_reads.put("id", meta_asm.get("id"));
	// 	[ meta_cor_reads, reads, meta_asm, asm ]
	// }
	// }.set{ inputs_reads_asm_ch }
	// inputs_reads_asm_ch.view{ "ASSEMBLIES+READS COMBINED CHANNEL: $it"}

	// Execute QUAST without a reference and a GFF
	QUAST ( reference_ch, [[], []], [[], []] );

    // if( params.aligner == 'bowtie' ){
	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ asm ] }.set{ reference_bowtie_ch }
	// 	BOWTIE_BUILD ( reference_bowtie_ch )

	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_cor_reads, reads ] }.set{ reads_ch_cor_meta }
	// 	BOWTIE_ALIGN ( reads_ch_cor_meta, BOWTIE_BUILD.out.index )

	// 	aligned_reads_ch.mix( BOWTIE_ALIGN.out.bam )
    //         .set { aligned_reads_ch }
    // } else if ( params.aligner == 'bwa' ) {
	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_asm, asm ] }.set{ reference_bwa_ch }
    //     BWA_INDEX ( reference_bwa_ch )

	// 	inputs_reads_asm_ch.map{ meta_cor_reads, reads, meta_asm, asm -> [ meta_cor_reads, reads ] }.set{ reads_ch_cor_meta }
    //     BWA_MEM ( reads_ch_cor_meta, BWA_INDEX.out.index, true )

	// 	aligned_reads_ch.mix( BWA_MEM.out.bam )
    //         .set { aligned_reads_ch }
    // }
    // aligned_reads_ch.view()

    emit:
	quast = QUAST.out.results   // queue channel: [ sample_id, file(bam_file) ]
	quast_tsv = QUAST.out.tsv
	versions = QUAST.out.versions

}