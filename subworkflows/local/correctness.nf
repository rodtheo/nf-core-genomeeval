include { ALE } from '../../modules/local/ALE'
include { REAPR } from '../../modules/local/REAPR'
include { HEADER_FASTA_REAPR } from '../../modules/local/header_fasta_reapr'

workflow CORRECTNESS_ASM { 

    take:
	asm 
    bam        // queue channel: [ sample_id, file(bam_file) ]
    bai    // queue channel: [ sample_id, file(bai_file) ]

    main:
    // Quality Check input reads
    // READ_QC ( reads )

    
	// reference.map{ meta, assemblies -> {
	// 	def meta_new = meta.clone();
	// 	meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
	// 	[ meta_new, assemblies['pri_asm'] ]
	//  } }.set{ reference_ch }
	

	// bam.combine(bai).map{ meta_bam, bam, meta_bai, bai -> {
	// 	def meta_cor_reads = meta_reads.clone();
	// 	meta_cor_reads.put("id", meta_asm.get("id"));
	// 	[ meta_cor_reads, reads, meta_asm, asm ]
	// }
	// }.set{ inputs_reads_asm_ch }
	// inputs_reads_asm_ch.view{ "ASSEMBLIES+READS COMBINED CHANNEL: $it"}

	left = asm.map{ meta_asm, asm -> [meta_asm['id'], meta_asm, asm] }
	right = bam.map{ meta_bam, bam -> [meta_bam['id'], meta_bam, bam] }

	joined_ch = left.join(right, failOnMismatch: true).view{ "### JOINED: $it"}

	joined_ch.multiMap{ meta_id, meta_asm, asm, meta_bam, bam ->
		asm: [meta_asm, asm]
		bam: [meta_bam, bam] }.set{ input_ale_ch }

	// SAMTOOLS_SORT( input_samtools_sort_ch.bam,  input_samtools_sort_ch.ref )

	// Align reads to reference
    Channel.empty()
        .set { ale_out_ch }

	ALE ( input_ale_ch )

	ale_out_ch.mix(ALE.out.ale).set{ ale_out_ch }
	
	Channel.empty()
        .set { reapr_out_ch }

	HEADER_FASTA_REAPR ( input_ale_ch.asm )

	REAPR ( HEADER_FASTA_REAPR.out.asm, input_ale_ch.bam )
	
	reapr_out_ch.mix( REAPR.out.reapr_summary ).set{ reapr_out_ch }

	// Execute QUAST without a reference and a GFF
	// QUAST ( reference_ch, [[], []], [[], []] );

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
    ale = ale_out_ch   // queue channel: [ sample_id, file(bam_file) ]
	reapr = reapr_out_ch
	versions = ALE.out.versions


}