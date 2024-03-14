include { QUAST } from '../../modules/nf-core/quast/main'
include { BUSCO } from '../../modules/nf-core/busco/main'
include { COMPLEASM } from '../../modules/local/compleasm/main'

workflow COMPLETENESS_ASM {

	take:
    reads        // queue channel; [ sample_id, [ file(read1), file(read2) ] ]
    reference    // file( "path/to/reference" )

    main:
    // Quality Check input reads
    // READ_QC ( reads )

    // Align reads to reference
    Channel.empty()
        .set { results_dir_ch }

	Channel.empty()
        .set { results_summary_ch }

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
	// inputs_reads_asm_ch.view{ "ASSEMBLIES+READS COMBINED CHANNEL: $it"}

	// Execute QUAST without a reference and a GFF
	// QUAST ( reference_ch, [[], []], [[], []] )

	lineages_ch = reference_ch.map{ meta_asm, asm -> meta_asm['busco_lineages'][0]}.first()
	lineages_ch.view{ "LINEAGES: $it" }


	if ( params.busco == true){
		BUSCO ( reference_ch, 'genome', lineages_ch, [], [] )

		results_dir_ch.mix (BUSCO.out.busco_dir ).set{ results_dir_ch }
		results_summary_ch.mix (BUSCO.out.short_summaries_txt).set{ results_summary_ch }

		

	} else {
		COMPLEASM ( reference_ch, lineages_ch )

		results_dir_ch.mix (COMPLEASM.out.busco_dir ).set{ results_dir_ch }
		results_summary_ch.mix (COMPLEASM.out.short_summaries_txt).set{ results_summary_ch }
		
	}	





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
	busco = results_dir_ch
	// busco_short_summaries_txt = BUSCO.out.short_summaries_txt
	busco_short_summaries_txt = results_summary_ch
	// versions = QUAST.out.versions
}