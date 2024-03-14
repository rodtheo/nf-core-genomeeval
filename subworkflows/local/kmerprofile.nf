include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq/main'
include { MERYL_COUNT as MERYL_COUNT_READS_01;
          MERYL_COUNT as MERYL_COUNT_GENOME_01  } from '../../modules/nf-core/meryl/count/main'
include { MERYL_GREATER_THAN } from '../../modules/local/meryl_greater'
include { MERYL_HISTOGRAM as MERYL_HISTOGRAM_READS_PRE;
          MERYL_HISTOGRAM as MERYL_HISTOGRAM_GENOME_PRE } from '../../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2 as GENOMESCOPE2_PRE } from '../../modules/nf-core/genomescope2/main'
include { MERFIN_COMPLETENESS } from '../../modules/local/merfin_completeness'
include { MERFIN_HIST as MERFIN_HIST_EVALUATE_POLISH;
          MERFIN_HIST as MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY } from '../../modules/local/merfin_hist'

process GET_PEAK {
    input:
    tuple val(meta), path(genomescope2model)

    output:
    path("*_peak.txt"), emit: peak
// cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))' > file

    // peak="\$(cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))')"
    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))' > ${prefix}_peak.txt
    """
}

workflow KMER_PROFILE {

	take:
    reads        // queue channel; [ sample_id, [ file(read1), file(read2) ] ]
    reference    // file( "path/to/reference" )

    main:

	// (0) Prepare nf empty variables with will receive results
    Channel.empty()
        .set { aligned_reads_ch }

	Channel.empty()
        .set { res_completeness_ch }

	Channel.empty()
        .set { res_logs_ch }
		
    
	// (1) Generate k-mers count distribution

	// (1.1) Count operation: count the occurrences of canonical k-mers in high-accuracy reads
	
	// Create a channel where paired-end data is mentioned as single-end
    // This is necessary in order to concatenate paired-end in a single file by cat/fastq nf-core module
	
    reads_ch_single = reads.map { meta, reads -> {
		def meta_new = meta.clone();
		meta_new.put('single_end', true);
		[meta_new, reads] }
	}
	reads_ch_single.view{ "++++++++++++ READS: $it" }
            
    CAT_FASTQ (
        reads_ch_single
    )
	
	MERYL_COUNT_READS_01 (
        CAT_FASTQ.out.reads
    )

	// (1.2) Rename output as meryldb - SKIPPED BECAUSE WE JOIN THE READS?
	// (1.3) Operation Union-sum: return k-mers that occur in any input, set the count to the sum of the counts - SKIPPED BECAUSE WE JOIN THE READS?
	// (1.4) Rename it as merged_meryldb - SKIPPED BECAUSE WE JOIN THE READS?

	// (1.5-ROD) Keep only k-mers with frequency > 1
	MERYL_GREATER_THAN (
        MERYL_COUNT_READS_01.out.meryl_db
    )

	// (1.6) Generate histogram 
	MERYL_HISTOGRAM_READS_PRE (
        // MERYL_COUNT_READS_01.out.meryl_db.mix(MERYL_COUNT_GENOME_01.out.meryl_db)
        // MERYL_GREATER_THAN.out.meryl_db
		MERYL_COUNT_READS_01.out.meryl_db
    )

	// (1.7) Genome profiling with GenomeScope2
	GENOMESCOPE2_PRE (
        MERYL_HISTOGRAM_READS_PRE.out.hist
    )

    peak_out = GET_PEAK (
        GENOMESCOPE2_PRE.out.model
    )

    peak_ch_val = peak_out.map{ it.text.trim() }.first()

	// (1.9) Genome assembly k-mer count
	reference.map{ meta, assemblies -> {
		def meta_new = meta.clone();
		meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
		[ meta_new, assemblies['pri_asm'] ]
	 } }.set{ reference_ch }

	MERYL_COUNT_GENOME_01 (
        reference_ch
    )

	// // (1.8) K-mer based evaluation with Merfin (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9745813/)
	MERFIN_COMPLETENESS (
        MERYL_COUNT_GENOME_01.out.meryl_db,
        MERYL_GREATER_THAN.out.meryl_db.first(),
        GENOMESCOPE2_PRE.out.lookup_table.first(),
        peak_ch_val
    )

	res_completeness_ch.mix(MERFIN_COMPLETENESS.out.completeness).set{ res_completeness_ch }	

    MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY (
        reference_ch,
        MERYL_GREATER_THAN.out.meryl_db.first(),
        GENOMESCOPE2_PRE.out.lookup_table.first(),
        peak_ch_val
    )

	res_logs_ch.mix( MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY.out.log).set{ res_logs_ch }


	




    

	// reference.map{ meta, assemblies -> {
	// 	def meta_new = meta.clone();
	// 	meta_new.put("id", meta.get("id")+"_"+meta.get("build"));
	// 	[ meta_new, assemblies['pri_asm'] ]
	//  } }.set{ reference_ch }
	

	// reference_ch.combine(reads).map{ meta_asm, asm, meta_reads, reads -> {
	// 	def meta_cor_reads = meta_reads.clone();
	// 	meta_cor_reads.put("id", meta_asm.get("id"));
	// 	[ meta_cor_reads, reads, meta_asm, asm ]
	// }
	// }.set{ inputs_reads_asm_ch }



    emit:
	res = aligned_reads_ch
	merfin_completeness = res_completeness_ch
	merfin_logs = res_logs_ch
    // quast = QUAST.out.results   // queue channel: [ sample_id, file(bam_file) ]
	// quast_tsv = QUAST.out.tsv
	// busco = BUSCO.out.busco_dir
	// busco_short_summaries_txt = BUSCO.out.short_summaries_txt
	// versions = QUAST.out.versions
}