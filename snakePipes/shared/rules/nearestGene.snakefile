sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
def get_outdir(folder_name):
    sample_name = re.sub('_sampleSheet.[a-z]{3}$','',os.path.basename(sampleSheet))
    return("{}_{}".format(folder_name, sample_name))

def calc_ng_input(pipeline): 
    if pipeline in ['chip-seq','ATAC-seq']:
        return("CSAW_{}".format(sample_name)+"/Filtered.results.{change_dir}.bed")
    elif pipeline in 'WGBS':
        return("{}".format(get_outdir("metilene_out"))+"/metilene.limma_filtered.{change_dir}.bed")
    else:
        return([])

def calc_fo(pipeline):
    if pipeline in ['chip-seq','ATAC-seq']:
        return("15")
    elif pipeline in 'WGBS':
        return("17")
    else:
       return("")

if pipeline in ['chip-seq','ATAC-seq']:
    change_direction = ["UP","DOWN","MIXED"]
elif pipeline in 'WGBS':
    change_direction= ["UP","DOWN"]

rule get_nearest_transcript:
    input:
        bed=calc_ng_input(pipeline)
    output:
        annotated_bed=temp("AnnotatedResults_{}".format(sample_name)+"/Filtered.results.{change_dir}_withNearestTranscript.bed")
    params:
        genes_bed=genes_bed,
        field_offset=calc_fo(pipeline) 
    log:
        err= "AnnotatedResults_{}".format(sample_name)+"/logs/bedtools_closest.{change_dir}.err",
    conda: CONDA_RNASEQ_ENV
    shell: "if [ -r {input.bed} ]; then bedtools closest -D b -a <( bedtools sort -i {input.bed} ) -b <( bedtools sort -i {params.genes_bed} ) | cut -f1-{params.field_offset},$(( {params.field_offset} + 1 ))-$(( {params.field_offset} + 4 )),$(( {params.field_offset} + 6 )),$(( {params.field_offset} + 13 )) > {output.annotated_bed};fi 2> {log.err}"

rule get_nearest_gene:
    input:
        bed="AnnotatedResults_{}".format(sample_name)+"/Filtered.results.{change_dir}_withNearestTranscript.bed",
        t2g="Annotation/genes.filtered.t2g",
        gene_symbol="Annotation/genes.filtered.symbol"
    output:
        annotated_bed="AnnotatedResults_{}".format(sample_name)+"/Filtered.results.{change_dir}_withNearestGene.txt"
    params:
        pipeline=pipeline
    log:
        err="AnnotatedResults_{}".format(sample_name)+"/logs/nearestGene.{change_dir}.err",
        out="AnnotatedResults_{}".format(sample_name)+"/logs/nearestGene.{change_dir}.out"
    conda: CONDA_RNASEQ_ENV
    script: "../rscripts/nearestGene.R"
