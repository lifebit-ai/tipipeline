#!/usr/bin/env nextflow
nextflow.preview.dsl=2
// general parameters
params.out_dir = ""
params.part=1
// Input path parameters Part 1
params.counts = "" // specify a path to a raw count matrix from cellranger 3.0
params.list_markers = "" // e.g : "MS4A1,GNLY,CD3E,CD14,FCER1A,FCGR3A,LYZ,PPBP,CD8A"
// Part 2
params.clusters="" // e.g : "1,7,9"
params.part1_object = "" // path to the .rds file from the first part of the pipeline

// R scripts parameters
// Part 1
params.prepro_script = "${baseDir}/scripts/pre_processing.R"
params.qc_script = "${baseDir}/scripts/qc_report.Rmd"
params.norm_script = "${baseDir}/scripts/norm.R"
params.clustering_script = "${baseDir}/scripts/clustering_report.Rmd"
params.save_script = "${baseDir}/scripts/save_part1.R"
// Part 2
params.cluster_selection = "${baseDir}/scripts/cluster_selection.R"
params.dimred_script = "${baseDir}/scripts/dynred_wrap.R"
params.dynverse_script = "${baseDir}/scripts/final_report.Rmd"


//#######################
// Process definition
//#######################
process pre_processing {
  label "r"
  container "qrouchon/r-scrnaseq-fqc"
  //publishDir "$outDir", mode: 'copy'
  input:
  file script
  file cmat
  output:
  file("preprocess.rds")
  script:
  """
  Rscript ${script} ${cmat}
  """
}

process filtering_qc_r {
  label "r"
  container "qrouchon/r-scrnaseq-fqc"
  publishDir "${params.out_dir}", mode: 'copy', pattern: 'qc_report.html'
  input:
  file script
  file rawMat
  output:
  file("filtered.rds")
  file("qc_report.html")
  script:
  """
  cp -L ${script} qc_report_nosymlink.Rmd
  cp -L ${rawMat} nosymlink.rds
  R -e 'Sys.setenv(RSTUDIO_PANDOC = "/usr/bin/pandoc");rmarkdown::render("qc_report_nosymlink.Rmd",output_file = "qc_report.html",params = list(input="nosymlink.rds"))'
  rm qc_report_nosymlink.Rmd nosymlink.rds
  """
}

process norm_r {
  label "r"
  container "qrouchon/r-norm"
  //publishDir "${params.out_dir}", mode: 'copy'
  input:
  file script
  file fmat
  output:
  file("norm.rds")
  script:
  """
  Rscript ${script} ${fmat}
  """
}

process clustering_r {
  label "r"
  container "qrouchon/r-norm-seurat"
  publishDir "${params.out_dir}", mode: 'copy', pattern: 'umap.html'
  input:
  file script
  file nmat
  val markers
  output:
  file("seurat_cluster.rds")
  file("umap.html")
  script:
  """
  cp -L ${script} clustering_report_nosymlink.Rmd
  cp -L ${nmat} nosymlink.rds
  R -e 'Sys.setenv(RSTUDIO_PANDOC = "/usr/bin/pandoc");rmarkdown::render("clustering_report_nosymlink.Rmd",output_file = "umap.html",params = list(input="nosymlink.rds",markers="${markers}"))'
  rm clustering_report_nosymlink.Rmd nosymlink.rds
  """
}

process save_part1 {
  label "r"
  container "qrouchon/r-norm-seurat"
  publishDir "${params.out_dir}", mode: 'copy'
  input:
  file script
  file filtered
  file normalized
  file cluster_info
  output:
  file("part1.rds")
  script:
  """
  Rscript ${script} ${filtered} ${normalized} ${cluster_info}
  """
}

process select_cluster {
  label "r"
  container "qrouchon/r-scrnaseq-fqc"
  //publishDir "${params.out_dir}", mode: 'copy'
  input:
  file script
  file part1
  val selected
  output:
  file("clusters_filtered.rds")
  file("clusters_norm.rds")
  script:
  """
  Rscript ${script} ${part1} ${selected}
  """
}

process dimred_wrap_expression {
  label "r"
  container "qrouchon/r-dynverse-lite"
  //publishDir "${params.out_dir}", mode: 'copy'
  input:
  file script
  file rdsfilt
  file rdsnorm
  output:
  file("ti_ready_dataset.h5")
  script:
  """
  Rscript ${script} ${rdsfilt} ${rdsnorm}
  """
}

process paga {
  echo
  label "r"
  container "dynverse/ti_paga"
  //publishDir "${params.out_dir}", mode: 'copy'
  input:
  file fh5
  output:
  file("paga_out.h5")
  script:
  """
  cp -L ${fh5} noLink.h5
  python /code/run.py --dataset noLink.h5 --output paga_out.h5
  rm noLink.h5
  """
}

process slingshot {
  label "r"
  container "dynverse/ti_slingshot"
  //publishDir "${params.out_dir}", mode: 'copy'
  input:
  file fh5
  output:
  file("slingshot_out.h5")
  script:
  """
  cp -L ${fh5} noLink.h5
  Rscript /code/run.R --dataset noLink.h5 --output slingshot_out.h5
  rm noLink.h5
  """
}

process dynplot {
  label "r"
  container "qrouchon/r-dynverse-lite"
  publishDir "${params.out_dir}", mode: 'copy'
  input:
  file script
  file rh5p
  file rh5s
  output:
  file("final_report.html")
  script:
  """
  cp -L ${script} finalreport_nosymlink.Rmd
  cp -L ${rh5p} paga.h5
  cp -L ${rh5s} slingshot.h5
  R -e 'Sys.setenv(RSTUDIO_PANDOC = "/usr/bin/pandoc");rmarkdown::render("finalreport_nosymlink.Rmd",output_file = "final_report.html",params = list(paga = "paga.h5",slingshot="slingshot.h5"))'
  rm finalreport_nosymlink.Rmd paga.h5 slingshot.h5
  """
}

//#################
// Workflow definition
//#################

log.info """\
         J&J tiPipeline for scRNAseq
         ===================================
         outdir       : ${params.out_dir}
         """
         .stripIndent()


if ( params.part == 1 ) {

  counts = Channel.fromPath(params.counts)
  prepro_script = Channel.fromPath(params.prepro_script)
  qc_script = Channel.fromPath(params.qc_script)
  norm_script = Channel.fromPath(params.norm_script)
  clustering_script = Channel.fromPath(params.clustering_script)
  save_script = Channel.fromPath(params.save_script)

  preproData = pre_processing(prepro_script,counts)
  filtData = filtering_qc_r(qc_script,preproData)
  normData = norm_r(norm_script,filtData[0])
  seuratObj = clustering_r(clustering_script,normData,params.list_markers)
  save_part1(save_script, filtData[0],normData,seuratObj[0])

} else if ( params.part == 2 ) {

  part1_objects = Channel.fromPath(params.part1_object)
  cluster_sel = Channel.fromPath(params.cluster_selection)
  dimred_script = Channel.fromPath(params.dimred_script)
  dynverse_script = Channel.fromPath(params.dynverse_script)

  selectedData = select_cluster(cluster_sel,part1_objects,params.clusters)
  dynObject = dimred_wrap_expression(dimred_script, selectedData[0], selectedData[1])
  dynplot(dynverse_script,paga(dynObject),slingshot(dynObject))

} else {
  error "you should choose a pipeline part to run (1 or 2)"
}
