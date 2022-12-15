#!/usr/bin/env nextflow
// USAGE: nextflow run b2btools.nf -resume -config b2btools.config -with-report -with-dag pipeline.PNG

params.targetSequences = "$launchDir/input_small.fasta"
params.groupBy = 10

allSequences = Channel.fromPath(params.targetSequences)

sequencesGrouped = allSequences.splitFasta( by: params.groupBy, file: true )

process createMultipleSequenceAlignment {
    input:
    path sequences

    output:
    path "*.msa", emit: multipleSequenceAlignment

    script:
    """
    clustalo -i $sequences -o ${sequences}.msa
    """
}

process buildPhilogeneticTree {
    input:
    path multipleSequenceAlignment

    output:
    path "*.tree", emit: tree

    script:
    """
    FastTree $multipleSequenceAlignment > ${multipleSequenceAlignment}.tree
    """
}

process renderTree {
    input:
    path tree

    output:
    path "*.png", emit: tree_plot

    script:
    """
    #!/usr/bin/env Rscript

    dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
    .libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path

    install.packages("ape")  # install like always
    library(ape)  # use library like always

    mytr <- read.tree("$tree")
    png("tree.png")
    plot(mytr)
    dev.off()
    """
}


process predictBiophysicalFeatures {
    tag "${sequences.baseName}"

    input:
    path sequences

    output:
    path '*.json', emit: predictions

    script:
    """
    python -m b2bTools -agmata -dynamine -disomine -efoldmine -file $sequences -output ${sequences}.json -identifier test
    """
}

process compressPredictions {
    publishDir "results", mode: 'copy'

    input:
    path predictions
    path multipleSequenceAlignment
    path tree
    path tree_plot

    output:
    path "*.tar.gz"

    script:
    """
    tar -czvf ${multipleSequenceAlignment.simpleName}.tar.gz $tree $tree_plot $multipleSequenceAlignment $predictions
    """
}

workflow multipleSequenceAlignmentAnalysis {
    take:
    allSequences

    main:
    createMultipleSequenceAlignment(allSequences)
    buildPhilogeneticTree(createMultipleSequenceAlignment.out.multipleSequenceAlignment)
    renderTree(buildPhilogeneticTree.out.tree)

    emit:
    multipleSequenceAlignment = createMultipleSequenceAlignment.out.multipleSequenceAlignment
    tree = buildPhilogeneticTree.out.tree
    tree_plot = renderTree.out.tree_plot
}

workflow b2bToolsAnalysis {
    take:
    sequencesGrouped

    main:
    predictBiophysicalFeatures(sequencesGrouped)

    emit:
    predictions = predictBiophysicalFeatures.out.predictions
}

workflow {
    // createMultipleSequenceAlignment(allSequences)
    // buildPhilogeneticTree(createMultipleSequenceAlignment.out.multipleSequenceAlignment)
    // renderTree(buildPhilogeneticTree.out.tree)
    multipleSequenceAlignmentAnalysis(allSequences)

    // predictBiophysicalFeatures(sequencesGrouped)
    b2bToolsAnalysis(sequencesGrouped)

    compressPredictions(
        b2bToolsAnalysis.out.predictions.collect(),
        multipleSequenceAlignmentAnalysis.out.multipleSequenceAlignment,
        multipleSequenceAlignmentAnalysis.out.tree,
        multipleSequenceAlignmentAnalysis.out.tree_plot
    )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}