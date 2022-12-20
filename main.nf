#!/usr/bin/env nextflow
// USAGE: nextflow run main.nf -resume -with-dag pipeline.png

params.targetSequences = "$launchDir/input_small.fasta"
params.groupBy = 1000

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

process buildLogo {
    input:
    path multipleSequenceAlignment

    output:
    path "*_logo.png", emit: logo

    script:
    """
    weblogo --sequence-type protein --title "MSA logo" --size large --format png_print < $multipleSequenceAlignment > ${multipleSequenceAlignment.simpleName}_logo.png
    """
}


process renderTree {
    input:
    path tree

    output:
    path "tree.png", emit: treePlot

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
    #!/usr/local/bin/python
    from b2bTools import SingleSeq
    import json

    single_seq = SingleSeq("$sequences")
    single_seq.predict(tools=['dynamine', 'efoldmine', 'disomine', 'agmata', 'psp'])

    all_predictions = single_seq.get_all_predictions()
    json.dump(all_predictions, open('b2b_results_${sequences.baseName}.json', 'w'), indent=4, sort_keys=True)
    """
}

process plotBiophysicalFeatures {
    tag "${predictions.baseName}"

    input:
    path predictions

    output:
    path "*_predictions.png", emit: plots

    script:
    """
    #!/usr/bin/python3

    import matplotlib.pyplot as plt
    import json

    with open('$predictions', 'r') as json_file:
        prediction_dict = json.loads(json_file.read())

    for id, prediction in enumerate(prediction_dict['results']):
        fig, axs = plt.subplots(2, 4)
        ax1 = axs[0, 0]
        ax2 = axs[0, 1]
        ax3 = axs[0, 2]
        ax4 = axs[0, 3]
        ax5 = axs[1, 0]
        ax6 = axs[1, 1]
        ax7 = axs[1, 2]
        ax8 = axs[1, 3]

        fig.set_figwidth(30)
        fig.set_figheight(10)

        fig.suptitle(f"Single Sequence Predictions for: {prediction['proteinID']}")

        x_position = range(len(prediction['sequence']))
        backbone_pred = prediction['backbone']
        coil_pred = prediction['coil']
        sheet_pred = prediction['sheet']
        ppII_pred = prediction['ppII']
        helix_pred = prediction['helix']
        sidechain_pred = prediction['sidechain']
        disomine_pred = prediction['disoMine']
        earlyFolding_pred = prediction['earlyFolding']

        ax1.plot(x_position, backbone_pred, label="Backbone")
        ax2.plot(x_position, sidechain_pred, label="Side chain")
        ax3.plot(x_position, coil_pred, label="Coil")
        ax4.plot(x_position, sheet_pred, label="Sheet")
        ax5.plot(x_position, ppII_pred, label="ppII")
        ax6.plot(x_position, helix_pred, label="Helix")
        ax7.plot(x_position, disomine_pred, label="Disorder")
        ax8.plot(x_position, earlyFolding_pred, label="Early folding")

        ax1.set_title('DynaMine backbone dynamics')
        ax1.set_ylim([-0.2, 1.2])
        ax1.set_xlabel('residue index')
        ax1.set_ylabel('prediction values')
        ax1.axhspan(1, 1.2, alpha=0.3, color='red')
        ax1.axhspan(0.8, 1, alpha=0.5, color='pink')
        ax1.axhspan(0.69, 0.8, alpha=0.5, color='orange')
        ax1.axhspan(-0.2, 0.69, alpha=0.5, color='yellow')
        ax1.grid(axis='y')
        ax1.set_xlim([0, len(prediction['sequence']) - 1])

        ax2.set_title('DynaMine sidechain dynamics')
        ax2.set_ylim([-0.2, 1.2])
        ax2.set_xlabel('residue index')
        ax2.set_ylabel('prediction values')
        ax2.grid(axis='y')
        ax2.set_xlim([0, len(prediction['sequence']) - 1])

        ax3.set_title('DynaMine conformational propensities: Coil')
        ax3.set_ylim([-0.2, 1.2])
        ax3.set_xlabel('residue index')
        ax3.set_ylabel('prediction values')
        ax3.grid(axis='y')
        ax3.set_xlim([0, len(prediction['sequence']) - 1])

        ax4.set_title('DynaMine conformational propensities: Sheet')
        ax4.set_ylim([-0.2, 1.2])
        ax4.set_xlabel('residue index')
        ax4.set_ylabel('prediction values')
        ax4.grid(axis='y')
        ax4.set_xlim([0, len(prediction['sequence']) - 1])

        ax5.set_title('DynaMine conformational propensities: ppII (polyproline II)')
        ax5.set_ylim([-0.2, 1.2])
        ax5.set_xlabel('residue index')
        ax5.set_ylabel('prediction values')
        ax5.grid(axis='y')
        ax5.set_xlim([0, len(prediction['sequence']) - 1])

        ax6.set_title('DynaMine conformational propensities: Helix')
        ax6.set_ylim([-0.2, 1.2])
        ax6.set_xlabel('residue index')
        ax6.set_ylabel('prediction values')
        ax6.grid(axis='y')
        ax6.set_xlim([0, len(prediction['sequence']) - 1])

        ax7.set_title('Early folding (EFoldMine)')
        ax7.set_ylim([-0.2, 1.2])
        ax7.set_xlabel('residue index')
        ax7.set_ylabel('prediction values')
        ax7.axhspan(-0.2, 0.169, alpha=0.5, color='yellow')
        ax7.axhspan(0.169, 1.2, alpha=0.5, color='orange')
        ax7.grid(axis='y')
        ax7.set_xlim([0, len(prediction['sequence']) - 1])

        ax8.set_title('Disorder (disoMine)')
        ax8.set_ylim([-0.2, 1.2])
        ax8.set_xlabel('residue index')
        ax8.set_ylabel('prediction values')
        ax8.axhspan(0.5, 1.2, alpha=0.5, color='orange')
        ax8.axhspan(-0.2, 0.5, alpha=0.5, color='yellow')
        ax8.grid(axis='y')
        ax8.set_xlim([0, len(prediction['sequence']) - 1])

        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, shadow=True, ncol=8)
        plt.subplots_adjust(hspace=0.4)
        plt.savefig(prediction['proteinID'] + '_predictions.png')
    """
}

process plotAgmata {
    tag "${predictions.baseName}"

    input:
    path predictions

    output:
    path "*_agmata_prediction.png", emit: plots

    script:
    """
    #!/usr/bin/python3

    import matplotlib.pyplot as plt
    import json

    with open('$predictions', 'r') as json_file:
        prediction_dict = json.loads(json_file.read())

    for id, prediction in enumerate(prediction_dict['results']):
        fig, ax = plt.subplots(1, 1)
        fig.set_figwidth(30)
        fig.set_figheight(5)
        fig.suptitle('Agmata aggregation propensity')

        agmata_pred = prediction['agmata']
        ax.plot(range(len(agmata_pred)), agmata_pred, label="AgMata")

        ax.set_xlim([0, len(agmata_pred) - 1])
        ax.set_xlabel('residue index')
        ax.set_ylabel('prediction values')
        ax.grid(axis='y')

        plt.savefig(prediction['proteinID'] + '_agmata_prediction.png')
    """
}

process fetchStructure {
    tag "$id"

    input:
    tuple val(id), val(sequence)

    output:
    path "*.pdb", emit: pdbStructures

    script:
    """
    curl -X POST --data "$sequence" https://api.esmatlas.com/foldSequence/v1/pdb/ > ${id}.pdb
    """
}

process compressPredictions {
    publishDir "results", mode: 'copy'

    input:
    path predictions
    // path plots
    // path agmata_plots

    // path multipleSequenceAlignment
    // path tree
    // path treePlot
    // path logo

    // path pdbStructures

    output:
    path "*.tar.gz"

    script:
    """
    tar -czvhf b2b_results.tar.gz $predictions
    """
}

workflow multipleSequenceAlignmentAnalysis {
    take:
    allSequences

    main:
    createMultipleSequenceAlignment(allSequences)
    buildPhilogeneticTree(createMultipleSequenceAlignment.out.multipleSequenceAlignment)
    renderTree(buildPhilogeneticTree.out.tree)

    buildLogo(createMultipleSequenceAlignment.out.multipleSequenceAlignment)

    emit:
    multipleSequenceAlignment = createMultipleSequenceAlignment.out.multipleSequenceAlignment
    tree = buildPhilogeneticTree.out.tree
    treePlot = renderTree.out.treePlot
    logo = buildLogo.out.logo
}

workflow b2bToolsAnalysis {
    take:
    sequencesGrouped

    main:
    predictBiophysicalFeatures(sequencesGrouped)
    // plotBiophysicalFeatures(predictBiophysicalFeatures.out.predictions)
    // plotAgmata(predictBiophysicalFeatures.out.predictions)

    emit:
    predictions = predictBiophysicalFeatures.out.predictions
    // plots = plotBiophysicalFeatures.out.plots
    // agmata_plots = plotAgmata.out.plots
}

workflow {
    // First sub-workflow
    // multipleSequenceAlignmentAnalysis(allSequences)
    // Second sub-workflow
    b2bToolsAnalysis(sequencesGrouped)

    // Third sub-workflow
    // fetchStructure(allSequences.splitFasta(record: [id: true, seqString: true ]).filter { record -> record.seqString.length() < 400 })

    // Main workflow
    compressPredictions(
        b2bToolsAnalysis.out.predictions.collect(),
        // b2bToolsAnalysis.out.plots.collect(),
        // b2bToolsAnalysis.out.agmata_plots.collect(),
        // multipleSequenceAlignmentAnalysis.out.multipleSequenceAlignment,
        // multipleSequenceAlignmentAnalysis.out.tree,
        // multipleSequenceAlignmentAnalysis.out.treePlot,
        // multipleSequenceAlignmentAnalysis.out.logo,
        // fetchStructure.out.pdbStructures.collect()
    )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
