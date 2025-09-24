# Main run method
import logging
def run_tool(args):
    logging.debug('Arguments {}'.format(args))

    from PatientObject import Patient
    from TreeObject import Tree
    from CellPopulationEngine import CellPopulationEngine
    from BuildTreeEngine import BuildTreeEngine


    # init a Patient
    patient_data = Patient(indiv_name=args.indiv_id)
    # try:  # if sif file is specified
    # Patient load cluster and mut ccf files
    parse_sif_file(args.sif, args.mutation_ccf_file, patient_data)
    load_clustering_results(args.cluster_ccf_file, patient_data)

    # Select the tree number
    tree_edges = load_tree_edges_file(args.tree_tsv, tree_num=args.tree_number)
    bt_engine = BuildTreeEngine(patient_data)
    tree = Tree()
    tree.init_tree_from_clustering(patient_data.ClusteringResults)
    tree.set_new_edges(tree_edges)
    patient_data.TopTree = tree
    bt_engine.set_top_tree(tree)
    # Computing Cell Population
    cp_engine = CellPopulationEngine(patient_data, seed=args.seed)
    constrained_ccf = cp_engine.compute_constrained_ccf()

    cell_ancestry = bt_engine.get_cell_ancestry()
    cell_abundance = cp_engine.get_cell_abundance_across_samples(constrained_ccf)

    # Output and visualization
    import PhylogicOutput
    phylogicoutput = PhylogicOutput.PhylogicOutput()
    # TODO write cell population MCMC trace to file
    phylogicoutput.write_all_cell_abundances(cp_engine.get_all_cell_abundances(), args.indiv_id)
    phylogicoutput.write_constrained_ccf_tsv(constrained_ccf, cell_ancestry, args.indiv_id)
    phylogicoutput.write_cell_abundances_tsv(cell_abundance, cell_ancestry, args.indiv_id)
    # phylogicoutput.generate_html_from_tree(args.mutation_ccf_file, args.cluster_ccf_file,
    #                                        args.indiv_id + '_build_tree_posteriors.tsv',
    #                                        args.indiv_id + '_constrained_ccf.tsv',
    #                                        sif=args.sif,
    #                                        drivers=patient_data.driver_genes,
    #                                        treatment_file=args.treatment_data,
    #                                        tumor_sizes_file=args.tumor_size,
    #                                        cnv_file=args.indiv_id + '.cnvs.txt')

def load_tree_edges_file(tree_tsv, tree_num):
    """Needs to return a list of tuples in the form of (parent, child)
    Example:
    [(None, 1),
     (1, 2),
     (1, 5),
     (2, 3),
     (2, 7),
     (3, 4),
     (3, 6),
     (5, 10),
     (6, 8),
     (6, 9)]
    """
    reader = open(tree_tsv, 'r')
    header = reader.readline()
    for i in range(0, tree_num):
        top_tree = reader.readline().split('\t')[-1].strip()
    edges = [tuple([None if x == 'None' else int(x) for x in item.split('-')]) for item in top_tree.split(',')]
    edges = [edges[-1]] + edges[:-1]
    return edges

def parse_sif_file(sif_file, mutation_ccf_file, patient_data):
    # TODO: duplicate code in Cluster
    with open(sif_file, 'r') as reader:
        for line in reader:
            if not line.strip() == "":
                # for now, assume input file order is of the type sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint
                values = line.strip('\n').split('\t')
                if line.startswith('sample_id'):
                    header = line.strip('\n').split('\t')
                    header = {x: i for i, x in enumerate(header)}
                else:
                    sample_id = values[header['sample_id']]
                    seg_fn = values[header['seg_fn']]
                    purity = float(values[header['purity']])
                    timepoint = float(values[header['timepoint']])
                    logging.debug("Adding sample {}".format(sample_id))

                    patient_data.addSample(mutation_ccf_file, sample_id,
                                           input_type="post-clustering",
                                           timepoint_value=timepoint,
                                           seg_file=seg_fn,
                                           purity=purity)

def load_clustering_results(cluster_info_file, patient_data):
    import ClusterObject
    clustering_results = {}
    ccf_headers = ['postDP_ccf_' + str(i / 100.0) for i in range(0, 101, 1)]
    sample_names = [sample.sample_name for sample in patient_data.sample_list]
    with open(cluster_info_file, 'r') as reader:
        for line in reader:
            values = line.strip('\n').split('\t')
            if line.startswith('Patient_ID'):
                header = dict((item, idx) for idx, item in enumerate(values))
            else:
                sample_id = values[header['Sample_ID']]
                cluster_id = int(values[header['Cluster_ID']])
                ccf = [float(values[header[i]]) for i in ccf_headers]
                if cluster_id not in clustering_results:
                    new_cluster = ClusterObject.Cluster(cluster_id, sample_names)
                    clustering_results[cluster_id] = new_cluster
                    logging.debug('Added cluster {} '.format(cluster_id))
                clustering_results[cluster_id].add_sample_density(sample_id, ccf)
    for cluster_id, cluster in clustering_results.items():
        cluster.set_blacklist_status()
        clustering_results[cluster_id] = cluster

    mutations_nd_hist = {}
    for sample in patient_data.sample_list:
        for mutation in sample.mutations:
            if mutation not in mutations_nd_hist:
                mutations_nd_hist[mutation] = []
            mutations_nd_hist[mutation].append(mutation.ccf_1d)
    for mutation, mutation_nd_hist in mutations_nd_hist.items():
        clustering_results[mutation.cluster_assignment].add_mutation(mutation, mutation_nd_hist)
    patient_data.ClusteringResults = clustering_results