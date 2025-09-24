
import logging
import TumorSample
from intervaltree import Interval, IntervalTree
class Patient:
    """CLASS info

     FUNCTIONS:
         public:

           add_sample() - add new samples from individual to self.sample_list by making a new TumorSample instance

        private:

           _auto_file_type() - guess input file type for ccf data if not manually provided

     PROPERTIES:
         regular variables:
            self.sample_list

         @property methods: none
    """

    def __init__(self, indiv_name='Indiv1',
                 ccf_grid_size=101,
                 impute_missing=False,
                 artifact_whitelist='',
                 use_indels=False,
                 min_coverage=8,
                 PoN_file=False):

        # DECLARATIONS
        self.indiv_name = indiv_name
        # :type : list [TumorSample]
        self.sample_list = []
        self.rna_sample_list = []

        self.samples_synchronized = False
        self.rna_samples_synchronized = False
        self.concordant_genes = []
        self.ccf_grid_size = ccf_grid_size

        self.PatientLevel_MutWhitelist = artifact_whitelist

        # Patient configuration settings
        # flag if to impute missing variants as ccf 0
        self.impute_missing = impute_missing
        # min cov is specified both here and passed to tumor sample.
        self.min_coverage = min_coverage
        self.use_indels = use_indels
        self.PoN_file = PoN_file

        # later filled data objects
        self.ND_mutations = []

        # storing of results
        # Clustering
        self.ClusteringResults = None
        self.MutClusters = None
        self.TruncalMutEvents = None
        self.MCMC_trace = None
        self.k_trace = None
        self.alpha_trace = None

        self.unclustered_muts = []

        self.concordant_cn_tree = {chrom: IntervalTree() for chrom in list(map(str, range(1, 23))) + ['X', 'Y']}

        # BuildTree
        self.TopTree = None
        self.TreeEnsemble = []
        self.alt_cn_states = []

    def initPatient(self):
        """ accepted input types abs; txt; sqlite3 .db # auto tab if .txt, .tsv or .tab ; abs if .Rdata; sqlite if .db """
        raise NotImplementedError

    def addSample(self, filen, sample_name, input_type='auto', seg_input_type='auto', grid_size=101, seg_file=None,
                  _additional_muts=None, purity=None, timepoint_value=None):
        """ accepted input types abs; txt; sqlite3 .db # auto tab if .txt, .tsv or .tab ; abs if .Rdata; sqlite if .db"""

        if _additional_muts == []:
            _additional_muts = None
        elif type(_additional_muts) is list:
            _additional_muts = _additional_muts[0]

        # make new sample and add to exiting list of samples
        logging.info("Adding Mutations from Sample: %s", sample_name)
        new_sample = TumorSample.TumorSample(filen, input_type, seg_input_type=seg_input_type, sample_name=sample_name,
                                 artifact_whitelist=self.PatientLevel_MutWhitelist,
                                 ccf_grid_size=grid_size, PoN=self.PoN_file, indiv=self.indiv_name,
                                 use_indels=self.use_indels, min_coverage=self.min_coverage,
                                 _additional_muts=_additional_muts, seg_file=seg_file,
                                 purity=purity, timepoint_value=timepoint_value)

        self.sample_list.append(new_sample)
        logging.info('Added sample ' + new_sample.sample_name)
        # turn of concordance flag when new sample is added
        self.samples_synchronized = False

        return True