import numpy as np
import Enums
class SomMutation:
    """CLASS info

        FUNCTIONS:

            public:

            private:

        PROPERTIES:

            regular variables:

            @property methods: none
    """

    def __init__(self, chrN, pos, ref, alt, ccf_1d, ref_cnt=None, alt_cnt=None, gene=None, prot_change=None,
                 mut_category=None, det_power=None, from_sample=None, type_=None,
                 clonality=None,
                 on_arm_event=None,
                 seg_tree=None,
                 clust_ccf=None,
                 cluster_assignment=None,
                 multiplicity=np.nan,
                 local_cn_a1=np.nan,
                 local_cn_a2=np.nan,
                 type_local_cn_event=np.nan,
                 clonal_cut_off=0.86):

        # try:
        # Link to sample object.
        self.from_sample = from_sample

        # Basic properties of mutation
        self.chrN = chrN
        self.pos = int(float(pos))
        self.ref = ref
        self.alt = alt
        self.ref_cnt = int(
            float(ref_cnt)) if ref_cnt is not None else None  # TODO: add support/handling for text CCF files
        self.alt_cnt = int(float(alt_cnt)) if alt_cnt is not None else None
        self.gene = gene if gene is not None and 'unknown' not in gene.lower() else None

        self.blacklist_status = False  # is blacklisted mutation
        self.graylist_status = False  # is blacklisted mutation

        # raw ccf before clustering
        self.ccf_1d = ccf_1d

        self.status = Enums.MutStatus.OK  # later if blacklisted or graylisted

        self.clean_local_cn(local_cn_a1, local_cn_a2)

        if type_ is None:
            # guess type for backwards compatibility, but prefer that type be passed explicitly.
            if "-" in self.ref:
                self.type = Enums.MutType.INS
            elif "-" in self.alt:
                self.type = Enums.MutType.DEL
            else:
                self.type = Enums.MutType.SNV
        else:
            self.type = type_

        self._var_str = ":".join(map(str, [self.chrN, self.pos, self.ref,
                                           self.alt]))  # variant string, set as private, gotten with self.var_str

        # try:
        #    self.ccf_1d = tuple(map(float, ccf_1d))  # make sure all ccf values are floats
        #    # TODO: CCF broadening - this needs to be implemented. Currently does nothing.
        #    if max(self.ccf_1d) > 0.25:  # hueristic - if 25% of the weight is in one bin, broaden it.
        #        self.ccf_1d = self._broaden_ccf(self.ccf_1d)
        # except:
        #    logging.debug("Cannot parse ccf_1d in SomMutation instantiation")
        #    self.ccf_1d = None

        # self.blacklist_status = False
        # self.graylist_status = False

        """
        TODO:auto set blacklist/graylist status.
        if self.from_sample is not None:
            if self.var_str in self.from_sample.known_blacklisted_mut:
                self.blacklist_status=True
        if self.from_sample
        """

        self.prot_change = prot_change if prot_change is not None and 'unknown' not in prot_change.lower() else None  # Unknown values are set as __UNKNOWN__ in ocotator, but Unknown is also something that should be skipped for gene names.

        # several values/placeholders for output options
        self.mut_category = mut_category
        self.color = None

        # move to output
        # This check is to avoid crashes when det_power is "NA" or a non valid str that still evaluates to true.
        self.det_power = float(det_power) if det_power and det_power not in self.from_sample.na_values else None
        self.power_color = 'white' if not self.det_power else 'rgb({0},{0},{0})'.format(
            int(255 * self.det_power))  # takes floor of the value

        ## TODO: remove this if chain
        # if mut_category:
        #    if 'Missense_Mutation' in mut_category:
        #        self.color = 'blue'
        #    elif 'Silent' in mut_category:
        #        self.color = 'orange'
        #    elif 'Nonsense_Mutation' in mut_category:
        #        self.color = 'red'
        #    elif 'Splice' in mut_category:
        #        self.color = 'purple'
        ##

        # Output of CCF clustering
        self.cluster_assignment = cluster_assignment  # moved to ND VirtualEvent
        # self.post_CCF_1d = None

        # if local_cn_a1 is not None:
        #    self.local_cn_a1 = float(local_cn_a1)
        # else:
        #    self.local_cn_a1 = np.nan
        # if local_cn_a2 is not None:
        #    self.local_cn_a2 = float(local_cn_a2)
        # else:
        #    self.local_cn_a2 = np.nan

    # except:
    #	err = str(sys.exc_info()[1])
    #	raise Exception('\n'.join(['Incorrect input format for variant ', err, 'please review input data!']))  # TODO:Print variant
    #	# This error should no longer come up that often.

    @property
    def var_str(self):
        return self._var_str

    @property
    def ccf_grid_size(self):  # Return the grid size of the CCF.
        return len(self.ccf_1d)

    # TODO make sure things can't be one at time changed from outside class

    @property
    def maf_annotation(self):
        raise NotImplementedError  # not implemented yet

    @classmethod  # a class factory helper to build from different inputs
    def from_list(cls, li, from_sample=None, ccfgrid_size=101):
        ccf_list = li[
                   4:ccfgrid_size + 4]  # skip first 4 element 'Chromosome','Start_position','Reference_Allele','Tumor_Seq_Allele2'
        variant = li[:4]
        in_list = variant + [ccf_list] + li[ccfgrid_size + 4:]
        return cls(*in_list, from_sample=from_sample)

    @classmethod  # a class factory helper to build from different inputs
    def from_dict(cls, required_li, optional_parm_dict, from_sample=None, ccfgrid_size=101):
        optional_parm_dict['from_sample'] = from_sample
        return cls(*required_li, **optional_parm_dict)

    @classmethod  # a class factory helper to build from different inputs
    def from_som_mutation_zero(cls, som_mut, from_sample=None, ccfgrid_size=101):
        ccf_1d = ['1.00'] + ['0.00'] * (ccfgrid_size - 1)  # assume ccf=0
        return cls(som_mut.chrN, som_mut.pos, som_mut.ref, som_mut.alt, ccf_1d, ref_cnt=som_mut.ref_cnt, alt_cnt=0,
                   gene=som_mut.gene, prot_change=som_mut.prot_change, mut_category=som_mut.mut_category,
                   from_sample=from_sample)

    # a method to define == operator to compare mutations , i.e. Mut1==Mut2 if their variant strings equal
    # TODO be careful to avoid bugs when a specifc sample is required (e.g. alt/ref data)
    def __eq__(self, other):
        if other is None:
            return False
        if self.var_str == other.var_str:
            return True
        else:
            return False

    # hash method to define mutation equality in sets, dictionaries etc.
    def __hash__(self):
        return hash(self.var_str)

    def __str__(self):
        return self.var_str

    def __repr__(self):
        # print "str_called"
        return self.var_str

    def clean_local_cn(self, cn1, cn2):
        self.local_cn_a1 = float(cn1) if not np.nan else np.nan
        self.local_cn_a2 = float(cn2) if not np.nan else np.nan

    # update allelic copy number AND assignment of mutation to a copy number event (arm level)
    def _phase_mutation(self, bam_file):
        # TODO: Phase mutation to an allele based on shifts of local hetsites and local copy number
        raise NotImplementedError

    # @classmethod #a class factory helper to build from different inputs
    def _from_clustered_maf(cls):
        raise NotImplementedError

    # @classmethod #a class factory helper to build from different inputs
    def _from_absolute_RData(cls):
        # TODO: implement function

        raise NotImplementedError