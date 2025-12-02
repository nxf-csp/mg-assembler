import os

class Sample:
    
    def __init__(self, sample_id, pool_id):
        self._sample_id = sample_id
        self._pool_id = pool_id
        self._sample_st = None
        self._sample_dir = os.path.join(pool_id, sample_id)
        self._asmb_dir = os.path.join(self._sample_dir, "assembly")
        self._bins_dir = os.path.join(self._sample_dir, "assembly", "bins")
        self._ant_dir = os.path.join(self._sample_dir, "annotations")
        self._rc_dir = os.path.join(self._sample_dir, "readcounts")
        self._reads_dir = os.path.join(self._sample_dir, "reads")
        self._reports_dir = os.path.join(self._sample_dir, "reports")
        
    @property
    def id(self):
        return self._sample_id

    @property
    def pool(self):
        return self._pool_id

    @property
    def dir(self):
        return self._sample_dir

    @property
    def assembly(self):
        return self._asmb_dir

    @property
    def bins(self):
        return self._bins_dir

    @property
    def annots(self):
        return self._ant_dir
    
    @property
    def rcounts(self):
        return self._rc_dir

    @property
    def reads(self):
        return self._reads_dir

    @property
    def reports(self):
        return self._reports_dir

    @property
    def status(self):
        return self._sample_st

    @status.setter
    def status(self, value):
        statuses = ["Assembled", "Assembled (warnings)", "Assembly failed",
                    "Binned", "Binned (warnings)", "Binning failed",
                    "Annotated", "Annotation failed"]
        if value not in statuses:
            raise ValueError("Unexpected sample status specification!")
        self._sample_st = value
    
    @property
    def log_file(self):
        return os.path.join(self._sample_dir, "assembly.log")

    def check_inputs(self):
        if self.__class__.__name__ == "Sample":
            raise AttributeError("Required attributes are not defined for the class 'Sample'! This method is only applicable to inheriting classes.")

        for k in self._inputs.keys():
            if not os.path.exists(self._inputs[k]):
                raise FileNotFoundError(f"Path '{self._inputs[k]}' does not exist!")

    def set_process_status(self, step, st):
        if self.__class__.__name__ == "Sample":
            raise AttributeError("Required attributes are not defined for the class 'Sample'! This method is only applicable to inheriting classes.")
 
        sts = ["COMPLETED", "FAILED", "FINISHED WITH WARNINGS"]
        if step in self._p_status.keys():
            if st in sts:
                self._p_status[step] = st
            else:
                raise ValueError("Unexpected process step status specification!")
        else:
            raise ValueError("Unexpected process step name!")
            
    def memoize_metric(self, metric, value):
        if self.__class__.__name__ == "Sample":
            raise AttributeError("Required attributes are not defined for the class 'Sample'! This method is only applicable to inheriting classes.")
 
        if metric in self._metrics.keys():
            self._metrics[metric] = value
        else:
            raise ValueError("Unexpected name of metric!")

    def formate(self):
        if self.__class__.__name__ == "Sample":
            raise AttributeError("Required attributes are not defined for the class 'Sample'! This method is only applicable to inheriting classes.")

        output_dict = {
            "mgid" : f"{self._pool_id}:{self._sample_id}",
            "sampleStatus" : self._sample_st,
            "processStatus" : self._p_status,
            "outputs" : {k : (v if os.path.exists(v) else None) for k, v in self._outputs.items()},
            "metrics" : self._metrics
        }
        return output_dict

class AssembledSample(Sample):
    
    def __init__(self, sample_id, pool_id, fastqs):
        super().__init__(sample_id, pool_id)
        self._inputs = fastqs    # {"R1" : [r1_1.fq, r1_2.fq...], "R2" : [r2_1.fq, r2_2.fq...]}
        self._swd = os.path.join(self._sample_dir, "mg_assembler")
        self._outputs = {
            "assemblyFasta" : os.path.join(self._asmb_dir, "contigs.fasta"), 
            "assemblyGraph" : os.path.join(self._asmb_dir, "assembly_graph.fastg"), 
            "contigPaths" : os.path.join(self._asmb_dir, "contigs.paths"), 
            "scaffoldPaths" : os.path.join(self._asmb_dir, "scaffolds.paths")
        }
        self._p_status = {
            "inputsCheck" : None,
            "readsProcessing" : None,
            "assembling" : None,
            "assemblyProcessing" : None
        }
        self._metrics = {
            "readsBeforeTrimming" : None,
            "readsAfterTrimming" : None,
            "prcntReadsPassed" : None,
            "contigsNum": None,
            "contigsSumLength" : None,
            "contigMinLength" : None,
            "contigMaxLength" : None,
            "NsNum" : None
        }

    @property
    def fastqs(self):
        return self._inputs
    
    @property  
    def swd(self):
        return self._swd

    @property
    def outputs(self):
        return self._outputs

    def check_inputs(self):
        for fq in [*self._inputs["R1"], *self._inputs["R2"]]:
            if not os.path.exists(fq):
                raise FileNotFoundError(f"Path '{fq}' does not exist!")

class BinnedSample(Sample):
    
    def __init__(self, sample_id, pool_id):
        super().__init__(sample_id, pool_id)
        self._p_status = {
            "inputsCheck" : None,
            "aligning" : None,
            "binning" : None,
            "binsProcessing" : None
        }
        self._inputs = {
            "assemblyFasta" : os.path.join(self._asmb_dir, "contigs.fasta"),
            "r1Paired" : os.path.join(self._reads_dir, "R1_paired.fastq.gz"),
            "r2Paired" : os.path.join(self._reads_dir, "R2_paired.fastq.gz"),
            "r1Unpaired" : os.path.join(self._reads_dir, "R1_unpaired.fastq.gz"),
            "r2Unpaired" : os.path.join(self._reads_dir, "R2_unpaired.fastq.gz")
        }
        self._outputs = {
            "contigToBin" : os.path.join(self._asmb_dir, "contig2bin.tsv"),
            "binsDir" : self._bins_dir
        }
        self._metrics = {
            "nBins" : None
        }
    
    @property
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs


class AnnotatedSample(Sample):
    
    def __init__(self, sample_id, pool_id):
        super().__init__(sample_id, pool_id)
        self._p_status = {
            "inputsCheck" : None,
            "taxonomyIdent&QC" : None,
            "geneCallsProcessing" : None,
            "annotating" : None,
            "readCounting" : None
        }
        self._inputs = {
            "assemblyFasta" : os.path.join(self._asmb_dir, "contigs.fasta"),
            "binsDir" : self._bins_dir,
            "r1Paired" : os.path.join(self._reads_dir, "R1_paired.fastq.gz"),
            "r2Paired" : os.path.join(self._reads_dir, "R2_paired.fastq.gz"),
        }
        #self._ref_db = "",
        self._swd = os.path.join(self._sample_dir, "mg_annotator")
        self._outputs = {
            "binsMetrics" : os.path.join(self._asmb_dir, "bin_metrics.tsv"),
            "annotationGff" : os.path.join(self._ant_dir, "cds_annot.gff"),
            "CDSaaFasta" : os.path.join(self._ant_dir, "cds_aa_seqs.faa"),
            "CDSnucFasta" : os.path.join(self._ant_dir, "cds_nuc_seqs.fna"),
            "CDS2UniRef90" : os.path.join(self._ant_dir, "cds2UniRef90_annot.tsv"),
            "readsPerContig" : os.path.join(self._rc_dir, "reads_per_contig.tsv"),
            "readsPerBin" : os.path.join(self._rc_dir, "reads_per_bin.tsv"),
            "readsPerCDS" : os.path.join(self._rc_dir, "reads_per_cds.tsv")
        }
        self._metrics = {
            "MAGs" : None,
            "CDS" : None,
            "readsTotal" : None,
            "mappedReads" : None,
            "mappedReadsPrcnt" : None,
            "assignedReads" : None,
            "assignedReadsPrcnt" : None
        }
    
    @property
    def inputs(self):
        return self._inputs

    '''
    @property
    def ref_db(self):
        return self._ref_db

    @ref_db.setter
    def ref_db(self, value):
        if not isinstance(value, str) or len(value) == 0:
            raise ValueError("DB path must be a non-empty string.")
        if not os.path.exists(value):
            raise FileNotFoundError(f"Path '{value}' does not exist!")
        self._ref_db = value
    '''
    
    @property  
    def swd(self):
        return self._swd

    @property
    def outputs(self):
        return self._outputs
