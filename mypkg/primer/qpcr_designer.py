import os
import primer3
from primer.pcr_components import Primer, Amplicon
from Bio.Seq import reverse_complement
from primer.designer import PrimerDesigner
import pysam


class qPCRdesigner():
    """
    """

    #TODO input as sequence
    def __init__(self, f_reference_fasta, chrom, start, end, region_id=None, 
                min_amplicon_length=80, max_amplicon_length=120, 
                n_probes=100, n_primers=100, bisulfite=False, cpg_default='methyl', 
                methylation_pattern=[], probe_kwargs={}, primer_kwargs={}):
        
        self.f_reference_fasta = f_reference_fasta
        self.chrom = chrom
        self.start = start - 1   
        self.end = end           

        if region_id is None:
            self.region_id = f"{chrom}:{start}-{end}"
        else:
            self.region_id = region_id

        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        self.bisulfite = bisulfite

        # template fetch 구간
        template_start = end - max_amplicon_length
        template_end   = self.start + max_amplicon_length + 1

        self.reference_template_sequence = f_reference_fasta.fetch(
            chrom, template_start, template_end
        ).upper()

        if bisulfite:
            self.template_sequence = self.bisulfite_conversion(
                f_reference_fasta, chrom, template_start, template_end,
                cpg_default=cpg_default, methylation_pattern=methylation_pattern
            )
        else:
            self.template_sequence = self.reference_template_sequence

        # ✔ target index 계산 (template_sequence 기준)
        self.target_start_index = self.start - template_start
        self.target_end_index   = self.end - template_start

        # ✔ target_sequence 추출
        self.target_sequence = self.template_sequence[
            self.target_start_index : self.target_end_index
        ]

        self.n_probes = n_probes
        self.n_primers = n_primers
        self.probe_kwargs = probe_kwargs
        self.primer_kwargs = primer_kwargs
        self.probe_designer = None
        self.primer_designer = None
        self.amplicon_list = []

        if self.n_probes == 0:
            self.design_primer_only()
        #### ####
        else:
            self.design_probe()
            self.design_primer()

    @staticmethod
    def bisulfite_conversion(f_reference_fasta, chrom, start, end, cpg_default='methyl', methylation_pattern=[]):
        #TODO apply to reverse strand
        sequence = f_reference_fasta.fetch(chrom, start, end+1).upper()
        converted_sequence = ''
        for i, base in enumerate(sequence):
            if i == len(sequence)-1:
                break
            base_genomic_position = start+i
            if sequence[i] == 'C' and sequence[i+1] == 'G':
                if base_genomic_position in methylation_pattern:
                    converted_sequence += 'C'
                elif -base_genomic_position in methylation_pattern:
                    converted_sequence += 'T'
                else:
                    if cpg_default == 'methyl':
                        converted_sequence += 'C'
                    elif cpg_default == 'unmethyl':
                        converted_sequence += 'T'
            elif sequence[i] == 'C' and sequence[i+1] != 'G':
                converted_sequence += 'T'
            else:
                converted_sequence += base
        
        return converted_sequence

    def design_probe(self):
        """
        """
        self.probe_kwargs['template_sequence'] = self.template_sequence
        self.probe_kwargs['reference_template_sequence'] = self.reference_template_sequence
        self.probe_kwargs['forward_primer'] = True
        self.probe_kwargs['reverse_primer'] = True
        self.probe_kwargs['probe'] = False
        self.probe_kwargs['target_start_index'] = self.target_start_index
        self.probe_kwargs['target_end_index'] = self.target_end_index
        self.probe_kwargs['min_amplicon_length'] = self.min_amplicon_length
        self.probe_kwargs['max_amplicon_length'] = self.max_amplicon_length
        self.probe_kwargs['n_primers'] = self.n_probes

        probe_designer = PrimerDesigner(**self.probe_kwargs)
        probe_designer.design_primer()
        self.probe_designer = probe_designer
        

    def design_primer_only(self):
        """
        """
        self.primer_kwargs['template_sequence'] = self.template_sequence
        self.primer_kwargs['reference_template_sequence'] = self.reference_template_sequence
        self.primer_kwargs['forward_primer'] = True
        self.primer_kwargs['reverse_primer'] = True
        self.primer_kwargs['probe'] = False
        self.primer_kwargs['target_start_index'] = self.target_start_index
        self.primer_kwargs['target_end_index'] = self.target_end_index
        self.primer_kwargs['min_amplicon_length'] = self.min_amplicon_length
        self.primer_kwargs['max_amplicon_length'] = self.max_amplicon_length
        self.primer_kwargs['n_primers'] = self.n_primers

        primer_designer = PrimerDesigner(**self.primer_kwargs)
        primer_designer.design_primer()
        self.primer_designer = primer_designer
        self.amplicon_list += self.primer_designer.amplicon_list
        

    def design_primer(self):
        """
        """
        if self.probe_designer == None:
            self.design_probe()
        else:
            pass
        
        for probe_rank, probe_amplicon in enumerate(self.probe_designer.amplicon_list):
            self.primer_kwargs['template_sequence'] = self.template_sequence
            self.primer_kwargs['reference_template_sequence'] = self.reference_template_sequence
            self.primer_kwargs['forward_primer'] = True
            self.primer_kwargs['reverse_primer'] = True
            self.primer_kwargs['probe'] = True
            self.primer_kwargs['target_start_index'] = self.target_start_index
            self.primer_kwargs['target_end_index'] = self.target_end_index
            self.primer_kwargs['opt_tm'] = probe_amplicon.probe.tm-self.min_primer_probe_tm_diff
            self.primer_kwargs['min_tm'] = probe_amplicon.probe.tm-self.max_primer_probe_tm_diff
            self.primer_kwargs['max_tm'] = probe_amplicon.probe.tm-self.min_primer_probe_tm_diff
            self.primer_kwargs['min_amplicon_length'] = self.min_amplicon_length
            self.primer_kwargs['max_amplicon_length'] = self.max_amplicon_length
            self.primer_kwargs['n_primers'] = self.n_primers
            probe_sequence = probe_amplicon.probe.sequence
            self.primer_kwargs['probe_sequence'] = probe_sequence
            

            primer_designer = PrimerDesigner(**self.primer_kwargs)
            primer_designer.design_primer()
            self.primer_designer = primer_designer
            self.amplicon_list += self.primer_designer.amplicon_list

