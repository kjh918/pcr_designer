import os
import primer3
from primer.pcr_components import Primer, Amplicon
from Bio.Seq import reverse_complement
from primer.designer import PrimerDesigner


class qPCRdesigner():
    """
    """

    #TODO input as sequence
    def __init__(self, f_reference_fasta, chrom, start, end, region_id=None, min_amplicon_length=80, max_amplicon_length=120, 
                n_probes=100, n_primers=100, min_primer_probe_tm_diff=6, max_primer_probe_tm_diff=8, 
                bisulfite=False, cpg_default='methyl', methylation_pattern=[], probe_kwargs={}, primer_kwargs={}):
        
        self.f_reference_fasta = f_reference_fasta
        self.chrom = chrom
        self.start = start
        self.end = end
        if region_id == None:
            self.region_id = f'{chrom}:{start}-{end}'
        else:
            self.region_id = region_id
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        self.bisulfite = bisulfite
        self.reference_template_sequence = f_reference_fasta.fetch(chrom, end-max_amplicon_length, start+max_amplicon_length+1).upper()
        if bisulfite == True:
            self.template_sequence = self.bisulfite_conversion(f_reference_fasta, chrom, end-max_amplicon_length, start+max_amplicon_length+1, cpg_default=cpg_default, methylation_pattern=methylation_pattern)
        else:
            self.template_sequence = f_reference_fasta.fetch(chrom, end-max_amplicon_length, start+max_amplicon_length+1).upper()

        self.target_start_index = max_amplicon_length-(end-start)-1
        self.target_end_index = max_amplicon_length-2
        self.target_sequence = self.template_sequence[self.target_start_index:self.target_end_index+1]

        self.n_probes = n_probes
        self.n_primers = n_primers
        self.min_primer_probe_tm_diff = min_primer_probe_tm_diff
        self.max_primer_probe_tm_diff = max_primer_probe_tm_diff
        self.probe_kwargs = probe_kwargs
        self.primer_kwargs = primer_kwargs

        self.probe_designer = None
        self.primer_designer = None
        self.amplicon_list = []

        #### ####
        self.design_probe()
        self.design_primer()

    @staticmethod
    def bisulfite_conversion(f_reference_fasta, chrom, start, end, cpg_default='methyl', methylation_pattern=[]):
        #TODO apply to reverse strand
        print('??????????')
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
        self.probe_kwargs['forward_primer'] = False
        self.probe_kwargs['reverse_primer'] = False
        self.probe_kwargs['probe'] = True
        self.probe_kwargs['target_start_index'] = self.target_start_index
        self.probe_kwargs['target_end_index'] = self.target_end_index
        self.probe_kwargs['min_amplicon_length'] = self.min_amplicon_length
        self.probe_kwargs['max_amplicon_length'] = self.max_amplicon_length
        self.probe_kwargs['n_primers'] = self.n_probes
        print(self.probe_designer)

        probe_designer = PrimerDesigner(**self.probe_kwargs)
        probe_designer.design_primer()
        print(123123)
        self.probe_designer = probe_designer

    def design_primer(self):
        """
        """
        if self.probe_designer == None:
            self.design_probe()
            print(1)
        else:
            pass
        print('okoko')
        print(self.primer_kwargs)
        print(self.probe_designer)

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
            
            print(self.primer_kwargs)


            primer_designer = PrimerDesigner(**self.primer_kwargs)
            primer_designer.design_primer()
            self.primer_designer = primer_designer
            self.amplicon_list += self.primer_designer.amplicon_list

