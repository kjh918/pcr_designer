import os
import primer3
from primer.pcr_components import Primer, Amplicon
from Bio.Seq import reverse_complement

class PrimerDesigner():
    """
    """

    DEFAULT_SALT_MONOVALENT = 50,
    DEFAULT_SALT_DIVALENT = 1.5,
    DEFAULT_DNTP_CONC = 0.6,
    DEFAULT_DNA_CONC = 50,

    #TODO design for reverse strand probe
    def __init__(self, template_sequence, target_start_index, target_end_index,
                min_amplicon_length=80, max_amplicon_length=100, n_primers=10, max_tm_difference=2,
                forward_primer=True, reverse_primer=True, probe=False, probe_sequence=None,
                opt_length=25, min_length=20, max_length=30,
                opt_tm=60, min_tm=50, max_tm=70,
                opt_gc=45, min_gc=35, max_gc=65,
                reference_template_sequence=None,
                primer3_seq_args=None, primer3_global_args=None) -> None:

        if reference_template_sequence != None:
            self.reference_template_sequence = reference_template_sequence
        else:
            self.reference_template_sequence = template_sequence
        self.template_sequence = template_sequence
        self.target_start_index = target_start_index
        self.target_end_index = target_end_index
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        self.n_primers = n_primers
        self.max_tm_difference = max_tm_difference
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.probe = probe
        self.probe_sequence = probe_sequence
        self.opt_tm = opt_tm
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.opt_length = opt_length
        self.min_length = min_length
        self.max_length = max_length
        self.opt_gc = opt_gc
        self.min_gc = min_gc
        self.max_gc = max_gc

        self.primer3_seq_args = None
        self.primer3_global_args = None
        
        self.primer3_result = None
        self.amplicon_list = []
    
        self.primer3_seq_args = {
            'SEQUENCE_ID': 'PRIMER',
            'SEQUENCE_TEMPLATE': self.template_sequence,
        }
        self.primer3_global_args = {
            'PRIMER_TASK': 'generic',
            'PRIMER_NUM_RETURN': n_primers,
            'PRIMER_PICK_LEFT_PRIMER': int(forward_primer),
            'PRIMER_PICK_RIGHT_PRIMER': int(reverse_primer),
            'PRIMER_PICK_INTERNAL_OLIGO': int(probe),
        }

        if probe == True:
            self.update_primer3_seq_args(
                {
                    'SEQUENCE_TARGET': [target_start_index, target_end_index-target_start_index+1],
                    'SEQUENCE_INTERNAL_EXCLUDED_REGION': [[0, target_end_index-self.max_length], [target_start_index+self.max_length, len(self.template_sequence)-(target_start_index+self.max_length)]]
                }
            )
            self.update_primer3_global_args(
                {
                    'PRIMER_INTERNAL_SALT_MONOVALENT': 50,
                    'PRIMER_INTERNAL_SALT_DIVALENT': 1.5,
                    'PRIMER_INTERNAL_DNTP_CONC': 0.6,
                    'PRIMER_INTERNAL_DNA_CONC': 50,
                    'PRIMER_INTERNAL_OPT_SIZE': opt_length,
                    'PRIMER_INTERNAL_MIN_SIZE': min_length,
                    'PRIMER_INTERNAL_MAX_SIZE': max_length,
                    'PRIMER_INTERNAL_OPT_TM': opt_tm,
                    'PRIMER_INTERNAL_MIN_TM': min_tm,
                    'PRIMER_INTERNAL_MAX_TM': max_tm,
                    'PRIMER_INTERNAL_OPT_GC_PERCENT': opt_gc,
                    'PRIMER_INTERNAL_MIN_GC': min_gc,
                    'PRIMER_INTERNAL_MAX_GC': max_gc,
                }
            )
        
        if forward_primer == True or reverse_primer == True:
            self.update_primer3_seq_args(
                {
                    'SEQUENCE_TARGET': [target_start_index, target_end_index-target_start_index+1]
                }
            )
            self.update_primer3_global_args(
                {   
                    'PRIMER_PAIR_MAX_DIFF_TM': max_tm_difference,
                    'PRIMER_OPT_SIZE': opt_length,
                    'PRIMER_MIN_SIZE': min_length,
                    'PRIMER_MAX_SIZE': max_length,
                    'PRIMER_OPT_TM': opt_tm,
                    'PRIMER_MIN_TM': min_tm,
                    'PRIMER_MAX_TM': max_tm,
                    'PRIMER_OPT_GC_PERCENT': opt_gc,
                    'PRIMER_MIN_GC': min_gc,
                    'PRIMER_MAX_GC': max_gc,
                    'PRIMER_PRODUCT_SIZE_RANGE': [min_amplicon_length, max_amplicon_length]
                }
            )

        if probe_sequence != None:
            probe_start = self.template_sequence.find(probe_sequence)
            self.update_primer3_seq_args(
                {
                    'SEQUENCE_INTERNAL_OLIGO': probe_sequence,
                    'SEQUENCE_EXCLUDED_REGION': [[probe_start-1, len(probe_sequence)+2]]#TODO move spacing between probe and primer to argument
                }
            )
            self.update_primer3_global_args(
                {
                    'PRIMER_PICK_INTERNAL_OLIGO': int(True),
                    'PRIMER_INTERNAL_SALT_MONOVALENT': 50,
                    'PRIMER_INTERNAL_SALT_DIVALENT': 1.5,
                    'PRIMER_INTERNAL_DNTP_CONC': 0.6,
                    'PRIMER_INTERNAL_DNA_CONC': 50,
                    'PRIMER_INTERNAL_MIN_SIZE': 0,
                    'PRIMER_INTERNAL_MAX_SIZE': 30,
                    'PRIMER_INTERNAL_MIN_TM': 0,
                    'PRIMER_INTERNAL_MAX_TM': 100,
                    'PRIMER_INTERNAL_MIN_GC': 0,
                    'PRIMER_INTERNAL_MAX_GC': 100,
                }
            )

        if primer3_seq_args != None:
            self.update_primer3_seq_args(primer3_seq_args)
        if primer3_global_args != None:
            self.update_primer3_global_args(primer3_global_args)

    def update_primer3_seq_args(self, primer3_seq_args):
        for key, value in primer3_seq_args.items():
            self.primer3_seq_args[key] = value
    
    def update_primer3_global_args(self, primer3_global_args):
        for key, value in primer3_global_args.items():
            self.primer3_global_args[key] = value

    def run_primer3(self):
        
        self.primer3_result = primer3.bindings.designPrimers(
            seq_args = self.primer3_seq_args,
            global_args = self.primer3_global_args
        )
        n_forward_primers = self.primer3_result['PRIMER_LEFT_NUM_RETURNED']
        n_reverse_primers = self.primer3_result['PRIMER_RIGHT_NUM_RETURNED']
        n_probes = self.primer3_result['PRIMER_INTERNAL_NUM_RETURNED']
        n_primer_pairs = self.primer3_result['PRIMER_PAIR_NUM_RETURNED']

        n_designed_primer = max(n_forward_primers, n_reverse_primers, n_probes, n_primer_pairs)
        amplicon_list = []
        forward_primer = None
        reverse_primer = None
        probe = None
        for primer3_rank in range(0, n_designed_primer):
            if self.probe == True:
                print(1)
                
                if self.primer3_result.get(f'PRIMER_LEFT_{primer3_rank}') != None:
                    forward_primer = Primer(template_sequence=self.template_sequence,
                                            reference_template_sequence=self.reference_template_sequence,
                                            sequence=self.primer3_result[f'PRIMER_LEFT_{primer3_rank}_SEQUENCE'],
                                            target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                                            strand='forward', primer_type='forward')    

                if self.primer3_result.get(f'PRIMER_RIGHT_{primer3_rank}') != None:
                    reverse_primer = Primer(template_sequence=self.template_sequence,
                                            reference_template_sequence=self.reference_template_sequence,
                                            sequence=self.primer3_result[f'PRIMER_RIGHT_{primer3_rank}_SEQUENCE'],
                                            target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                                            strand='reverse', primer_type='reverse')
                
                if self.primer3_result.get(f'PRIMER_INTERNAL_{primer3_rank}') != None:
                    probe = Primer(template_sequence=self.template_sequence,
                                reference_template_sequence=self.reference_template_sequence,
                                sequence=self.primer3_result[f'PRIMER_INTERNAL_{primer3_rank}_SEQUENCE'],
                                target_start_index=self.target_start_index, target_end_index=self.target_end_index,
                                strand='forward', primer_type='probe')
                
                amplicon = Amplicon(template_sequence=self.template_sequence,
                                    reference_template_sequence=self.reference_template_sequence,
                                    target_start_index=self.target_start_index, 
                                    target_end_index=self.target_end_index,
                                    forward_primer=forward_primer,
                                    reverse_primer=reverse_primer,
                                    probe=probe)
                print(amplicone)
                if probe != None:
                    probe_start = probe.template_sequence.find(probe.sequence)
                    probe_end = probe_start + len(probe.sequence)
                    if (probe_start<=self.target_start_index) and (probe_end>=(self.target_end_index)):
                        amplicon_list.append(amplicon)
                else:
                    amplicon_list.append(amplicon)
        
        self.amplicon_list = amplicon_list

    def design_primer(self):
        self.run_primer3()