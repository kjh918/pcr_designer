import primer3
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction

def get_start_end_index(template_sequence, sequence):
    try:
        start_index = template_sequence.index(sequence)
    except:
        try:
            start_index = template_sequence.index(reverse_complement(sequence))
        except ValueError:
            print(f'{sequence} not in {template_sequence}')
    end_index = start_index + len(sequence) - 1
    return (start_index, end_index)


class Primer():
    
    template_sequence: str
    sequence: str
    strand: str
    tm: float
    gc_percent: float

    DEFAULT_SALT_MONOVALENT = 50
    DEFAULT_SALT_DIVALENT = 1.5
    DEFAULT_DNTP_CONC = 0.6
    DEFAULT_DNA_CONC = 50

    def __init__(
        self,
        template_sequence,
        sequence,
        strand,
        primer_type,
        target_start_index,
        target_end_index,
        reference_template_sequence=None,
        chrom=None,
        start=None,
        end=None,
        salt_monovalent_conc=DEFAULT_SALT_MONOVALENT,
        salt_divalent_conc=DEFAULT_SALT_DIVALENT,
        dntp_conc=DEFAULT_DNTP_CONC,
        dna_conc=DEFAULT_DNA_CONC,
    ):

        if reference_template_sequence is not None:
            self.reference_template_sequence = reference_template_sequence
        else:
            self.reference_template_sequence = template_sequence

        self.template_sequence = template_sequence
        self.sequence = sequence
        self.strand = strand
        self.primer_type = primer_type
        self.target_start_index = target_start_index
        self.target_end_index = target_end_index
        self.length = len(sequence)
        self.start_index, self.end_index = get_start_end_index(
            self.template_sequence, self.sequence
        )
        self.chrom = chrom
        self.start = start
        self.end = end

        # PCR reaction condition
        self.salt_monovalent_conc = salt_monovalent_conc
        self.salt_divalent_conc = salt_divalent_conc
        self.dntp_conc = dntp_conc
        self.dna_conc = dna_conc

        # 🔁 primer3 v2: calcTm -> calc_tm
        self.tm = primer3.calc_tm(
            self.sequence,
            mv_conc=salt_monovalent_conc,
            dv_conc=salt_divalent_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )

        self.gc_percent = gc_fraction(self.sequence, ambiguous='ignore') * 100

        # 🔁 primer3 v2: calcHairpin -> calc_hairpin
        primer3_hairpin_result = primer3.calc_hairpin(
            self.sequence,
            mv_conc=salt_monovalent_conc,
            dv_conc=salt_divalent_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
        self.hairpin = primer3_hairpin_result.structure_found
        self.hairpin_tm = primer3_hairpin_result.tm
        # dg, dh, ds 단위는 기존과 동일하게 ThermoResult에서 제공됩니다. 필요에 따라 1000으로 나누어 사용.
        self.hairpin_dg = primer3_hairpin_result.dg / 1000
        self.hairpin_dh = primer3_hairpin_result.dh / 1000
        self.hairpin_ds = primer3_hairpin_result.ds / 1000

        # 🔁 primer3 v2: calcHomodimer -> calc_homodimer
        primer3_homodimer_result = primer3.calc_homodimer(
            self.sequence,
            mv_conc=salt_monovalent_conc,
            dv_conc=salt_divalent_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
        self.homodimer = primer3_homodimer_result.structure_found
        self.homodimer_tm = primer3_homodimer_result.tm
        self.homodimer_dg = primer3_homodimer_result.dg / 1000
        self.homodimer_dh = primer3_homodimer_result.dh / 1000
        self.homodimer_ds = primer3_homodimer_result.ds / 1000

    def check_three_prime_is(self, sequence):
        if self.strand == 'forward':
            three_primer_sequence = self.reference_template_sequence[
                self.end_index + 1 - len(sequence) : self.end_index + 1
            ]
        elif self.strand == 'reverse':
            three_primer_sequence = reverse_complement(
                self.reference_template_sequence[
                    self.start_index : self.start_index + len(sequence)
                ]
            )
        if three_primer_sequence == sequence:
            return True
        else:
            return False

    def check_cpg_count(self, min_cpg_count, max_cpg_count):
        cpg_count = self.count_cpg()
        if min_cpg_count <= cpg_count <= max_cpg_count:
            return True
        else:
            return False

    def check_non_cpg_cytosine_count(
        self, min_non_cpg_cytosine_count, max_non_cpg_cytosine_count
    ):
        non_cpg_cytosine_count = self.count_non_cpg_cytosine()
        if min_non_cpg_cytosine_count <= non_cpg_cytosine_count <= max_non_cpg_cytosine_count:
            return True
        else:
            return False

    def count_cpg(self):
        try:
            if self.strand == 'forward':
                return self.reference_template_sequence[
                    self.start_index : self.end_index + 2
                ].count('CG')
            elif self.strand == 'reverse':
                return reverse_complement(
                    self.reference_template_sequence[
                        self.start_index - 1 : self.end_index + 1
                    ]
                ).count('CG')
        except:
            if self.strand == 'forward':
                return self.reference_template_sequence[
                    self.start_index : self.end_index + 1
                ].count('CG')
            elif self.strand == 'reverse':
                return reverse_complement(
                    self.reference_template_sequence[
                        self.start_index : self.end_index + 1
                    ]
                ).count('CG')

    def count_non_cpg_cytosine(self):
        if self.strand == 'forward':
            return (
                self.reference_template_sequence[
                    self.start_index : self.end_index + 1
                ].count('C')
                - self.count_cpg()
            )
        elif self.strand == 'reverse':
            return (
                reverse_complement(
                    self.reference_template_sequence[
                        self.start_index : self.end_index + 1
                    ]
                ).count('C')
                - self.count_cpg()
            )

    def to_dict(
        self,
        ignore_attributes=[
            'template_sequence',
            'reference_template_sequence',
            'primer_type',
            'chrom',
            'start',
            'end',
            'target_start_index',
            'target_end_index',
        ],
    ):
        primer_dict = {}
        for key, value in self.__dict__.items():
            if key in ignore_attributes:
                continue
            else:
                primer_dict[f'{self.primer_type}_{key}'] = value

        return primer_dict


class Amplicon():

    reference_template_sequence: str
    template_sequence: str
    target_start_index: int
    target_end_index: int
    chrom: str
    start: int
    end: int
    forward_primer: Primer
    reverse_primer: Primer
    probe: Primer
    amplicon_sequence: str

    def __init__(
        self,
        template_sequence,
        target_start_index,
        target_end_index,
        reference_template_sequence=None,
        chrom=None,
        start=None,
        end=None,
        forward_primer=None,
        reverse_primer=None,
        probe=None,
    ):

        if reference_template_sequence is not None:
            self.reference_template_sequence = reference_template_sequence
        else:
            self.reference_template_sequence = template_sequence

        self.template_sequence = template_sequence
        self.target_start_index = target_start_index
        self.target_end_index = target_end_index
        self.chrom = chrom
        self.start = start
        self.end = end

        self.forward_primer = forward_primer
        if forward_primer is not None:
            (
                self.forward_start_index,
                self.forward_end_index,
            ) = get_start_end_index(self.template_sequence, self.forward_primer.sequence)

        self.reverse_primer = reverse_primer
        if reverse_primer is not None:
            (
                self.reverse_start_index,
                self.reverse_end_index,
            ) = get_start_end_index(self.template_sequence, self.reverse_primer.sequence)

        self.probe = probe
        if probe is not None:
            (
                self.probe_start_index,
                self.probe_end_index,
            ) = get_start_end_index(self.template_sequence, self.probe.sequence)

        self.amplicon_sequence = self.cal_amplicon_sequence()

    def cal_amplicon_sequence(self):
        if self.forward_primer is not None and self.reverse_primer is not None:
            amplicon_sequence = self.template_sequence[
                self.forward_start_index : self.reverse_end_index + 1
            ]
            return amplicon_sequence
        else:
            return None

    def to_dict(self):
        amplicon_dict = {}
        amplicon_dict['reference_template_sequence'] = self.reference_template_sequence
        amplicon_dict['template_sequence'] = self.template_sequence
        amplicon_dict['target_start_index'] = self.target_start_index
        amplicon_dict['target_end_index'] = self.target_end_index
        try:
            amplicon_dict['amplicon_sequence'] = self.amplicon_sequence
            amplicon_dict['amplicon_length'] = len(self.amplicon_sequence)
        except:
            amplicon_dict['amplicon_sequence'] = None
            amplicon_dict['amplicon_length'] = None

        if self.forward_primer is not None:
            for key, value in self.forward_primer.to_dict().items():
                amplicon_dict[key] = value
        if self.reverse_primer is not None:
            for key, value in self.reverse_primer.to_dict().items():
                amplicon_dict[key] = value
        if self.probe is not None:
            for key, value in self.probe.to_dict().items():
                amplicon_dict[key] = value

        return amplicon_dict
