from dataclasses import dataclass

@dataclass(frozen=True)
class Variant:
    """
    Single variant model.
    index: 0-based coordinate on reference_template_sequence (ref allele).
    """
    index: int
    ref: str
    alt: str
