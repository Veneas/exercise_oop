# Class 01

# =-=-=-=-=-=-=-= Task 1 =-=-=-=-=-=-=-=

# Constructor + safety checks and warnings
class GenomicFeature:
    def __init__(self, chromosome, start, end, strand):
        if not isinstance(chromosome, str):
            raise ValueError("Chromosome is not a string")
        if not isinstance(start, int) or not isinstance (end, int):
            raise ValueError("'Start' and 'end' must be an integers")
        if start < 1 or end < 1:
            raise ValueError("'Start' and 'end' must be higher than 1")
        if start > end:
            raise ValueError("'Start' cannot be a higher number than an 'end'")
        if strand not in ["+", "-"]:
            raise ValueError("Strand is not defined '+' nor '-'")

        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand

# Methods
    def length(self):
        return self.end - self.start + 1

    def overlaps(self, other):
        return(
            self.chromosome == other.chromosome
            and self.start <= other.end
            and other.start <= self.end
        )

    def describe(self):
        return(
            f"{type(self).__name__} " 
            f"{self.chromosome}:{self.start}-{self.end}"
            f"({self.strand})"
        )



# =-=-=-=-=-=-=-= Task 2 =-=-=-=-=-=-=-=

class Exon(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, exon_number):
        super().__init__(chromosome, start, end, strand)
        self.exon_number = exon_number

    # Override
    def describe(self):
        return(
            f"{super().describe()}"
            f" exon #{self.exon_number}"
        )



# =-=-=-=-=-=-=-= Task 3 =-=-=-=-=-=-=-=
class Gene(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, name, exons):
        super().__init__(chromosome, start, end, strand)
        self.name = name
        self.exons = []
    def add_exon(self, exon):
        self.exons.append(exon)

    def total_exon_length(self):
        return sum(exon.length() for exon in self.exons)

    def describe(self):
        return(
            f"Gene {self.name}",
            f"{self.chromosome}:{self.start}-{self.end}",
            f"({self.strand}),"
            f"{len(self.exons)} exon(s)"
        )


class Variant(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, ref_allele, alt_allele):
        super().__init__(chromosome, start, end, strand)
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele


if __name__ == "__main__":
    a = GenomicFeature("chr1", 1000, 5000, "+")
    b = GenomicFeature("chr1", 4800, 6000, "+")
    c = GenomicFeature("chr2", 1000, 5000, "+")

    print(a.describe())     # GenomicFeature chr1:1000-5000(+)
    print(a.length())       # 4001
    print(a.overlaps(b))    # True (4800 - 5000 shared)
    print(a.overlaps(c))    # False (different chromosome) -> chr1 vs chr2

    # Raising Value error as the start > end
    #GenomicFeature("chr1", 5000, 1000, "+")

    features = [
        GenomicFeature("chr1", 1000, 5000, "+"),
        Exon("chr1", 1000, 1200, "+", 1),
        Exon("chr1", 3000, 3300, "+", 2)
    ]

    for feature in features:
        print(feature.describe())



