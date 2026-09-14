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
    def __init__(self, chromosome, start, end, strand, name):
        super().__init__(chromosome, start, end, strand)
        self.name = name
        self.exons = []
    def add_exon(self, exon):
        self.exons.append(exon)

    def total_exon_length(self):
        return sum(exon.length() for exon in self.exons)

    def describe(self):
        return(
            f"Gene {self.name} "
            f"{self.chromosome}:{self.start}-{self.end}"
            f"({self.strand}), "
            f"{len(self.exons)} exon(s)"
        )


class Variant(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, ref_allele, alt_allele):
        super().__init__(chromosome, start, end, strand)
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
    def variant_type(self):
        if len(self.ref_allele) == 1 and len(self.alt_allele) == 1:
            return "SNP"
        elif len(self.ref_allele) < len(self.alt_allele):
            return "insertion"
        elif len(self.ref_allele) > len(self.alt_allele):
            return "deletion"
        else:
            return "MNV"
    def describe(self):
        return (
            f"Variant "
            f"{self.chromosome}:{self.start}-{self.end}"
            f"({self.strand}) "
            f"{self.ref_allele} > {self.alt_allele} "
            f"({self.variant_type()})"
        )

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
    
    # Loading the .tsv file 
    genes = {}          # Dict for genes
    variants = []       # List for variants
    report_items = []   # Combined list of genes and variants for the report

    with open("oop_data.tsv", "r") as f:
        # Skipping the header line
        next(f)
        # Reading the rest line by line
        for line in f:
            # Splitting the row into columns 
            fields = line.strip().split("\t")

            record_type = fields[0]
            chromosome = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            strand = fields[4]
            field_a = fields[5]

            # Creating objects 
            if record_type == "gene":
                gene = Gene(
                    chromosome, start, end, strand, field_a
                )
                # Saving gene in dict 
                genes[field_a] = gene
                # Add to report 
                report_items.append(gene)
            # creating exon objects 
            elif record_type == "exon":
                exon_number = int(fields[6])
                parent_gene_name = fields[5]
                exon = Exon(
                    chromosome, start, end, strand, exon_number
                )
                # Adding exon to the parent gene
                genes[parent_gene_name].add_exon(exon)
            # creating variant objects
            elif record_type == "variant":
                variant = Variant(
                    chromosome, start, end, strand, field_a, fields[6]
                )
                # Saving variant in list 
                variants.append(variant)
                # Add to report 
                report_items.append(variant)
        # Report section 
        print("\n Report:")
        print("--------")
        # Polymorphism ->
        # Every object gets same method call, but each class uses its own describe()
        for item in report_items:
            print(item.describe())
            # Only genes have total_exon_length() method, so checking the type of the object
            if isinstance(item, Gene):
                print(f"Total exon length: {item.total_exon_length()}")
        
        # Variant location analysis 
        print("\n Variant location analysis:")
        print("--------------------------")
        for variant in variants:
            print(f"{variant.describe()}")
            genes_found = False 
            # Checking every gene 
            for gene in genes.values():
                # Is variant inside this gene?
                if variant.overlaps(gene):
                    genes_found = True
                    # Checking which exons of this gene overlap with variant
                    overlapping_exons = [exon for exon in gene.exons if variant.overlaps(exon)]

                    if overlapping_exons:
                        for exon in overlapping_exons:
                            print(f" -> located in gene {gene.name}, exon # {exon.exon_number}")
                    else:
                        print(f" -> located in gene {gene.name}, but NOT in any exon")
            if not genes_found:
                print(" -> intergenic")


