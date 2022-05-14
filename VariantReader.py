"""
Write an specialization of the VariantReader abstract class that is able to read genomic variants from a VCF file:
  * The name of the new class must be VcfVariantReader
  * The API must implement the spec defined by the VariantReader class
  * Returned variants must contain, parsed from the VCF, the following:
    - Chromosome
    - Position
    - Reference allele
    - Alternate allele
    - FILTER string
    - INFO string
    - For each of the samples in the VCF, the genotype (and only the genotype, no further sample data is required)
  * You can create additional private methods for either the VariantReader or the VcfVariantReader class, if necessary
  * You can safely assume the VCF will *only* contain *one variant per line*, i.e. no parsing of multiallelic positions
  is required
  * The Variant data model is provided in a separate file "variantschema.proto" and defined using Google protocol
  buffers https://developers.google.com/protocol-buffers
  * No comprehensive validation of the VCF file format is expected
  * Bear in mind that efficiency is critical, as a real multi-sample VCF might contain tens of millions of variants.
  * The tests written below must complete successfully at the end of the exercise
  * Please, write the code before the "TESTS" section
"""
from abc import abstractmethod


class VariantReader:
    """
    A VariantReader will enable to read Variant objects from specific file formats
    """

    def __init__(self, filename: str) -> None:
        self._filename = filename

    @abstractmethod
    def pre(self) -> None:
        """
        This method can be implemented to parse and process metadata from file headers. Might not be strictly necessary
        to implement if no metadata is needed to be used from the header or no header is present at all
        """
        pass

    @abstractmethod
    def read(self) -> Variant:
        """
        This method reads one single Variant object each time is called.
        :return: a Variant object each time the method is called while there remain variants to read in the file
        :rtype Variant
        """
        pass

    @abstractmethod
    def post(self) -> None:
        """
        This method can be implemented to perform some post-processing tasks once all variants in the file have been
        parsed. Might not be strictly necessary to implement if no post-processing actions are required.
        """
        pass


""" *******************************************************************************************************************
****************************************  TESTS  **********************************************************************
******************************************************************************************************************* """

with VcfVariantReader("test.vcf.gz") as vcf_variant_reader:
    vcf_variant_reader.pre()
    variant_list = [variant for variant in vcf_variant_reader.read()]
    vcf_variant_reader.post()

# Three variants properly read from the VCF
assert len(variant_list) == 3
# All three are Variant objects
assert all([isinstance(variant, Variant) for variant in variant_list])


def assert_contains(
    variant_list,
    chromosome,
    position,
    reference,
    alternate,
    filter,
    info,
    genotype_sample1,
    genotype_sample2,
    genotype_sample3,
):
    for variant in variant_list:
        if (
            variant.chromosome == chromosome
            and variant.position == position
            and variant.reference == reference
            and variant.alternate == alternate
            and variant.filter == filter
            and variant.info == info
            and variant.genotypes
            and len(variant.genotypes) == 3
            and variant.genotypes["SAMPLE1"]
            and variant.genotypes["SAMPLE1"] == genotype_sample1
            and variant.genotypes["SAMPLE2"]
            and variant.genotypes["SAMPLE2"] == genotype_sample2
            and variant.genotypes["SAMPLE3"]
            and variant.genotypes["SAMPLE3"] == genotype_sample3
        ):
            assert True
            return

    assert False


# Check all relevant info for the first variant was properly parsed and loaded
assert_contains(
    variant_list,
    chromosome="chr1",
    position=10146,
    reference="AC",
    alternate="A",
    filter="badReads;MQ",
    info="FR=0.7127;MMLQ=14.0;TCR=1;HP=4;WE=10154;Source=Platypus;WS=10136;PP=188.0;TR=5;NF=4;TCF=5;NR=1"
    ";TC=6;MGOF=666;SbPval=0.59;MQ=8.81;QD=39.6;SC=CCTAACCCTAACCCCTAACCC;BRF=0.99;HapScore=2",
    genotype_sample1="1/1",
    genotype_sample2="0/1",
    genotype_sample3="0/1",
)

# Check all relevant info for the second variant was properly parsed and loaded
assert_contains(
    variant_list,
    chromosome="chr1",
    position=10430,
    reference="C",
    alternate="CCT",
    filter="PASS",
    info="FR=0.6232;MMLQ=7.0;TCR=7;HP=2;WE=10441;Source=Platypus;WS=10420;PP=230.0;TR=6;NF=1;TCF=1;NR=5;TC"
    "=8;MGOF=206;SbPval=0.58;MQ=19.98;QD=40.0;SC=AACCCTAACCCTAACCCTAAC;BRF=0.96;HapScore=8",
    genotype_sample1="1/1",
    genotype_sample2="1/0",
    genotype_sample3="./.",
)

# Check all relevant info for the third variant was properly parsed and loaded
assert_contains(
    variant_list,
    chromosome="chr1",
    position=10431,
    reference="T",
    alternate="A",
    filter="PASS",
    info="FR=0.6231;MMLQ=7.0;TCR=7;HP=3;WE=10441;Source=Platypus;WS=10420;PP=186.0;TR=3;NF=0;TCF=1;NR=3;TC"
    "=8;MGOF=206;SbPval=1.0;MQ=19.98;QD=79.8226;SC=ACCCTAACCCTAACCCTAACC;BRF=0.96;HapScore=8"
    ";OLD_CLUMPED=chr1:10431:TAA/AAC",
    genotype_sample1="1|1",
    genotype_sample2="1|0",
    genotype_sample3=".|.",
)

print("all tests pass")
