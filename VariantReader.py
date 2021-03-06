#!/usr/bin/env python
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
from variantschema_pb2 import Variant
import gzip

GZIP_EXTENSIONS = ('.gz', '.gzip')

class VariantReader:
    """
    A VariantReader will enable to read Variant objects from specific file formats
    """

    def __init__(self, filename: str) -> None:
        self._filename    = filename
        self._filehandler = None

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
        File is already closed automatically and each Variant automatically removed through the iterator returned from the `read` method.
        """
        pass

    def __enter__(self):
        """
        Implemented context manager as described by the doc: https://docs.python.org/2/reference/datamodel.html#with-statement-context-managers
        Essential for the `with` keyword.
        """
        if list(filter(self._filename.endswith, GZIP_EXTENSIONS)):
            self._filehandler = gzip.open(self._filename, "rt", encoding='utf8')
        else:
            self._filehandler = open(self._filename, "r", encoding='utf8')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Same as __enter__ for context manager
        """
        self._filehandler.close()


class VcfVariantReader(VariantReader):
    # Fixed column list as specified by the VCF documentation.
    VCF_COLUMNS      = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    VCF_COL_IDX_MAP  = { colname: index for index, colname in enumerate(VCF_COLUMNS) }
    VCF_HEADER_START = '#'
    VCF_COLUMN_SEP   = '\t'
    VCF_GENOTYPE_SEP = ':'

    def __init__(self, filename: str) -> None:
        super().__init__(filename)
        self._headerlines   = []
        self._headerColumns = []
        self._sampleNames   = []

    def pre(self) -> None:
        """
        Reads the header and stores it for the header() method, although not necessary for this test.
        Advances the file handler to the first variant.
        Validates the header, althouth not necessary for this test.
        Processes the columns metadata.
        """
        self._readHeader()
        self._validateHeader()
        self._processHeaderColumns()


    def _readHeader(self) -> None:
        """
        Stores the header.
        Assumes the header rigorously precedes the variants.
        """
        while line := self._filehandler.readline():
            if line.startswith( self.VCF_HEADER_START ):
                self._headerlines.append(line)
                # If last header line
                if line.startswith( self.VCF_COLUMNS[0] ):
                    break

    def _validateHeader(self) -> None:
        """
        Assumes the header is well formed
        """
        pass

    def _processHeaderColumns(self) -> None:
        """
        Finds and stores the VCF base column names and sample names.
        """
        self._headerColumns = self._headerlines[-1] \
                                .lstrip( self.VCF_HEADER_START ).strip() \
                                .split ( self.VCF_COLUMN_SEP )
        self._sampleNames = self._headerColumns[ len(self.VCF_COLUMNS) : ]

    def header(self) -> str:
        """
        Returns the VCF header as a string
        """
        return ''.join(self._headerlines)

    def _readVariantLine(self) -> str:
        """
        Reads the next variant line in the vcf file.
        Assumes the _filehandler has already passed the header with `self.pre()`
        :return: variant line
        :rtype str
        """
        return self._filehandler.readline()

    def read(self) -> Variant:
        """
        Yields each single Variant object as an iterator.
        :return: a Variant object each time the method is called while there remain variants to read in the file
        :rtype Variant
        """
        while varline := self._readVariantLine():
            variantFields = varline.split( self.VCF_COLUMN_SEP )
            variant = Variant(chromosome= str(variantFields[ self.VCF_COL_IDX_MAP['#CHROM'] ] ),
                              position  = int(variantFields[ self.VCF_COL_IDX_MAP['POS']   ] ),
                              reference = str(variantFields[ self.VCF_COL_IDX_MAP['REF']   ] ),
                              alternate = str(variantFields[ self.VCF_COL_IDX_MAP['ALT']   ] ),
                              filter    = str(variantFields[ self.VCF_COL_IDX_MAP['FILTER']] ),
                              info      = str(variantFields[ self.VCF_COL_IDX_MAP['INFO']  ] ),
                              genotypes = self.getVariantGenotypesMap( *variantFields[ len(self.VCF_COLUMNS) : ] )
                      )
            yield(variant)

    def getVariantGenotypesMap(self, *genotypeStrings):
        """
        Maps the header's list of samples with the genotype from the provided list of genotype strings.
        """
        return { str(sampleName): str(genotypeString.split( self.VCF_GENOTYPE_SEP ).pop(0))
                 for sampleName, genotypeString in zip(self._sampleNames, genotypeStrings)
               }

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
