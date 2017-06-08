import vcf
from rdflib import Namespace, Graph, URIRef, BNode, Literal
from rdflib.namespace import DCTERMS, RDFS, RDF, DC
import urllib
import sys
from SPARQLWrapper import SPARQLWrapper, JSON
from hashlib import md5

chrom = dict()
chrom["X"] = "NC_000023.11"
chrom["Y"] = "NC_000024.10"
chrom["14"] = "NC_000014.9"
chrom["17"] = "NC_000017.11"
chrom["1"] = "NC_000001.11"
chrom["19"] = "NC_000019.10"
chrom["6"] = "NC_000006.12"
chrom["8"] = "NC_000008.11"
chrom["2"] = "NC_000002.12"
chrom["7"] = "NC_000007.14"
chrom["20"] = "NC_000020.11"
chrom["3"] = "NC_000003.12"
chrom["16"] = "NC_000016.10"
chrom["21"] = "NC_000021.9"
chrom["22"] = "NC_000022.11"
chrom["15"] = "NC_000015.10"
chrom["18"] = "NC_000018.10"
chrom["4"] = "NC_000004.12"
chrom["9"] = "NC_000009.12"
chrom["13"] = "NC_000013.11"
chrom["10"] = "NC_000010.11"
chrom["5"] = "NC_000005.10"
chrom["11"] = "NC_000011.10"
chrom["12"] = "NC_000012.12"
chrom["MT"] = "NC_012920.1"

vcfGraph = Graph()
vcfGraph.bind("dcterms", DCTERMS)
vcfGraph.bind("wikidata_prop", URIRef("http://www.wikidata.org/prop/direct/"))
wikidataprop = Namespace("http://www.wikidata.org/prop/direct/")

# vcf_reader = vcf.Reader(open('/Users/andra/Downloads/CGC_flagship.missense_variants_snpEff_snpSift_GoNLv5.vcf', 'r'))
vcf_reader = vcf.Reader(open('../../CGC_flagship.missense_variants_snpEff_snpSift_GoNLv5.header.vcf', 'r'))
cgc_uri = "http://purl.org/fair/cgc"
oncoxl_uri = "http://oncoxl.fair-dtls.surf-hosted.nl"
sample = vcf_reader.samples[0]
sample_md5 = md5(sample.encode()).hexdigest()

analysis_uri = URIRef(oncoxl_uri + "/rdf/analysis/" + sample_md5)
sample_uri = URIRef(cgc_uri + "/sample/" + sample_md5)

# Annotate analysis
vcfGraph.add((analysis_uri, RDF.type, URIRef("http://edamontology.org/operation_2478")))
vcfGraph.add((analysis_uri, URIRef("http://www.w3.org/ns/prov#used"), URIRef("http://edamontology.org/topic_3673")))
vcfGraph.add((analysis_uri, URIRef("http://www.w3.org/ns/prov#used"), URIRef("https://zenodo.org/record/495587#.WTkUtcmxWL8")))

# Get Ensembl gene ID URI
# sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
#
# ensemblURI = dict()
# ensembl_geneQuery = """SELECT ?item ?itemLabel ?ensemblGeneID
#            WHERE
#            {
#                ?item wdt:P594 ?ensemblGeneID ;
#                      wdt:P703 wd:Q15978631 .
#                SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
#            }
#            """
# print(ensembl_geneQuery)
# sparql.setQuery(ensembl_geneQuery)
# sparql.setReturnFormat(JSON)
# results = sparql.query().convert()
# for result in results["results"]["bindings"]:
#     print(result["item"]["value"], result["ensemblGeneID"]["value"])
#     ensemblURI[result["ensemblGeneID"]["value"]] = result["item"]["value"]

for record in vcf_reader:
    chrom_nr = chrom[record.CHROM]
    variant_hgvs = chrom_nr+":g."+str(record.POS)+str(record.REF)+">"+str(record.ALT[0])
    print("hgvs: " + variant_hgvs)
    variant_uri = URIRef(oncoxl_uri + "/rdf/variant/" + urllib.parse.quote_plus(variant_hgvs))
    vcfGraph.add((variant_uri, RDF.type, URIRef("http://purl.obolibrary.org/obo/SO_0001060")))
    vcfGraph.add((variant_uri, DCTERMS.identifier, Literal(variant_hgvs)))
    vcfGraph.add((variant_uri, URIRef("http://www.wikidata.org/prop/direct/P3331"), Literal(variant_hgvs)))
    vcfGraph.add((variant_uri, URIRef("http://www.wikidata.org/prop/direct/P2576"), Literal("hg19")))
    chromosomeIRI = URIRef(oncoxl_uri + "/rdf/chromosome/"+chrom_nr)
    vcfGraph.add((chromosomeIRI, RDF.type, URIRef("https://www.wikidata.org/wiki/Q37748")))
    vcfGraph.add((chromosomeIRI, DCTERMS.identifier, Literal(chrom_nr)))
    vcfGraph.add((variant_uri, DCTERMS.isPartOf, chromosomeIRI))

    # Genomic START
    vcfGraph.add((variant_uri, wikidataprop.P644, Literal(record.POS)))
    vcfGraph.add((variant_uri, wikidataprop.P645, Literal(record.POS)))
    # print(record.)
    vcfInfo = record.INFO['ANN'][0].split("|")

    gene_uri = URIRef("http://rdf.ebi.ac.uk/resource/ensembl/"+vcfInfo[4])
    # print(record.INFO['ANN'][0])

    vcfGraph.add((gene_uri, DCTERMS.identifier, Literal(vcfInfo[4])))
    vcfGraph.add((variant_uri, URIRef("http://www.wikidata.org/prop/direct/P3433"), gene_uri))

    transcript_uri = URIRef("http://rdf.ebi.ac.uk/resource/ensembl.transcript/"+vcfInfo[6])
    vcfGraph.add((transcript_uri, DCTERMS.identifier, Literal(vcfInfo[6])))
    vcfGraph.add((transcript_uri, URIRef("http://purl.obolibrary.org/obo/so#transcribed_from"), gene_uri))

    # Add samples
    measurement_uri = URIRef(oncoxl_uri + "/rdf/measurement/" + sample_md5 + "/" + urllib.parse.quote_plus(variant_hgvs))
    vcfGraph.add((measurement_uri, RDF.type, URIRef("http://www.ebi.ac.uk/efo/EFO_0001444")))
    vcfGraph.add((measurement_uri, URIRef("http://semanticscience.org/resource/SIO_000300"), Literal(record.samples[0]['GT'])))
    vcfGraph.add((measurement_uri, URIRef("http://semanticscience.org/resource/SIO_000628"), variant_uri))

    vcfGraph.add((sample_uri, URIRef("http://purl.obolibrary.org/obo/GENO_0000222"), measurement_uri))

    # Link measurement and analysis
    vcfGraph.add((analysis_uri, URIRef("http://semanticscience.org/resource/SIO_000628"), measurement_uri))
    vcfGraph.add((analysis_uri, URIRef("http://www.w3.org/ns/prov#used"), sample_uri))

vcfGraph.serialize(destination='UMC_gene.ttl', format='turtle')
