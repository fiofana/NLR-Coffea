# Installation
!pip install biopython

from Bio import Entrez
import xml.etree.ElementTree as ET

# Making a dictionary based on NLR-Annotator predictions
predicted_genes = {}

# Change directory name according to the chromosome
with open('posicoes_1c.txt', mode='r', encoding='utf-8') as prediction: 
    for line in prediction:
        line = line.strip()
        line = line.split('-')
        key = line[0]
        value = line[1:]
        predicted_genes[key] = value

# Verification: checking if list is correct
print(predicted_genes)

# Verification: the last ID (nlr_X) has to have X equal to seq. length
print('The dictionary has', len(predicted_genes.keys()), 'lists of predicted genes')

# Using Esearch to look for gene IDs in the predicted positions
def get_gene_ids(chromosome, gene_start, gene_end):
    gene_ids = []
    query = f"{chromosome}[Chromosome] AND ({gene_start}:{gene_end}[Base Position]) AND Coffea arabica[Organism]"
    # For troubleshooting
    # print(query)

    try:
        handle = Entrez.esearch(db="gene", term=query, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            gene_ids.extend(record["IdList"])

    except Exception as e:
        print(f"Error retrieving gene IDs for {chromosome}: {e}")

    return gene_ids

# Input
verified_ids = []
verified_genes = {}
only_predicted = {}

for key, value in predicted_genes.items():
    # Attributing variables for each element in the gene list for the current key
    chromosome, predicted_gene_start, predicted_gene_end, predicted_direction = value
    id_result = get_gene_ids(chromosome,  predicted_gene_start,  predicted_gene_end)

    if id_result:
        # Extending the value list with gene IDs  
        extended_value = value + [id_result]
        verified_genes[key] = extended_value

    else:
        only_predicted[key] = value

# Using Efetch to gather annotations using gene ID anotações usando o gene ID 
def get_annotations(gene_ids):
    annotations = []
    for gene_id in gene_ids:
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        annotation = handle.read()
        handle.close()
        annotations.append(annotation)
    return annotations

# Input
full_annotations = {}

for key, value in verified_genes.items():
    # Attributing variables for each element in the gene list for the current key
    chromosome, predicted_gene_start, predicted_gene_end, predicted_direction, real_ids = value
    annot_result = get_annotations(real_ids)
    annot_extended = value + [annot_result]
    full_annotations[key] = annot_extended

# Verification - identified genes number is maintained
print('There were retrieved', len(full_annotations.keys()), 'annotations')

# In case the verification goes wrong:

# Checking if the keys are the same
if set(full_annotations.keys()) == set(verified_genes.keys()):
    print("Both dictionaries have the same keys")
else:
    # Difference between the sets (keys exclusive to one)
    keys_difference = set(full_annotations.keys()) ^ set(verified_genes.keys())
    print(f"Different keys: {keys_difference}")

# Extracting relevant info from annotation 
def extract_gene_info(annotations):
    gene_info = []

    # Iterate over each annotation in case it's a list
    for annotation in annotations:
        # Parsing XML string in an ElementTree object
        annotation_tree = ET.fromstring(annotation)

        for entrezgene in annotation_tree.findall(".//Entrezgene"):
            gene_data = {}

            gene_id_element = entrezgene.find(".//Gene-track_geneid")
            if gene_id_element is not None:
                gene_data["gene_id"] = gene_id_element.text

            gene_start_element = entrezgene.find(".//Seq-interval_from")
            if gene_start_element is not None:
                gene_data["gene_start"] = gene_start_element.text

            gene_end_element = entrezgene.find(".//Seq-interval_to")
            if gene_end_element is not None:
                gene_data["gene_end"] = gene_end_element.text

            direction_element = entrezgene.find(".//Seq-interval_strand/Na-strand")
            if direction_element is not None:
                gene_data["direction"] = direction_element.attrib["value"]

            # Corresponding Protein
            protein_id_element = entrezgene.find(".//Gene-commentary/Gene-commentary_type[@value='peptide']/../Gene-commentary_accession")
            if protein_id_element is not None:
                gene_data["protein_id"] = protein_id_element.text

            gene_info.append(gene_data)

    return gene_info

# Input
verified_annotations = {}

for key, value in full_annotations.items():
    # Attributing variables for each element in the gene list for the current key
    chromosome, predicted_gene_start, predicted_gene_end, predicted_direction, real_ids, annotations = value
    
    # Extracting gene information
    gene_info = extract_gene_info(annotations)
    
    # Creating a complete dictionary to store the combined information
    combined_info = {
        "chromosome": chromosome,
        "predicted_gene_start": predicted_gene_start,
        "predicted_gene_end": predicted_gene_end,
        "predicted_direction": predicted_direction,
        "annotation": gene_info
    }
    verified_annotations[key] = combined_info

# Verification - identified genes number is maintained
print('There were fully validated', len(verified_annotations.keys()), 'genes')


