import time
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from src.algorithm.refined_embedding_blast import embedding_based_blast
from src.embedding.sequence_embeddings import generate_protein_embeddings, generate_database_embeddings

def run_ncbi_blast(query, database_file):
    start_time = time.time()
    result_handle = NCBIWWW.qblast("blastp", "nr", query, entrez_query=f'"{database_file}"[PACC]')
    blast_record = NCBIXML.read(result_handle)
    end_time = time.time()
    
    alignments = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            alignments.append((alignment.title, hsp.query, hsp.sbjct, hsp.score))
    
    return alignments, end_time - start_time

def run_embedding_blast(query, database_file):
    start_time = time.time()
    database = {record.id: str(record.seq) for record in SeqIO.parse(database_file, "fasta")}
    database_embeddings = generate_database_embeddings(database_file, "protein")
    results = embedding_based_blast(query, database, database_embeddings, generate_protein_embeddings)
    end_time = time.time()
    
    return results, end_time - start_time

def compare_results(ncbi_results, embedding_results):
    ncbi_hits = set(result[0] for result in ncbi_results)
    embedding_hits = set(result[0] for result in embedding_results)
    
    common_hits = ncbi_hits.intersection(embedding_hits)
    ncbi_only = ncbi_hits - embedding_hits
    embedding_only = embedding_hits - ncbi_hits
    
    print(f"Common hits: {len(common_hits)}")
    print(f"NCBI BLAST only: {len(ncbi_only)}")
    print(f"Embedding BLAST only: {len(embedding_only)}")
    
    return len(common_hits) / len(ncbi_hits) if ncbi_hits else 0