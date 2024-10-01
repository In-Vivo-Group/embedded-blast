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

if __name__ == "__main__":
    query = "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG"
    database_file = "path_to_database.fasta"
    
    print("Running NCBI BLAST...")
    ncbi_results, ncbi_time = run_ncbi_blast(query, database_file)
    print(f"NCBI BLAST completed in {ncbi_time:.2f} seconds")
    
    print("\nRunning Embedding-based BLAST...")
    embedding_results, embedding_time = run_embedding_blast(query, database_file)
    print(f"Embedding-based BLAST completed in {embedding_time:.2f} seconds")
    
    print("\nComparing results...")
    sensitivity = compare_results(ncbi_results, embedding_results)
    print(f"Sensitivity (proportion of NCBI BLAST hits found): {sensitivity:.2f}")
    
    print(f"\nSpeed comparison: Embedding-based BLAST was {ncbi_time / embedding_time:.2f}x faster")