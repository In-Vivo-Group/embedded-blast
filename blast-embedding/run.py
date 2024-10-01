import os
from src.benchmarking.blast_comparison import run_embedding_blast, run_ncbi_blast, compare_results

def main():
    # Load the query sequence from an environment variable or use a default
    query = os.environ.get('QUERY_SEQUENCE', "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG")
    
    # Use the sample database provided in the Docker image
    database_file = "data/sample_database.fasta"
    
    print("Running Embedding-based BLAST...")
    embedding_results, embedding_time = run_embedding_blast(query, database_file)
    print(f"Embedding-based BLAST completed in {embedding_time:.2f} seconds")
    
    print("\nRunning NCBI BLAST...")
    ncbi_results, ncbi_time = run_ncbi_blast(query, database_file)
    print(f"NCBI BLAST completed in {ncbi_time:.2f} seconds")
    
    print("\nComparing results...")
    sensitivity = compare_results(ncbi_results, embedding_results)
    print(f"Sensitivity (proportion of NCBI BLAST hits found): {sensitivity:.2f}")
    
    print(f"\nSpeed comparison: Embedding-based BLAST was {ncbi_time / embedding_time:.2f}x faster")

if __name__ == "__main__":
    main()