import numpy as np
from sklearn.neighbors import KDTree
from src.embedding.sequence_embeddings import generate_protein_embeddings, generate_nucleotide_embeddings

def generate_query_subseq_embeddings(query, window_size, step_size, embed_func):
    subseqs = [query[i:i+window_size] for i in range(0, len(query)-window_size+1, step_size)]
    return [embed_func(subseq) for subseq in subseqs]

def find_seeds(query_embeddings, database_embeddings, threshold):
    db_ids, db_embeds = zip(*database_embeddings.items())
    tree = KDTree(db_embeds)
    seeds = []
    for i, q_embed in enumerate(query_embeddings):
        distances, indices = tree.query([q_embed], k=10)  # Get top 10 closest matches
        for dist, idx in zip(distances[0], indices[0]):
            if dist < threshold:
                seeds.append((i, db_ids[idx], dist))
    return seeds

def embedding_based_blast(query, database_embeddings, embed_func, window_size=3, step_size=1, threshold=0.1):
    query_embeddings = generate_query_subseq_embeddings(query, window_size, step_size, embed_func)
    seeds = find_seeds(query_embeddings, database_embeddings, threshold)
    
    # Next steps would be to extend these seeds and perform the rest of the BLAST algorithm
    # ...

    return seeds  # For now, just return the seeds

# Example usage
if __name__ == "__main__":
    from src.embedding.sequence_embeddings import generate_database_embeddings

    query = "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG"
    database_embeddings = generate_database_embeddings("path_to_database.fasta", "protein")
    results = embedding_based_blast(query, database_embeddings, generate_protein_embeddings)

    print(f"Number of seeds found: {len(results)}")
    for seed in results[:10]:  # Print first 10 seeds
        print(f"Query position: {seed[0]}, Database sequence: {seed[1]}, Distance: {seed[2]}")