import numpy as np
from sklearn.neighbors import KDTree
from src.embedding.sequence_embeddings import generate_protein_embeddings, generate_nucleotide_embeddings

def generate_query_subseq_embeddings(query, window_size, step_size, embed_func):
    subseqs = [query[i:i+window_size] for i in range(0, len(query)-window_size+1, step_size)]
    return [embed_func(subseq) for subseq in subseqs]

def find_seeds(query_embeddings, database_embeddings, threshold):
    db_ids, db_embeds = zip(*database_embeddings.items())
    db_embeds = np.array(db_embeds)
    if db_embeds.ndim == 1:
        db_embeds = db_embeds.reshape(-1, 1)
    tree = KDTree(db_embeds)
    seeds = []
    for i, q_embed in enumerate(query_embeddings):
        distances, indices = tree.query([q_embed], k=5)  # Get top 5 closest matches
        for dist, idx in zip(distances[0], indices[0]):
            if dist < threshold:
                seeds.append((i, db_ids[idx], dist))
    return seeds

def extend_seed(query, db_seq, seed_pos, gap_open=-10, gap_extend=-1):
    # Simple extension algorithm (can be improved)
    left_score, left_query, left_db = 0, "", ""
    right_score, right_query, right_db = 0, "", ""
    
    # Extend left
    i, j = seed_pos, 0
    while i > 0 and j > 0:
        if query[i-1] == db_seq[j-1]:
            left_score += 1
            left_query = query[i-1] + left_query
            left_db = db_seq[j-1] + left_db
            i -= 1
            j -= 1
        else:
            break

    # Extend right
    i, j = seed_pos, 0
    while i < len(query) and j < len(db_seq):
        if query[i] == db_seq[j]:
            right_score += 1
            right_query += query[i]
            right_db += db_seq[j]
            i += 1
            j += 1
        else:
            break

    return left_query + right_query, left_db + right_db, left_score + right_score

def evaluate_alignment(query_seq, db_seq, gap_open=-10, gap_extend=-1):
    # Simple scoring function (can be improved)
    score = sum(1 if q == d else -1 for q, d in zip(query_seq, db_seq))
    gaps = query_seq.count('-') + db_seq.count('-')
    score += gap_open * gaps + gap_extend * (len(query_seq) - len(db_seq.replace('-', '')))
    return score

def embedding_based_blast(query, database, database_embeddings, embed_func, window_size=3, step_size=1, threshold=0.1):
    query_embeddings = generate_query_subseq_embeddings(query, window_size, step_size, embed_func)
    seeds = find_seeds(query_embeddings, database_embeddings, threshold)
    
    alignments = []
    for seed in seeds:
        query_pos, db_id, dist = seed
        db_seq = database[db_id]
        extended_query, extended_db, score = extend_seed(query, db_seq, query_pos)
        alignment_score = evaluate_alignment(extended_query, extended_db)
        alignments.append((db_id, extended_query, extended_db, alignment_score))
    
    return sorted(alignments, key=lambda x: x[3], reverse=True)

# Example usage
if __name__ == "__main__":
    from src.embedding.sequence_embeddings import generate_database_embeddings
    from Bio import SeqIO

    query = "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG"
    database = {record.id: str(record.seq) for record in SeqIO.parse("path_to_database.fasta", "fasta")}
    database_embeddings = generate_database_embeddings("path_to_database.fasta", "protein")
    results = embedding_based_blast(query, database, database_embeddings, generate_protein_embeddings)

    print(f"Number of alignments found: {len(results)}")
    for i, (db_id, q_seq, db_seq, score) in enumerate(results[:5]):  # Print top 5 alignments
        print(f"Alignment {i+1}:")
        print(f"Database sequence: {db_id}")
        print(f"Score: {score}")
        print(f"Query:  {q_seq}")
        print(f"Match:  {''.join('|' if q == d else ' ' for q, d in zip(q_seq, db_seq))}")
        print(f"Target: {db_seq}")
        print()