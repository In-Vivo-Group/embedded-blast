import numpy as np
from transformers import BertModel, BertTokenizer
from gensim.models import Word2Vec
from Bio import SeqIO

def generate_protein_embeddings(sequence, model_type='ProtBERT'):
    if model_type == 'ProtBERT':
        tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
        model = BertModel.from_pretrained("Rostlab/prot_bert")
        inputs = tokenizer(sequence, return_tensors='pt')
        outputs = model(**inputs)
        embeddings = outputs.last_hidden_state.mean(dim=1).squeeze().detach().numpy()
    elif model_type == 'ProtVec':
        # Assuming ProtVec model is pre-trained and loaded
        protvec_model = Word2Vec.load('path_to_protvec_model')
        words = [sequence[i:i+3] for i in range(len(sequence)-2)]
        embeddings = np.mean([protvec_model.wv[word] for word in words if word in protvec_model.wv], axis=0)
    return embeddings

def generate_nucleotide_embeddings(sequence, model_type='DNA2Vec'):
    if model_type == 'DNA2Vec':
        # Assuming DNA2Vec model is pre-trained and loaded
        dna2vec_model = Word2Vec.load('path_to_dna2vec_model')
        words = [sequence[i:i+11] for i in range(len(sequence)-10)]
        embeddings = np.mean([dna2vec_model.wv[word] for word in words if word in dna2vec_model.wv], axis=0)
    return embeddings

def generate_database_embeddings(fasta_file, sequence_type='protein'):
    embeddings = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if sequence_type == 'protein':
            embeddings[record.id] = generate_protein_embeddings(str(record.seq))
        else:
            embeddings[record.id] = generate_nucleotide_embeddings(str(record.seq))
    return embeddings