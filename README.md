# BLAST in Embedding Space

## Project Overview

This project develops a novel approach to sequence similarity searching in bioinformatics by implementing a BLAST-like algorithm that operates on sequence embeddings rather than raw sequence data. By leveraging modern machine learning techniques and embedding methods, we aim to improve both the speed and sensitivity of sequence similarity searches compared to traditional BLAST (Basic Local Alignment Search Tool) algorithms.

## Table of Contents

1. [Background](#background)
2. [Methods](#methods)
3. [Project Goals](#project-goals)
4. [Project Structure](#project-structure)
5. [Prerequisites](#prerequisites)
6. [Installation](#installation)
7. [Usage](#usage)
8. [Building and Running with Docker](#building-and-running-with-docker)
9. [Contributing](#contributing)
10. [License](#license)
11. [Contact](#contact)

## Background

BLAST is a fundamental tool in bioinformatics for comparing biological sequence information, such as DNA sequences of genes or amino acid sequences of proteins. Traditional BLAST algorithms work directly on the sequence data, using heuristics to find regions of local similarity between sequences. While effective, these methods can be computationally intensive for large databases.

Recent advancements in machine learning, particularly in natural language processing, have shown that embedding techniques can capture complex relationships in sequential data. This project applies similar principles to biological sequences, hypothesizing that conducting similarity searches in embedding space could lead to faster and potentially more sensitive results.

## Methods

Our approach consists of several key components:

1. **Sequence Embedding**: We use advanced embedding techniques, primarily the ProtBERT model, to convert amino acid sequences into high-dimensional vector representations. These embeddings aim to capture the functional and structural properties of the proteins.

2. **Embedding-based Seeding**: Instead of using k-mers as in traditional BLAST, we identify seed regions by finding similar subsequences in the embedding space using efficient nearest neighbor search algorithms (KD-Tree).

3. **Alignment Extension**: We extend the seeds to form larger alignments, adapting traditional dynamic programming approaches to work with embedded representations.

4. **Scoring**: We develop a scoring system that combines similarity in embedding space with biologically relevant scoring matrices.

5. **Database Indexing**: We create an efficient index of the embedded database sequences to enable rapid searching.

## Project Structure

The repository is organized as follows:

```
blast-embedding/
├── src/
│   ├── embedding/
│   │   └── sequence_embeddings.py
│   ├── algorithm/
│   │   └── refined_embedding_blast.py
│   └── benchmarking/
│       └── blast_comparison.py
├── data/
│   └── sample_database.fasta
├── tests/
├── docs/
├── notebooks/
├── results/
├── Dockerfile
├── requirements.txt
├── run.py
├── embedded_blast.ipynb
└── README.md
```

## Prerequisites

- Python 3.9+
- Docker (for containerized usage)
- 8GB+ RAM recommended for running embedding models

## Installation

To set up the development environment:

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/blast-embedding.git
   cd blast-embedding
   ```

2. Create a virtual environment:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   ```

3. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

To run the benchmarking script locally:

```
python run.py
```

This will run the embedding-based BLAST and compare it with NCBI BLAST using a default query sequence.

## Building and Running with Docker

Docker provides an isolated environment to run the project, ensuring consistency across different systems. Follow these steps to build and run the project using Docker:

1. Build the Docker image:
   ```
   docker build -t blast-embedding .
   ```
   This command builds a Docker image named 'blast-embedding' based on the instructions in the Dockerfile.

2. Run the Docker container:
   ```
   docker run -it --rm blast-embedding
   ```
   This command starts a container from the 'blast-embedding' image, runs the benchmarking script with a default query sequence, and removes the container after execution.

3. To use a custom query sequence:
   ```
   docker run -it --rm -e QUERY_SEQUENCE="YOURSEQUENCEHERE" blast-embedding
   ```
   Replace "YOURSEQUENCEHERE" with your actual protein sequence.

4. To run an interactive shell in the container:
   ```
   docker run -it --rm --entrypoint /bin/bash blast-embedding
   ```
   This allows you to explore the container's file system and run commands manually.

5. To mount a local directory and save results:
   ```
   docker run -it --rm -v /path/to/local/directory:/app/results blast-embedding
   ```
   Replace "/path/to/local/directory" with the actual path on your host machine.

### Customization

To use your own database:
1. Replace the `data/sample_database.fasta` file with your FASTA format database.
2. Rebuild the Docker image:
   ```
   docker build -t blast-embedding .
   ```

## Contributing

We welcome contributions to the BLAST in Embedding Space project! Please follow these steps to contribute:

1. Fork the repository
2. Create a new branch (`git checkout -b feature/AmazingFeature`)
3. Make your changes
4. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
5. Push to the branch (`git push origin feature/AmazingFeature`)
6. Open a Pull Request

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact

Alexander Titus - Send me a note on [LinkedIn](https://www.linkedin.com/in/alexandertitus/)

Project Link: [https://github.com/In-Vivo-Group/embedded-blast](https://github.com/In-Vivo-Group/embedded-blast)

## Acknowledgments

- This project builds upon the work of many researchers in the fields of bioinformatics and machine learning.
- We thank the developers of the ProtBERT model and other open-source tools used in this project.

## Notes

This is a research project and the embedding-based BLAST is a proof-of-concept. It may not be as comprehensive or accurate as established BLAST implementations. The project is designed to explore new approaches to sequence similarity search and may evolve significantly as research progresses.
