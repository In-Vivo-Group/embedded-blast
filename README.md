# BLAST in Embedding Space

## Project Overview

This project aims to develop a novel method for conducting the BLAST (Basic Local Alignment Search Tool) algorithm over embedding space rather than raw sequence data. By leveraging modern machine learning techniques and embedding methods, we seek to improve the speed and sensitivity of sequence similarity searches in bioinformatics.

## Table of Contents

1. [Project Goals](#project-goals)
2. [Project Structure](#project-structure)
3. [Prerequisites](#prerequisites)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Building and Running with Docker](#building-and-running-with-docker)
7. [Contributing](#contributing)
8. [Roadmap](#roadmap)
9. [License](#license)
10. [Contact](#contact)

## Project Goals

1. Develop a new BLAST-like algorithm that operates on sequence embeddings
2. Improve search speed and sensitivity compared to traditional BLAST
3. Create a user-friendly interface for the new algorithm
4. Publish findings in a peer-reviewed journal

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
└── README.md
```

- `src/`: Contains the source code for the project
  - `embedding/`: Code for sequence embedding methods
  - `algorithm/`: Implementation of the new BLAST algorithm
  - `benchmarking/`: Scripts for performance evaluation
- `data/`: Sample datasets and benchmark data
- `tests/`: Unit tests and integration tests
- `docs/`: Project documentation
- `notebooks/`: Jupyter notebooks for exploration and analysis
- `results/`: Benchmark results and performance analyses
- `Dockerfile`: Instructions for building the Docker image
- `requirements.txt`: List of Python dependencies
- `run.py`: Main script to run the benchmarking

## Prerequisites

- Python 3.9+
- Docker (for containerized usage)

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

1. Build the Docker image:
   ```
   docker build -t blast-embedding .
   ```

2. Run the Docker container:
   ```
   docker run -it --rm blast-embedding
   ```

   This will run the benchmarking script with a default query sequence.

3. To use a custom query sequence, you can pass it as an environment variable:
   ```
   docker run -it --rm -e QUERY_SEQUENCE="YOURSEQUENCEHERE" blast-embedding
   ```

### Customization

To use your own database, replace the `data/sample_database.fasta` file with your FASTA format database before building the Docker image.

## Contributing

We welcome contributions to the BLAST in Embedding Space project! Please follow these steps to contribute:

1. Fork the repository
2. Create a new branch (`git checkout -b feature/AmazingFeature`)
3. Make your changes
4. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
5. Push to the branch (`git push origin feature/AmazingFeature`)
6. Open a Pull Request

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## Roadmap

1. **Understanding Current BLAST** (Weeks 1-2)
   - Review BLAST literature
   - Analyze BLAST stages

2. **Choose Embedding Method** (Weeks 3-5)
   - Research embedding techniques
   - Implement and test methods

3. **Define Embedding Space** (Weeks 6-7)
   - Determine dimensionality
   - Select distance metrics

4. **Adapt BLAST Stages** (Weeks 8-13)
   - Develop seeding, extension, and scoring for embeddings

5. **Implementation and Optimization** (Weeks 14-17)
   - Develop full prototype
   - Optimize for efficiency

6. **Validation and Benchmarking** (Weeks 18-20)
   - Prepare datasets
   - Run comprehensive benchmarks

7. **Iteration and Refinement** (Weeks 21-24)
   - Analyze results
   - Refine algorithm

8. **Practical Considerations** (Weeks 25-27)
   - Design embedding database
   - Develop user interface

9. **Documentation and Publication** (Weeks 28-30)
   - Write documentation
   - Prepare scientific paper

## License

TBD

## Notes

This is a research project and the embedding-based BLAST is a proof-of-concept. It may not be as comprehensive or accurate as established BLAST implementations. The project is designed to explore new approaches to sequence similarity search and may evolve significantly as research progresses.