# Use an official Python runtime as the base image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container
COPY . /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir --only-binary=:all: \
        numpy==1.21.0 \
        scipy==1.7.0 \
        gensim==4.1.2 \
        scikit-learn==0.24.2 \
        transformers==4.10.0 \
        biopython==1.79 \
        torch==1.9.0

# Run the application
CMD ["python", "run.py"]