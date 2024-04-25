# Identification performance comparison

## Requirements
* Java (version depends on Nextflow)
* Nextflow
* Docker
* make

## Usage
1. Build and download the Docker images: `make docker-imgs`
2. Run the workflows
    ```shell
    nextflow run -profile docker src/main.nf --raws <FOLDER_WITH_RAW_FILES> ... 
    ```

## Development / Contribution
* Nextflow is going into `src/*`
* Each process should use a container
* Make sure to put the each container download or build into the Makefile
* Python scripts are going into `bin/` with the shebang `#!/usr/bin/env python`. Make it executable.
* Python dependencies are going into `requirements.txt`
* Any binary dependency for the Python Docker container are is going into `environment.yml`