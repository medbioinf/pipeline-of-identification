docker-imgs:
	docker build -t medbioinf/ident-comparison-python:latest -f docker/python/Dockerfile .
	docker build -t medbioinf/ident-comparison-comet:latest -f docker/comet/Dockerfile .
	docker build -t medbioinf/ident-comparison-percolator:latest -f docker/percolator/Dockerfile .
	docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses