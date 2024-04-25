docker-imgs:
	docker build -t medbioinforub/ident-comparison-python:latest -f docker/python/Dockerfile .
	docker build -t medbioinforub/ident-comparison-comet:latest -f docker/comet/Dockerfile .
	docker pull chambm/pwiz-skyline-i-agree-to-the-vendor-licenses