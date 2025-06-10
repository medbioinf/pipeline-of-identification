docker-imgs:
	docker build -t medbioinf/ident-comparison-python:latest -f docker/python/Dockerfile .

	docker pull proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.25073-842baef

	docker pull quay.io/medbioinf/openms:3.1.0

	docker pull quay.io/medbioinf/comet-ms:v2024.01.0
	docker pull quay.io/medbioinf/sage:v0.14.7
	docker pull quay.io/medbioinf/xtandem:2017.2.1.4
	docker pull quay.io/medbioinf/msamanda:3.0.22.071
	docker pull quay.io/medbioinf/msgfplus:v2024.03.26
	
	docker pull quay.io/medbioinf/mzid-merger:1.4.26

	docker pull quay.io/medbioinf/percolator:3.6.5
