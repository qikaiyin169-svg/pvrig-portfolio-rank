.PHONY: setup fetch replay reproduce test export

PYTHON ?= python3

setup:
	$(PYTHON) -m pip install -e .

fetch:
	$(PYTHON) scripts/fetch_inputs.py

replay:
	$(PYTHON) scripts/reproduce.py --mode replay

reproduce:
	$(PYTHON) scripts/reproduce.py --mode full

test:
	$(PYTHON) -m unittest discover -s tests -v

export:
	$(PYTHON) scripts/export_results.py
