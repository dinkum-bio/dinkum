.PHONY: dist

all:
	./simple.py
	./simple2.py
	./simple3.py

test:
	pytest

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

dist:
	python -m build
