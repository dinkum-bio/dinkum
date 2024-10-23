.PHONY: dist

all:
	./simple.py
	./simple2.py
	./simple3.py

test:
	pytest

nbtest:
	py.test --nbval-lax notebooks/getting-started.ipynb
	py.test --nbval-lax notebooks/1-positive-feedback.ipynb
	py.test --nbval-lax notebooks/2-simple-oscillation.ipynb
	py.test --nbval-lax notebooks/4-double-negative-gate.ipynb
	py.test --nbval-lax notebooks/5-intermediate-custom-logic.ipynb
	py.test --nbval-lax notebooks/6-decay-example.ipynb
	py.test --nbval-lax notebooks/9-advanced-examples.ipynb

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

dist:
	python -m build
