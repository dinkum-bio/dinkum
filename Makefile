.PHONY: dist

all:
	./simple.py
	./simple2.py
	./simple3.py

test:
	pytest --cov --cov-report=html:coverage_report

nbtest:
	py.test --nbval-lax notebooks/0-getting-started.ipynb --current-env
	py.test --nbval-lax notebooks/1-positive-feedback.ipynb --current-env
	py.test --nbval-lax notebooks/2-simple-oscillation.ipynb --current-env
	py.test --nbval-lax notebooks/4-double-negative-gate.ipynb --current-env
	py.test --nbval-lax notebooks/5-intermediate-custom-logic.ipynb --current-env
	py.test --nbval-lax notebooks/6-decay-example.ipynb --current-env
	py.test --nbval-lax notebooks/6-multi-level-activation.ipynb --current-env
	py.test --nbval-lax notebooks/9-advanced-examples.ipynb --current-env

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

dist:
	python -m build
