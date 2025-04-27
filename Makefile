.PHONY: dist nbconvert

all: test

format:
	black src

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

nbconvert:
	cd notebooks && \
	   rm -f *.nbconvert.ipynb && \
	   jupyter nbconvert --to notebook --execute [0-9]-*.ipynb --allow-errors

update_notebooks: 
	cd notebooks && for i in *.nbconvert.ipynb; do mv $$i $$(basename $$i .nbconvert.ipynb).ipynb; done

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

dist:
	python -m build
