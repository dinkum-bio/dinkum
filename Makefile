.PHONY: dist nbconvert book

all: test

format:
	black src

book:
	jupyter-book build book

test:
	pytest --cov --cov-report=html:coverage_report

nbtest:
	py.test --nbval-lax --current-env notebooks/0-getting-started.ipynb \
	    notebooks/1-positive-feedback.ipynb \
	    notebooks/2-simple-oscillation.ipynb \
        notebooks/4-double-negative-gate.ipynb \
	    notebooks/5-intermediate-custom-logic.ipynb \
	    notebooks/6-decay-example.ipynb \
	    notebooks/6-multi-level-activation.ipynb \
	    notebooks/7-fit-functions.ipynb \
	    notebooks/9-advanced-examples.ipynb

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
