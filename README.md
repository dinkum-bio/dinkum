# dinkum

Directed Interaction NetworKs are fair dinkum!

## What is dinkum?

`dinkum` is a piece of software for simple modeling of gene regulatory
networks, initially based on the GeNeTool software described in
[Faure, Peter, and Davidson, 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3667423/). It
supports simple definition of genes, tissues, ligands/receptors, and
regulatory influences, and is intended to run in Jupyter
Notebooks. It's primarily intended for teaching purposes, and was
developed for the 2024
[Gene Regulatory Networks for Development](https://www.mbl.edu/education/advanced-research-training-courses/course-offerings/gene-regulatory-networks-development)
course at the Marine Biological Laboratory.

To get started with dinkum, see [notebooks/0-getting-started.ipynb](notebooks/0-getting-started.ipynb). Here is a full list of example Jupyter notebooks:

* [0-getting-started.ipynb](notebooks/0-getting-started.ipynb)
* [1-positive-feedback.ipynb](notebooks/1-positive-feedback.ipynb)
* [2-simple-oscillation.ipynb](notebooks/2-simple-oscillation.ipynb)
* [4-double-negative-gate.ipynb](notebooks/4-double-negative-gate.ipynb)
* [5-intermediate-custom-logic.ipynb](notebooks/5-intermediate-custom-logic.ipynb)
* [6-decay-example.ipynb](notebooks/6-decay-example.ipynb)
* [6-multi-level-activation.ipynb](notebooks/6-multi-level-activation.ipynb)
* [7-fit-functions.ipynb](notebooks/7-fit-functions.ipynb)
* [9-advanced-examples.ipynb](notebooks/9-advanced-examples.ipynb)

## Why 'dinkum'?

Dinkum is a backronym constructed from "directed interaction
networks". It's also named in honor of one of the course directors,
who is Australian; it turns out that dinkum is one of the few
Australian-specific slang words that is not rude.

## Installing dinkum

dinkum is available on the Python Package Index, PyPI, as `dinkum-bio`.
It requires Python 3.11 or later.

To install:
```
pip install dinkum-bio
```

## Developing dinkum

Dinkum is developed on github under
[dinkum-bio/dinkum](https://github.com/dinkum-bio/dinkum/). It is
released under the GNU Affero General Public License v3 open source
license.

You can run the tests with `make test`.

CTB 10/2024
