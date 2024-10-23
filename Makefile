ver= $(shell grep -i version pyproject.toml | cut -z -f2 -d\')

default:
	python -m build
clean:
	rm dist/*

install:
	pip install --upgrade dist/pmike-$(ver)-*whl
remove:
	pip uninstall pmike
