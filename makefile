# Makefile for deleting files with "chrM" or "chrY" in their names

all build_ext clean

all: build_ext

build_ext:
	python3 setup.py build_ext --inplace

clean:
	@echo "Deleting files containing 'chrM' or 'chrY' in their names..."
	find . -type f \( -name "*chrM*" -o -name "*chrY*" \) -exec rm -v {} \; > /dev/null 2>&1
	@echo "Cleanup complete!"
