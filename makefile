# Makefile for deleting files with "chrM" or "chrY" in their names

.PHONY: clean

clean:
	@echo "Deleting files containing 'chrM' or 'chrY' in their names..."
	find . -type f \( -name "*chrM*" -o -name "*chrY*" \) -exec rm -v {} \; > /dev/null 2>&1
	@echo "Cleanup complete!"
