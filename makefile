# Make will use bash instead of sh
SHELL := /usr/bin/env bash

.PHONY: help
help:
	@echo ' '
	@echo '    make build   	Builds the code base incrementally (fast) for dev.'
	@echo '    make check   	Checks the code base for security vulnerabilities.'
	@echo '    make fix   		Fixes linting issues as reported by clippy.'
	@echo '    make format   	Formats call code according to cargo fmt style.'
	@echo '    make test   	Runs all tests across all crates.'

# "---------------------------------------------------------"
# "---------------------------------------------------------"

.PHONY: build
build:
	@source build/scripts/build.sh

.PHONY: check
check:
	@source build/scripts/check.sh


.PHONY: fix
fix:
	@source build/scripts/fix.sh


.PHONY: format
format:
	@source build/scripts/format.sh


.PHONY: release
release:
	@source build/scripts/release.sh

.PHONY: test
test:
	@source build/scripts/test.sh

.PHONY: sbom
sbom:
	 @source build/scripts/sbom.sh

.PHONY: vendor
vendor:
	@source build/scripts/vendor.sh