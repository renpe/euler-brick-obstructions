# Reproduction Makefile for Paper 4
#
# Verification targets that run the analysis scripts against the
# published `pub` schema in PostgreSQL and report on the empirical
# claims of the paper.
#
# Usage:
#   make verify-theorem      # Theorem 4.2 (verified exponent-one blocker)
#   make verify-perfect      # Theorem 4.9 (no perfect cuboid in D_all)
#   make verify-consistency  # Theorem 8.1 (internal database consistency)
#   make verify-all          # all of the above

PYTHON ?= python3

.PHONY: verify-theorem verify-perfect verify-consistency verify-all

verify-theorem:
	$(PYTHON) analysis/theorem_check.py

verify-perfect:
	$(PYTHON) analysis/perfect_check.py

verify-consistency:
	$(PYTHON) analysis/validate_xyz.py

verify-all: verify-theorem verify-perfect verify-consistency
