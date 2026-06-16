PYTHON ?= python3
HF_REPO ?= jang1563/genelab-benchmark
HF_UPLOAD_PLAN ?= output/hf_upload_plan.json
HF_TASK ?= all

.PHONY: help
help:
	@echo "SpaceBio-Bench release helpers"
	@echo ""
	@echo "  make release-qa             Run local public release validation"
	@echo "  make hpc-public-qa          Run the public-surface slice of the HPC gate"
	@echo "  make public-docs-qa         Validate README, HF card, citation, and manifest"
	@echo "  make hf-upload-plan         Write a dry-run HF upload plan JSON"
	@echo "  make hf-card-plan           Dry-run only the HF dataset card upload"
	@echo "  make hf-manifest-plan       Dry-run only HF release metadata upload"
	@echo "  make validate-hf-remote     Post-upload smoke-check remote HF files"
	@echo ""
	@echo "Variables:"
	@echo "  PYTHON=$(PYTHON)"
	@echo "  HF_REPO=$(HF_REPO)"
	@echo "  HF_TASK=$(HF_TASK)"
	@echo "  HF_UPLOAD_PLAN=$(HF_UPLOAD_PLAN)"

.PHONY: public-docs-qa
public-docs-qa:
	$(PYTHON) scripts/validate_release_manifest.py
	$(PYTHON) scripts/validate_public_docs_consistency.py
	$(PYTHON) -m unittest tests/test_release_manifest.py tests/test_public_docs_consistency.py

.PHONY: release-qa
release-qa: public-docs-qa
	$(PYTHON) -m unittest tests/test_hf_upload_workflow.py
	$(PYTHON) -m py_compile scripts/upload_to_hf.py scripts/validate_release_manifest.py scripts/validate_public_docs_consistency.py

.PHONY: hpc-public-qa
hpc-public-qa:
	PYTHON_BIN=$(PYTHON) bash scripts/hpc_release_validate.sh --public-only

.PHONY: hf-upload-plan
hf-upload-plan:
	$(PYTHON) scripts/upload_to_hf.py --task $(HF_TASK) --repo $(HF_REPO) --dry-run --write-upload-plan $(HF_UPLOAD_PLAN)

.PHONY: hf-card-plan
hf-card-plan:
	$(PYTHON) scripts/upload_to_hf.py --repo $(HF_REPO) --card-only --dry-run

.PHONY: hf-manifest-plan
hf-manifest-plan:
	$(PYTHON) scripts/upload_to_hf.py --repo $(HF_REPO) --manifest-only --dry-run

.PHONY: validate-hf-remote
validate-hf-remote:
	$(PYTHON) scripts/upload_to_hf.py --task $(HF_TASK) --repo $(HF_REPO) --dry-run --remote-diff --validate-remote
