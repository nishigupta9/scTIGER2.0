## Analysis Scenarios

### 1. Inclusion of the Time-Dependence (Causal) Validation Step

In practice, the choice depends on whether scTIGER2.0 is being used as a **discovery tool** or a **decision-support tool**.

- If your goal is **hypothesis generation or exploratory screening**, especially with downstream computational filtering, the causal validation step can be turned **off** (`-nV` flag).
- If your goal is to **perturb a target experimentally and interpret the outcome**, the step should be left **on**.

---

## Example Scenarios

### A. Exploratory Target Discovery from Single-Cell RNA-seq Data  
**Causal Validation: OFF**

This is the broadest discovery setting. The goal is to identify genes, pathways, or regulators associated with a phenotype without making strong mechanistic claims.

**Use this when:**
- Analyzing a new dataset  
- Biology is poorly characterized  
- Sensitivity matters more than specificity  
- You plan downstream GRN inference or computational filtering  

**Example:**  
Run scTIGER on diseased vs. healthy cell states to identify candidate regulators, then evaluate hits using in silico approaches such as protein binding, structural modeling, pathway enrichment, or network prioritization. Here, scTIGER acts as an upstream screening tool.

---

### B. Early-Stage Feasibility Studies or Pilot Analyses  
**Causal Validation: OFF**

Used to assess whether the method produces biologically meaningful signals before applying stricter criteria.

**Use this when:**
- Dataset quality is uncertain  
- Benchmarking preprocessing choices  
- Evaluating parameter sensitivity  
- Work is still exploratory  

**Example:**  
Compare outputs across gene sets, filtering thresholds, or cell populations to check whether known regulators emerge.

---

### C. Designing a Perturbation Experiment  
**Causal Validation: ON**

Use when results will directly guide experimental intervention.

**Use this when:**
- Selecting genes for CRISPR knockout/knockdown  
- Designing perturb-seq experiments  
- Experimental resources are limited  
- False positives are costly  

---

### D. Testing a Specific Mechanistic Hypothesis  
**Causal Validation: ON**

Used to evaluate whether a candidate regulator plausibly drives a phenotype.

**Use this when:**
- Validating literature-driven hypotheses  
- Working with a small candidate set  
- Mechanistic interpretability is important  
- Making directional biological claims  

**Example:**  
Test whether a transcription factor implicated in prior studies alters a differentiation trajectory when perturbed.

---

