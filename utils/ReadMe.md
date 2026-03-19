## Analysis Scenarios

The most important choice for running scTIGER2.0 is whether or not to use the time-dependence (causal) validation step. In practice, the choice depends on whether scTIGER2.0 is being used as a **discovery tool** or a **decision-support tool**.

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

```
./run_scTIGER2.py \
  -goi Gene1+Gene2+Gene3 \
  -exp ./data.csv \
  -nV \
  -top 100 \
  -o exploratory_output
```

Here, causal validation is turned off and we increase the number of top genes from 50 to 100. Increasing the number of top genes increases run time, but casts a wider net for identification of potential causal edges. 

**Expected Outcome**
Results priotitize recall of regulatory edges, but not necessarily high confidence edges. 

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

```
./run_scTIGER2.py \
  -goi Gene1+Gene2+Gene3 \
  -exp ./data.csv \
  -nV \
  -top 100 \
  -a 0.1 \
  -p 50 \
  -zero 0.2 \
  -o pilot_output
```

Here causal validation is turned off. We change several other flags as exploration may benefit from relaxed parameters. The number of top genes in increased from 50 to 100 to cast a wider net of potential causal edges. The alpha value is increased from 0.05 to 0.10 to filter out fewer edges in negative binomial fitting. The permutation number is decreased from 100 to 50 as high confidence edges are not required to explore the dataset. The zero threshold is also decreased from 0.3 to 0.2 to filter out fewer genes before correlation testing. 


**Expected Outcome**
Results prioritize broad edge detection with relaxed filtering, enabling rapid exploration of candidate regulatory relationships. Designed to be a starting place for you to refine the parameters of scTIGER2.0 to fir your analysis needs. 

---

### C. Designing a Perturbation Experiment  
**Causal Validation: ON**

Use when results will directly guide experimental intervention.

**Use this when:**
- Selecting genes for CRISPR knockout/knockdown  
- Designing perturb-seq experiments  
- Experimental resources are limited  
- False positives are costly  

```
./run_scTIGER2.py \
  -goi Gene1+Gene2+Gene3 \
  -exp ./data.csv \
  -top 100 \
  -a 0.01 \
  -p 200 \
  -o perturbation_output
```

Here we increase the number of top genes from 50 to 100 to cast a wider net and potentially catch off-target effects. The alpha value is decreased from 0.05 to 0.01 increase stringency and only keep high confidence edges. The number of permutations are therefore raised from 100 to 200 to aid negative binomial fitting. 

**Expected Outcome**
Results priotitize specificity of identify high confidence targets for your experimental work. 

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

```
./run_scTIGER2.py \
  -goi Gene1+Gene2+Gene3 \
  -exp ./data.csv \
  -o mechanism_output
```

**Expected Outcome**
Results priotitize specificity of idetifying relationships between your candidiate genes. 

---

