## 🌐 Live Interactive Report

[Open the complete AI Genomics Platform report](https://educationus3003-afk.github.io/AI-Genomics-Platform/)# AI Genomics Platform

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Research](https://img.shields.io/badge/status-research%20prototype-orange)
![License](https://img.shields.io/badge/license-MIT-green)
![Genomics](https://img.shields.io/badge/domain-genomic%20variant%20interpretation-purple)

An integrated research platform for genomic variant annotation, database
evidence integration, ACMG/AMP evidence assessment, AI-assisted explanation,
visualization, and automated HTML reporting.

> **Research-use only:** This software is not validated for clinical diagnosis,
> treatment selection, genetic counselling, or patient management.

---

## Overview

AI Genomics Platform accepts genomic variant data and organizes evidence from
annotation tools, population resources, clinical databases, disease knowledge,
and prediction resources into transparent research reports.

The demonstration dataset currently covers:

| Gene | Associated condition | Inheritance |
|---|---|---|
| `CFTR` | Cystic fibrosis | Autosomal recessive |
| `HBB` | Beta-thalassemia / sickle-cell disorders | Autosomal recessive |
| `MECP2` | Rett syndrome | X-linked dominant |
| `FBN1` | Marfan syndrome | Autosomal dominant |

## Main capabilities

- VCF normalization and validation
- Ensembl VEP and SnpEff annotation
- ClinVar evidence retrieval
- gnomAD population-frequency lookup
- UniProt protein annotation
- OMIM-compatible disease integration
- ClinGen gene-disease and dosage resources
- AlphaMissense, REVEL and CADD integration structure
- Research-oriented ACMG/AMP evidence assessment
- AI-assisted evidence summaries
- SNV/indel and CNV workflow documentation
- Disease-specific figures and HTML reports

## Workflow

```text
VCF
 │
 ├── Validation and normalization
 │
 ├── VEP / SnpEff functional annotation
 │
 ├── Database integration
 │   ├── ClinVar
 │   ├── gnomAD
 │   ├── ClinGen
 │   ├── OMIM
 │   ├── UniProt
 │   └── Prediction resources
 │
 ├── ACMG/AMP evidence assessment
 │
 ├── AI-assisted explanation
 │
 └── Tables, figures and HTML report

## Live HTML Report

View the complete interactive report through GitHub Pages:

https://educationus3003-afk.github.io/AI-Genomics-Platform/
