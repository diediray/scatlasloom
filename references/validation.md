# Validation scope

This release includes reusable R preparation, rendering and scoring tools, plus an agent workflow requiring explicit user review of clustering and annotation before downstream analysis or formal figure production. Visual acceptance and scientific acceptance are recorded separately.

The tools have been exercised on simulated layouts and local single-cell case studies. Local datasets, patient-level tables, private feedback, and unpublished biological results are not included. These checks establish operation in the tested cases, not universal annotation accuracy or journal acceptance.

The compact renderer currently supports 2–9 clusters and up to six selected marker rows per cluster. Larger hierarchies use the separate `render_atlas.R` interface or require a reviewed implementation adaptation. It is not an automatic clustering or annotation pipeline.

The synthetic example in `examples/` is a rendering smoke test. Its coordinates and values are generated, its features are artificial, and it provides no biological evidence.

For real data, review cell membership, lineage labels, marker support, contradictory signals, sample support and QC with the user. Preserve uncertainty. After scientific inputs change, review the affected version again.

The hierarchy gallery is a separate synthetic example with 80,000 generated coordinates, 28 subtypes, five parent groups, six marker features, six program features, and eight signature features. All memberships and values are artificial; this example validates rendering and alignment only.
