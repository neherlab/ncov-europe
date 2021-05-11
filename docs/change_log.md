# Change Log

As of April 2021, we use major version numbers (e.g. v2) to reflect backward incompatible changes to the workflow that likely require you to update your Nextstrain installation.
We also use this change log to document new features that maintain backward compatibility, indicating these features by the date they were added.

## v5 (7 May 2021)

[See the corresponding pull request](https://github.com/nextstrain/ncov/pull/615) for more details about this release.

### Major changes

- Drop support for old sequence/metadata inputs
- Use nextalign for alignment instead of mafft

### Minor changes

- Drop unused haplotype status rule and script
- Remove unused nucleotide mutation frequencies rule
- Use augur distance for mutation counts

## v4 (5 May 2021)

[See the corresponding pull request](https://github.com/nextstrain/ncov/pull/605) for more details about changes in this release.

### Major changes

- Change the default build name from "global" to "default-build" and use a default subsampling scheme that selects all input sequences
- Warn about duplicate sequences found when merging sequences from multiple inputs instead of throwing an error (set `combine_sequences_for_subsampling: warn_about_duplicates: false` in your configuration file to revert this behavior)

### Features

- Define a new subsampling scheme named `all` that selects all input sequences
- Add a new top-level configuration parameter `default_build_name` to allow overriding new default name of "default-build"
- Support compressed sequence inputs for alignment with mafft and nextalign (requires mafft upgrade)
- Sanitize strain names in sequences and metadata from different sources (e.g., `hCoV-19/` from GISAID or `SARS-CoV-2/` from GenBank, etc.)

## New features since last version update

- 20 April 2021: Surface emerging lineage as a colorby. This replaces the rather stale color by "Emerging Clade" with a new color by "Emerging Lineage". This focuses on PANGO lineages that are of interest triangulated by [CoVariants](https://covariants.org/), [PANGO](https://cov-lineages.org/) international lineage reports, [CDC](https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/variant-surveillance/variant-info.html) VUIs and VOCs and [PHE](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/975742/Variants_of_Concern_VOC_Technical_Briefing_8_England.pdf) VUIs and VOCs. The intention is for the listing at `emerging_lineages.tsv` to be updated frequently with new lineages added and no longer interesting lineages dropped. [#609](https://github.com/nextstrain/ncov/pull/609)

- 12 April 2021: Calculate current clade frequency and logistic growth rate across nodes in the phylogeny. This produces a new `logistic_growth.json` file and uses this file to add a coloring the final Auspice JSON. Implementation choices are discussed in PR [#595](https://github.com/nextstrain/ncov/pull/595).

- 12 April 2021: Annotate Pangolin lineages per build in a `pangolineages.json` file and final Auspice JSON by adding `run_pangolin: true` to the top-level workflow config (`builds.yaml`). **Note: this annotation only works when running the workflow with Snakemake's `--use-conda` flag or if your environment has [Pangolin installed](https://github.com/cov-lineages/pangolin).** [#593](https://github.com/nextstrain/ncov/pull/593)


## v3 (12 April 2021)

- Use Augur 11.2.0's metadata-only output arguments to aggregate subsampled sequences and metadata [#592][]
- Use Augur 11.3.0's `io.py` module to combine and deduplicate uncompressed or compressed sequences when handling multiple input datasets [#592][]

[#592]: https://github.com/nextstrain/ncov/pull/592

## v2 (9 April 2021)

This release reflects the state of the workflow when we instituted our workflow versioning system.

## July 2020 Update

In order to make our repository more accessible for first-time users, and those who want to create their own customized build, we made a number of changes in July 2020.
Here, we outline those changes, and also list 'breaking changes' that may effect advanced users who were already running their own builds.

You can check the lists below to see if there's any action you could take to make updating to the latest version of this repository is smooth!

### Summary of changes

The goal of this release is to make it easier for new users to get their own `ncov` build up and running. This consists of a new, extensive tutorial hosted on github pages; a simplified repository structure; and more didactic file names.

Summary of breaking changes:
- All nextstrain-specific files now live under `nextstrain_profiles`
- `config` is now named `defaults`
- `envs`, `rules` and `schemas` now live under `workflow`

#### Additions and updates
- README is updated
- Adds in-depth tutorial to `docs/` and accompanying github pages setup files
- `data/` now contains example data that should remain packaged with the repository
- `narratives/` now contains only a `template_narrative.md` for new users (not to worry -- all finished Situation Reports are still hosted over in the `nextstrain/narratives` repo)
- `my_profiles` now contains extensive `example` (formerly `king-county`) and `example_advanced_customization` (formerly `swiss`) profiles

#### File name changes
- The main `config.yaml` file is now more precisely labeled as `parameters.yaml`
- `reference.gb` > `reference_seq.gb`
- `ordering.tsv` > `color_ordering.tsv`
- `config` > `defaults`, with nextstrain-specific files (e.g., `description.md`) moved to `nextstrain_profiles/nextstrain/`
- nextstrain-specific profiles labeled with `nextstrain_*`


#### Structural changes
To improve repo clarity and approachability for new users, the following files have been reorganized.
- Workflow-related files now live in a `workflow` directory, rather than the top level
  - `envs` > `workflow/envs`
  - `schemas` > `workflow/schemas`
  - `rules` > `workflow/snakemake_rules`
- Nextstrain-specific profiles now live under `nextstrain_profiles`

(Reiterated from last section because they're important :)
- `config` is now `defaults`, with nextstrain-specific files (e.g., `description.md`) moved to `nextstrain_profiles/nextstrain/`
- `config.yaml` is now `parameters.yaml`

### Breaking Changes

If you've got custom builds running and regularly pull from the repository, you should check here to see if any of the changes you've made might make the merge harder.
Usually, just copying and renaming folders and/or files is enough to avoid too many merge conflicts!

First, you should check the file locations that have changed, above.
If you reference any of these files specifically in your pipelines or profile, be sure to update those links!

Be sure to also check any links to files that were previously in `profiles`.
For example, if you were referencing the cluster information in `profiles/nextstrain-scicore/cluster.json` for example) - we've now got two profiles folders - you'd find `nextstrain-scicore` in the `nextstrain_profiles` folder now.

If you modified the `swiss`, `default`, or `king-county` profiles, this may cause a merge conflict, as these folders have been renamed.
To keep your builds as they are, we recommend copying these folders, and then after the merge, restoring the name you prefer.

We now have a `data` folder which comes with the repository and contains example build data.
If you had previously created a `data` folder, this may cause a merge conflict - usually this can be easily resolved.


### Developer Changes

For those who help run and maintain the Nextstrain.org builds, there's a few changes that might be different from what you're used to.

Just focusing on the files mostly accessed by the 'build shepherds,' here's what's good to know:
- The `config` folder is now `defaults`
- `orderings.tsv` is now `color_orderings.tsv`
- `rules/nextstrain_exports.smk` (the rules that do the 'finalizing' of the builds) isnow `workflow/snakemake_rules/export_for_nextstrain.smk`

The 'finishing builds' commands are also slightly different because the profiles have moved - they'll look familiar though!
You can already find these changed commands in the `export_for_nextstrain.smk` file, but here they are as well:

```
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall --profile nextstrain_profiles/nextstrain
#   snakemake --profile nextstrain_profiles/nextstrain clean_export_regions
#   snakemake --profile nextstrain_profiles/nextstrain export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake --profile nextstrain_profiles/nextstrain all_regions
# to produce the final Auspice files!
```
