localrules: filter, partition_sequences, aggregate_alignments, mask, adjust_metadata_regions, clades, colors, recency, export, incorporate_travel_history, fix_colorings, all_regions, export_all_regions, export_gisaid, incorporate_travel_history_gisaid, incorporate_travel_history_zh, export_zh, dated_json, fix_colorings_zh, fix_colorings_gisaid, finalize, pangolin, rename_legacy_clades, adjust_metadata_regions_ecdc, reassign_metadata

ruleorder: adjust_metadata_regions_ecdc > adjust_metadata_regions

rule reassign_metadata:
    input:
        metadata = rules.download.output.metadata
    output:
        metadata = "results/metadata_with_reassignments.tsv"
    run:
        import pandas as pd

        d = pd.read_csv(input.metadata, sep='\t')
        ind = d['country'] == 'Israel'
        d.loc[ind,'region'] = 'Europe'
        d.loc[ind,'region_exposure'] = 'Europe'
        d.to_csv(output.metadata, sep='\t')

rule adjust_metadata_regions_ecdc:
    message:
        """
        Adjusting metadata for build '{wildcards.build_name}'
        """
    input:
        metadata = rules.reassign_metadata.output.metadata
    output:
        metadata = "results/{build_name}/metadata_adjusted.tsv"
    params:
        region = lambda wildcards: config["builds"][wildcards.build_name]["region"]
    log:
        "logs/adjust_metadata_regions_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        {python:q} scripts/adjust_regional_meta.py \
            --region {params.region:q} \
            --metadata {input.metadata} \
            --output {output.metadata} 2>&1 | tee {log}
        """
