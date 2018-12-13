#!/usr/bin/env bash

# Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

Rscript /opt/code_docker_lib/merge_gxg_gtex_variants.R "$@"

