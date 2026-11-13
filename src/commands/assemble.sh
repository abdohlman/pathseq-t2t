cmd_assemble() {
  local sample_id="" input_unaligned="" input_decoys=""
  local threads="" dont_overwrite=0 keep_intermediate=0
  local assembly_dir="" min_contig_len="2500"
  local trim_galore_args="--illumina --stringency 5 --length 75 --max_n 2 --trim-n"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sample-id)           sample_id="${2:-}"; shift 2 ;;
      --outdir)              OUTDIR="${2:?--outdir requires a directory path}"; shift 2; _set_outdirs ;;
      --assembly-dir)        assembly_dir="${2:-}"; shift 2 ;;
      --input-unaligned)     input_unaligned="${2:-}"; shift 2 ;;
      --input-decoys)        input_decoys="${2:-}"; shift 2 ;;
      --threads)             threads="${2:-}"; shift 2 ;;
      --min-contig-len)      min_contig_len="${2:-}"; shift 2 ;;
      --trim-galore-args)    trim_galore_args="${2:-}"; shift 2 ;;
      --dont-overwrite)      dont_overwrite=1; shift ;;
      --keep-intermediate)   keep_intermediate=1; shift ;;
      -h|--help)
        cat <<'HLP'
Usage: pathseq-t2t assemble \
  --input-unaligned <bam> \
  [--input-decoys <bam>] \
  [--sample-id <id>] \
  [--outdir <dir>] \
  [--assembly-dir <dir>] \
  [--threads <int>] \
  [--min-contig-len <int>] \
  [--trim-galore-args "<args>"] \
  [--dont-overwrite] \
  [--keep-intermediate]

What it does:
  (1) BAM -> FASTQ (streams unaligned + decoys; preserves pairing via samtools collate)
      Note: singleton/unpaired reads are discarded.
  (2) trim_galore (paired only; default args mirror manuscript)
  (3) MEGAHIT assembly
  (4) bowtie2-build + bowtie2 realignment to contigs
  (5) samtools sort/index mapped BAM
  (6) MetaBAT2 binning (min contig length default 2500bp)

Outputs (default):
  <outdir>/assembly/<sample>/*

Notes:
  - If --sample-id is provided, --assembly-dir defaults to <outdir>/assembly/<sample-id>.
  - If --sample-id is omitted, --assembly-dir is required.
HLP
        return 0 ;;
      --) shift; break ;;
      -*) die "Unknown option for assemble: $1" ;;
      *)  die "Unexpected argument to assemble: $1" ;;
    esac
  done

  _require_metabat2
  local -a trim_galore_args_ary=()
  _split_cli_args "${trim_galore_args}" trim_galore_args_ary

  # Threads auto-detect
  if [[ -z "${threads}" ]]; then
    if command -v nproc >/dev/null 2>&1; then threads="$(nproc)"
    elif [[ "$OSTYPE" == "darwin"* ]] && command -v sysctl >/dev/null 2>&1; then threads="$(sysctl -n hw.ncpu)"
    else threads="8"; fi
  fi

  # Naming policy
  local base=""
  if [[ -n "${sample_id}" ]]; then
    base="${sample_id}"
    [[ -n "${assembly_dir}" ]] || assembly_dir="${OUTDIR%/}/assembly/${base}"
  else
    require_nonempty "${assembly_dir}" "--assembly-dir (required when --sample-id is not provided)"
    base="$(basename "${assembly_dir%/}")"
  fi

  # Outdir
  mkdir -p "${assembly_dir}"

  # Layout
  local fqdir="${assembly_dir}/fastq"
  local tgdir="${assembly_dir}/trim_galore"
  local mhdir="${assembly_dir}/megahit"
  local bt2dir="${assembly_dir}/bowtie2"
  local bindir="${assembly_dir}/metabat2_bins"
  mkdir -p "${bt2dir}" "${bindir}"

  # Output paths
  local fq1="${fqdir}/${base}.R1.fq.gz"
  local fq2="${fqdir}/${base}.R2.fq.gz"
  local trim1="${tgdir}/${base}.pe_val_1.fq.gz"
  local trim2="${tgdir}/${base}.pe_val_2.fq.gz"
  local mh_contigs="${mhdir}/${base}.contigs.fa"
  local bt2_base="${bt2dir}/${base}.contigs"
  local map_bam="${bt2dir}/${base}.contigs.bowtie2.sorted.bam"
  local depth_txt="${bt2dir}/${base}.contigs.depth.txt"

  # Pre-flight: check which outputs already exist
  local have_fq=0 have_trim=0 have_contigs=0 have_bam=0 have_depth=0
  [[ -s "${fq1}" && -s "${fq2}" ]] && have_fq=1
  [[ -s "${trim1}" && -s "${trim2}" ]] && have_trim=1
  { [[ -s "${mh_contigs}" ]] || [[ -s "${mh_contigs}.gz" ]]; } && have_contigs=1
  [[ -s "${map_bam}" && -s "${map_bam}.bai" ]] && have_bam=1
  [[ -s "${depth_txt}" ]] && have_depth=1

  # Staleness guard: if contigs are newer than BAM or depth, invalidate them
  # so downstream steps regenerate from the current assembly.
  if (( dont_overwrite && have_contigs )); then
    local _contigs_ref="${mh_contigs}"
    [[ -s "${_contigs_ref}" ]] || _contigs_ref="${mh_contigs}.gz"
    if (( have_bam )) && [[ "${_contigs_ref}" -nt "${map_bam}" ]]; then
      log "assemble --dont-overwrite: contigs newer than BAM — will regenerate bowtie2 + depth"
      have_bam=0 have_depth=0
    elif (( have_depth )) && [[ "${_contigs_ref}" -nt "${depth_txt}" ]]; then
      log "assemble --dont-overwrite: contigs newer than depth file — will regenerate depth"
      have_depth=0
    fi
  fi

  # Validate inputs if early steps need to run.
  # Steps 1-2 run when: contigs missing, OR contigs exist but BAM needs
  # regenerating and trimmed reads were cleaned up.
  local _need_early_steps=0
  if (( ! (dont_overwrite && have_contigs) )); then
    _need_early_steps=1
  elif (( dont_overwrite && have_contigs && !have_bam && !have_trim && !have_fq )); then
    _need_early_steps=1
  fi
  if (( _need_early_steps )); then
    require_nonempty "${input_unaligned}" "--input-unaligned"
    require_file "${input_unaligned}"
    ubam_check_or_die "${input_unaligned}" "assemble: input unaligned"
    if [[ -n "${input_decoys}" ]]; then
      require_file "${input_decoys}"
      ubam_check_or_die "${input_decoys}" "assemble: input decoys"
    fi
    _require_samtools_116
    _require_trim_galore_0610
    if (( !have_contigs )); then
      _require_megahit_129
    fi
    _require_bowtie2
  fi

  # ---------- (1) BAM -> FASTQ (stream cat|collate|fastq) ----------
  if (( dont_overwrite && (have_fq || have_trim || (have_contigs && have_bam)) )); then
    log "assemble --dont-overwrite: skipping BAM -> FASTQ"
  else
    mkdir -p "${fqdir}"
    log "assemble: BAM -> FASTQ (threads=${threads})"
    if [[ -n "${input_decoys}" ]]; then
      log_cmd samtools merge -@ "${threads}" -u -o - "${input_unaligned}" "${input_decoys}" \
        | log_cmd samtools collate -@ "${threads}" -u -O - \
        | log_cmd samtools fastq -@ "${threads}" -n -F 0x900 - \
            -1 >( _compress_to "${fq1}" "${threads}" ) \
            -2 >( _compress_to "${fq2}" "${threads}" ) \
            -s /dev/null \
            -0 /dev/null
    else
      log_cmd samtools collate -@ "${threads}" -u -O "${input_unaligned}" \
        | log_cmd samtools fastq -@ "${threads}" -n -F 0x900 - \
            -1 >( _compress_to "${fq1}" "${threads}" ) \
            -2 >( _compress_to "${fq2}" "${threads}" ) \
            -s /dev/null \
            -0 /dev/null
    fi
    [[ -s "${fq1}" && -s "${fq2}" ]] || die "FASTQ extraction failed (missing/empty paired FASTQs)."
  fi

  # ---------- (2) trim_galore ----------
  if (( dont_overwrite && (have_trim || (have_contigs && have_bam)) )); then
    log "assemble --dont-overwrite: skipping trim_galore"
  else
    mkdir -p "${tgdir}"
    log "assemble: trim_galore"
    log_cmd trim_galore --paired --gzip -j "${threads}" --basename "${base}.pe" -o "${tgdir}" \
      "${trim_galore_args_ary[@]}" "${fq1}" "${fq2}"
    [[ -s "${trim1}" && -s "${trim2}" ]] || die "trim_galore paired outputs missing: ${trim1}, ${trim2}"
  fi

  # ---------- (3) MEGAHIT ----------
  if (( dont_overwrite && have_contigs )); then
    log "assemble --dont-overwrite: skipping MEGAHIT (found contigs)"
  else
    log "assemble: MEGAHIT"
    if [[ -d "${mhdir}" ]]; then
      rm -rf "${mhdir}"
    fi
    log_cmd megahit -1 "${trim1}" -2 "${trim2}" -o "${mhdir}" --out-prefix "${base}" -t "${threads}"
    [[ -s "${mh_contigs}" ]] || die "MEGAHIT contigs missing: ${mh_contigs}"
  fi

  # Decompress contigs if only gzipped copy exists (sbatch wrapper may have gzipped them)
  if [[ ! -s "${mh_contigs}" && -s "${mh_contigs}.gz" ]]; then
    gunzip -c "${mh_contigs}.gz" > "${mh_contigs}"
  fi

  # ---------- (4) bowtie2-build + realignment ----------
  if (( dont_overwrite && have_bam )); then
    log "assemble --dont-overwrite: skipping bowtie2-build + align (found ${map_bam})"
  else
    log "assemble: Bowtie2 index + align"
    local -a bt2_idx_files=(
      "${bt2_base}.1.bt2" "${bt2_base}.2.bt2" "${bt2_base}.3.bt2"
      "${bt2_base}.4.bt2" "${bt2_base}.rev.1.bt2" "${bt2_base}.rev.2.bt2"
    )
    local bt2_idx_ok=1
    local idxf
    for idxf in "${bt2_idx_files[@]}"; do
      [[ -s "${idxf}" ]] || bt2_idx_ok=0
    done
    if (( dont_overwrite && bt2_idx_ok )); then
      log "assemble --dont-overwrite: skipping bowtie2-build (index files found)"
    else
      log_cmd bowtie2-build "${mh_contigs}" "${bt2_base}"
    fi

    log_cmd bowtie2 -x "${bt2_base}" -1 "${trim1}" -2 "${trim2}" -p "${threads}" \
      | log_cmd samtools view -@ "${threads}" -b - \
      | log_cmd samtools sort -@ "${threads}" -o "${map_bam}"
    log_cmd samtools index -@ "${threads}" "${map_bam}"
  fi
  bam_check_or_die "${map_bam}" "assemble: contig-mapped BAM"

  # ---------- (5) Depth estimation ----------
  if (( dont_overwrite && have_depth )); then
    log "assemble --dont-overwrite: skipping depth estimation (found ${depth_txt})"
  else
    log "assemble: depth estimation"
    log_cmd jgi_summarize_bam_contig_depths --outputDepth "${depth_txt}" "${map_bam}"
  fi
  [[ -s "${depth_txt}" ]] || die "Depth file missing: ${depth_txt}"

  # ---------- (6) MetaBAT2 ----------
  log "assemble: MetaBAT2"
  mkdir -p "${bindir}"
  if (( dont_overwrite )) && compgen -G "${bindir}/bin*.fa" >/dev/null \
     && [[ ! "${depth_txt}" -nt "$(ls -t "${bindir}"/bin*.fa | tail -1)" ]]; then
    log "assemble --dont-overwrite: skipping MetaBAT2 binning (found existing bin FASTAs)"
  else
    log_cmd metabat2 -i "${mh_contigs}" -a "${depth_txt}" -o "${bindir}/bin" -m "${min_contig_len}" -t "${threads}"
  fi

  # Write metabat2 summary (completion marker + bin stats)
  local metabat2_summary="${bindir}/metabat2_summary.tsv"
  if compgen -G "${bindir}/bin*.fa" >/dev/null; then
    local -a _bin_fas=("${bindir}"/bin*.fa)
    local _num_bins="${#_bin_fas[@]}"
    awk -v nb="${_num_bins}" '
      /^>/ { contigs++; next }
      { len += length($0) }
      END { printf "num_bins\ttotal_contigs\ttotal_length\n%d\t%d\t%d\n", nb, contigs+0, len+0 }
    ' "${_bin_fas[@]}" > "${metabat2_summary}"
  else
    printf 'num_bins\ttotal_contigs\ttotal_length\n0\t0\t0\n' > "${metabat2_summary}"
  fi
  log "assemble: metabat2 summary -> ${metabat2_summary}"

  # Cleanup
  if (( ! keep_intermediate )); then
    # Remove transient FASTQ + trim outputs
    rm -f "${fq1}" "${fq2}" "${trim1}" "${trim2}" 2>/dev/null || true

    # Bowtie2 index files are rebuildable and not needed after mapped BAM exists.
    rm -f "${bt2_base}".*.bt2 "${bt2_base}".*.bt2l "${bt2_base}".rev.*.bt2 "${bt2_base}".rev.*.bt2l 2>/dev/null || true

    # MEGAHIT intermediates are large and not required for downstream steps.
    rm -rf "${mhdir}/intermediate_contigs" 2>/dev/null || true
    rm -f "${mhdir}/checkpoints.txt" "${mhdir}/done" 2>/dev/null || true

    # Remove empty transient dirs for a cleaner assembly layout.
    rmdir "${fqdir}" 2>/dev/null || true
    rmdir "${tgdir}" 2>/dev/null || true
  fi

  log "assemble done → ${assembly_dir}"
}
