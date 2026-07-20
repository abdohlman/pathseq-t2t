cmd_qcfilter() {
  # ---------- Defaults & arg parsing ----------
  local input_unaligned="" input_decoys="" paired_out="" unpaired_out=""
  local filter_metrics_unaligned="" filter_metrics_decoys=""
  local unaln_final_paired=0 unaln_final_unpaired=0
  local decoys_final_paired="" decoys_final_unpaired=""
  local ram_gb="16"
  local threads=""                       # auto-detect below if not provided
  local TMPDIR_OPT="${TMPDIR:-/tmp}"     # from environment if set
  local HOSTDIR_OPT="${HOSTDIR:-}"       # from environment if set (required)
  local min_clipped_read_length="60"
  local sample_id=""                     # optional override
  local dont_overwrite=0
  local keep_intermediate=0
  local PICARD_JAR="${PICARD_JAR:-}"     # optional; else use 'picard' on PATH
  local psfilterspark_args=""            # extra args to PathSeqFilterSpark

  # 1) OUTDIR defaults up front (user can still override via --outdir)
  : "${OUTDIR:=./pst2t_out}"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      # IDs / dirs to use as defaults for input/output
      --sample-id)           sample_id="${2:-}";       shift 2 ;;
      --outdir)              OUTDIR="${2:?--outdir requires a directory path}" ; shift 2 ; _set_outdirs ;;

      # Inputs (explicit; otherwise inferred from --outdir + --sample-id/base)
      --input-unaligned)     input_unaligned="${2:-}"; shift 2 ;;
      --input-decoys)        input_decoys="${2:-}";    shift 2 ;;

      # Metrics (explicit; otherwise inferred from --outdir + basenames)
      --metrics-unaligned)   filter_metrics_unaligned="${2:-}"; shift 2 ;;
      --metrics-decoys)      filter_metrics_decoys="${2:-}";    shift 2 ;;

      # Outputs (explicit; otherwise inferred from --outdir + sample/base)
      --paired-out)          paired_out="${2:-}";      shift 2 ;;
      --unpaired-out)        unpaired_out="${2:-}";    shift 2 ;;

      # Reference
      --hostdir)             HOSTDIR_OPT="${2:-}";     shift 2 ;;

      # Resources / behavior
      --threads)             threads="${2:-}";         shift 2 ;;
      --ram-gb)              ram_gb="${2:-}";          shift 2 ;;
      --tmpdir)              TMPDIR_OPT="${2:-}";      shift 2 ;;
      --min-clipped-read-length) min_clipped_read_length="${2:-}"; shift 2 ;;
      --psfilterspark-args)  psfilterspark_args="${2:-}"; shift 2 ;;
      --picard-jar)          PICARD_JAR="${2:-}";      shift 2 ;;

      # Behavior toggles
      --dont-overwrite)      dont_overwrite=1;         shift ;;
      --keep-intermediate)   keep_intermediate=1;      shift ;;

      -h|--help)
        cat <<'HLP'
Usage: pathseq-t2t qcfilter [ARGS] \

Core options:
  [--outdir <dir>]                  Root output directory (default: ./pst2t_out)
  [--sample-id <string>]            Sample ID to use for naming outputs if needed
  [--hostdir <dir>]                 PathSeq host resource dir (requires pathseq_host.bfi + pathseq_host.fa.img)
  [--input-unaligned <bam>]         Unaligned BAM from prefilter (if omitted, inferred when --sample-id is provided)
  [--input-decoys  <bam> ]          Decoys BAM from prefilter (optional; if omitted, inferred when --sample-id is provided)
  [--paired-out <bam>]              Output merged paired reads (default: <outdir>/bams/<ID>.qcfilt_paired.bam)
  [--unpaired-out <bam>]            Output merged unpaired reads (default: <outdir>/bams/<ID>.qcfilt_unpaired.bam)
  [--metrics-unaligned <tsv>]       Explicit filter metrics path for unaligned (default: inferred)
  [--metrics-decoys <tsv>]          Explicit filter metrics path for decoys (default: inferred)

Performance / env:
  [--threads <int>]                 CPU threads (default: auto-detect)
  [--ram-gb <int>]                  Java heap for GATK, in GB (default: 16)
  [--tmpdir <dir>]                  Temp directory (default: $TMPDIR or /tmp)
  [--min-clipped-read-length <int>] Minimum clipped read length (default: 60)
  [--psfilterspark-args "<args>"]   Extra args passed to GATK PathSeqFilterSpark
  [--picard-jar </path/picard.jar>] Use a specific Picard JAR instead of 'picard' on PATH
  [--dont-overwrite]                Skip steps whose expected outputs already exist
  [--keep-intermediate]             Keep intermediate files (by default they are removed)

Notes:
- Outputs are written under $OUTDIR/{filter_stats,bams}.
- Required: provide HOSTDIR via --hostdir <dir> or \$HOSTDIR.
- If --sample-id is provided, missing inputs/outputs/metrics are filled from defaults.
- If --sample-id is omitted, provide explicit --input-unaligned, --paired-out, --unpaired-out, and --metrics-unaligned.
- If --sample-id is omitted and --input-decoys is provided, also provide --metrics-decoys.
- By default, pathseq-t2t qcfilter is re-run each time. Use --dont-overwrite to skip if all expected outputs already exist.
- Requires: samtools >=1.16, GATK >=4, Java 17.
HLP
        return 0 ;;
      --) shift; break ;;
      -*) die "Unknown option for qcfilter: $1" ;;
      *)  die "Unexpected argument to qcfilter: $1" ;;
    esac
  done

  # ---------- Validation & tool version checks ----------
  #require_nonempty "${input_unaligned}" "--input-unaligned"
  #require_nonempty "${input_decoys}" "--input-decoys"
  # require_file "${input_unaligned}"
  # require_file "${input_decoys}"

  _require_samtools_116
  _require_java17
  _require_gatk4
  _require_picard
  _require_hostdir "${HOSTDIR_OPT}"
  local -a psfilterspark_args_ary=()
  _split_cli_args "${psfilterspark_args}" psfilterspark_args_ary

  # Threads auto-detect
  if [[ -z "${threads}" ]]; then
    if command -v nproc >/dev/null 2>&1; then
      threads="$(nproc)"
    elif [[ "$OSTYPE" == "darwin"* ]] && command -v sysctl >/dev/null 2>&1; then
      threads="$(sysctl -n hw.ncpu)"
    else
      threads="8"
    fi
  fi

  # Ensure OUTDIR subfolders are defined even if user didn't pass --outdir explicitly
  if declare -F _set_outdirs >/dev/null; then
    _set_outdirs
  else
    OUTDIR_BAMS="${OUTDIR}/bams"
    OUTDIR_FILTER="${OUTDIR}/filter_stats"
  fi
  mkdir -p "${OUTDIR_FILTER}" "${OUTDIR_BAMS}"

  # ---------- Naming policy ----------
  # - with --sample-id: fill missing paths from defaults
  # - without --sample-id: require explicit input/output/metrics paths
  local sample_base=""
  if [[ -n "${sample_id}" ]]; then
    sample_base="${sample_id}"
    [[ -n "${input_unaligned}" ]] || input_unaligned="${OUTDIR_BAMS}/${sample_base}.prefilter.unaligned.bam"
    [[ -n "${input_decoys}"   ]] || input_decoys="${OUTDIR_BAMS}/${sample_base}.prefilter.decoys.bam"

    [[ -n "${paired_out}" ]] || paired_out="${OUTDIR_BAMS}/${sample_base}.qcfilt_paired.bam"
    [[ -n "${unpaired_out}" ]] || unpaired_out="${OUTDIR_BAMS}/${sample_base}.qcfilt_unpaired.bam"

    [[ -n "${filter_metrics_unaligned}" ]] || filter_metrics_unaligned="${OUTDIR_FILTER}/${sample_base}.prefilter.unaligned.filter_metrics.txt"
    [[ -n "${filter_metrics_decoys}" ]] || filter_metrics_decoys="${OUTDIR_FILTER}/${sample_base}.prefilter.decoys.filter_metrics.txt"
  else
    require_nonempty "${input_unaligned}" "--input-unaligned (required when --sample-id is not provided)"
    require_nonempty "${paired_out}" "--paired-out (required when --sample-id is not provided)"
    require_nonempty "${unpaired_out}" "--unpaired-out (required when --sample-id is not provided)"
    require_nonempty "${filter_metrics_unaligned}" "--metrics-unaligned (required when --sample-id is not provided)"
    if [[ -n "${input_decoys}" ]]; then
      require_nonempty "${filter_metrics_decoys}" "--metrics-decoys (required when --input-decoys is provided without --sample-id)"
    fi
    sample_base="$(basename "${paired_out%.bam}")"
  fi

  # Validate inputs: unaligned is required; decoys is optional (skip step 1A if missing)
  require_file "${input_unaligned}"
  local have_decoys=0
  if [[ -f "${input_decoys}" ]]; then
    have_decoys=1
  fi

  # Basenames used for intermediates
  local base_unaligned base_decoys base_merge
  if [[ -n "${sample_id}" ]]; then
    base_unaligned="${sample_id}.prefilter.unaligned"
  else
    base_unaligned="${sample_base}.unaligned"
  fi
  if [[ ${have_decoys} -eq 1 ]]; then
    if [[ -n "${sample_id}" ]]; then
      base_decoys="${sample_id}.prefilter.decoys"
    else
      base_decoys="${sample_base}.decoys"
    fi
  else
    base_decoys=""
  fi
  base_merge="${sample_base}"
  if [[ ${have_decoys} -eq 0 ]]; then
    filter_metrics_decoys=""
  fi

  # Derived intermediate/output paths
  local bam_input_unaligned="${input_unaligned}"
  local bam_input_decoys="${input_decoys}"
  local bam_input_decoys_rvt=""
  if [[ ${have_decoys} -eq 1 ]]; then
    bam_input_decoys_rvt="${OUTDIR_BAMS}/${base_decoys}.rvt.bam"
  fi

  local bam_paired_unaligned_filt="${OUTDIR_BAMS}/${base_unaligned}.paired.bam"
  local bam_unpaired_unaligned_filt="${OUTDIR_BAMS}/${base_unaligned}.unpaired.bam"
  local bam_paired_decoys_filt="" bam_unpaired_decoys_filt=""
  if [[ ${have_decoys} -eq 1 ]]; then
    bam_paired_decoys_filt="${OUTDIR_BAMS}/${base_decoys}.paired.bam"
    bam_unpaired_decoys_filt="${OUTDIR_BAMS}/${base_decoys}.unpaired.bam"
  fi

  local bam_paired_merge="${paired_out}"
  local bam_unpaired_merge="${unpaired_out}"

  ensure_parent_dir "${bam_paired_merge}"
  ensure_parent_dir "${bam_unpaired_merge}"

  # ---------- Failure-only cleanup of Spark parts ----------
  local PART_PAIRED_DIR="${OUTDIR_BAMS}/${base_unaligned}.qcfilt_paired.bam.parts"
  local PART_UNPAIRED_DIR="${OUTDIR_BAMS}/${base_unaligned}.qcfilt_unpaired.bam.parts"
  on_exit_qcfilter() {
    local ec=$?
    if [[ $ec -ne 0 ]]; then
      for d in "${PART_PAIRED_DIR}" "${PART_UNPAIRED_DIR}"; do
        if [[ -d "$d" ]]; then
          log "Cleaning leftover Spark parts (failure): $d"
          rm -rf "$d"
        fi
      done
    fi
  }
  trap on_exit_qcfilter EXIT

  # ---------- Pre-flight: validate input BAMs ----------
  ubam_check_or_die "${bam_input_unaligned}" "qcfilter: bam_input_unaligned"
  if [[ ${have_decoys} -eq 1 ]]; then
    bam_check_or_die "${bam_input_decoys}" "qcfilter: bam_input_decoys"
  fi

  # ---------- STEP 1A. PathSeq on DECOYS (optional) ----------
  if [[ ${have_decoys} -eq 1 ]]; then
    if [[ ${dont_overwrite} -eq 1 && -f "${filter_metrics_decoys}" && -f "${bam_paired_decoys_filt}" && -f "${bam_unpaired_decoys_filt}" ]]; then
      log "STEP 1A (--dont-overwrite): outputs exist, skipping decoys QC-filtering."
    else
      log "Reverting decoys bam"
      if [[ -n "${PICARD_JAR}" ]]; then
        java -jar "${PICARD_JAR}" RevertSam \
          --INPUT "${bam_input_decoys}" \
          --OUTPUT "${bam_input_decoys_rvt}" \
          --REMOVE_DUPLICATE_INFORMATION false \
          --RESTORE_ORIGINAL_QUALITIES false \
          --VERBOSITY ERROR
      else
        picard RevertSam \
          --INPUT "${bam_input_decoys}" \
          --OUTPUT "${bam_input_decoys_rvt}" \
          --REMOVE_DUPLICATE_INFORMATION false \
          --RESTORE_ORIGINAL_QUALITIES false \
          --VERBOSITY ERROR
      fi
      [[ -f "${bam_input_decoys_rvt}" ]] || die "RevertSam failed for decoys BAM."

      log "Running QC-filtering on decoys reads"
      local -a decoys_filter_cmd=(
        gatk --java-options "-Xmx${ram_gb}G" PathSeqFilterSpark
        --input "${bam_input_decoys_rvt}"
        --tmp-dir "${TMPDIR_OPT}"
        --spark-master "local[${threads}]"
        --bam-partition-size 0
        --is-host-aligned true
        --kmer-file "${HOSTDIR_OPT}/pathseq_host.bfi"
        --filter-bwa-image "${HOSTDIR_OPT}/pathseq_host.fa.img"
        --min-clipped-read-length "${min_clipped_read_length}"
        --filter-metrics "${filter_metrics_decoys}"
        --paired-output "${bam_paired_decoys_filt}"
        --unpaired-output "${bam_unpaired_decoys_filt}"
      )
      ((${#psfilterspark_args_ary[@]})) && decoys_filter_cmd+=("${psfilterspark_args_ary[@]}")
      log_cmd "${decoys_filter_cmd[@]}"

      [[ -f "${filter_metrics_decoys}" ]] || die "Decoys PathSeqFilterSpark failed (no metrics)."

      # Validate PathSeq outputs ONLY if their .sbi exists (i.e., tool wrote them)
      if [[ -f "${bam_paired_decoys_filt}.sbi" ]]; then
        ubam_check_or_die "${bam_paired_decoys_filt}" "qcfilter: decoys paired"
      fi
      if [[ -f "${bam_unpaired_decoys_filt}.sbi" ]]; then
        ubam_check_or_die "${bam_unpaired_decoys_filt}" "qcfilter: decoys unpaired"
      fi
    fi
  else
    log "No decoys/excluded BAM present; skipping STEP 1A."
  fi

  # ---------- STEP 1B. PathSeq on UNALIGNED ----------
  for d in "${PART_PAIRED_DIR}" "${PART_UNPAIRED_DIR}"; do
    if [[ -d "$d" ]]; then
      log "Removing stale Spark parts from previous run: $d"
      rm -rf "$d"
    fi
  done

  if [[ ${dont_overwrite} -eq 1 && -f "${bam_paired_unaligned_filt}" && -f "${bam_unpaired_unaligned_filt}" ]]; then
    log "STEP 1B (--dont-overwrite): outputs exist, skipping unaligned QC-filtering."
  else
    log "Running QC-filtering on unaligned reads"
    local -a unaligned_filter_cmd=(
      gatk --java-options "-Xmx${ram_gb}G" PathSeqFilterSpark
      --input "${bam_input_unaligned}"
      --tmp-dir "${TMPDIR_OPT}"
      --spark-master "local[${threads}]"
      --bam-partition-size 4000000
      --is-host-aligned true
      --kmer-file "${HOSTDIR_OPT}/pathseq_host.bfi"
      --filter-bwa-image "${HOSTDIR_OPT}/pathseq_host.fa.img"
      --min-clipped-read-length "${min_clipped_read_length}"
      --filter-metrics "${filter_metrics_unaligned}"
      --paired-output "${bam_paired_unaligned_filt}"
      --unpaired-output "${bam_unpaired_unaligned_filt}"
    )
    ((${#psfilterspark_args_ary[@]})) && unaligned_filter_cmd+=("${psfilterspark_args_ary[@]}")
    time "${unaligned_filter_cmd[@]}"

    if [[ -f "${bam_paired_unaligned_filt}.sbi" ]]; then
      ubam_check_or_die "${bam_paired_unaligned_filt}"   "qcfilter: unaligned paired"
    else
      die "Expected qcfilter output missing .sbi: ${bam_paired_unaligned_filt}.sbi"
    fi

    if [[ -f "${bam_unpaired_unaligned_filt}.sbi" ]]; then
      ubam_check_or_die "${bam_unpaired_unaligned_filt}" "qcfilter: unaligned unpaired"
    else
      die "Expected qcfilter output missing .sbi: ${bam_unpaired_unaligned_filt}.sbi"
    fi
  fi

  # Collect final number of paired/unpaired reads
  if [[ -n "${filter_metrics_unaligned:-}" && -f "${filter_metrics_unaligned}" ]]; then
    IFS=$'\t' read -r unaln_final_paired unaln_final_unpaired < <(
      awk 'NF && $1 ~ /^[0-9]+$/ {print $6 "\t" $7; exit}' "${filter_metrics_unaligned}"
    )
    log "Unaligned metrics: FINAL_PAIRED_READS=${unaln_final_paired:-0} FINAL_UNPAIRED_READS=${unaln_final_unpaired:-0}"
  fi


  # Read FINAL_PAIRED/UNPAIRED_READS from decoy filter metrics (if present)
  local decoys_final_paired=0 decoys_final_unpaired=0
  if [[ ${have_decoys} -eq 1 && -n "${filter_metrics_decoys:-}" && -f "${filter_metrics_decoys}" ]]; then
    read -r decoys_final_paired decoys_final_unpaired < <(
      awk 'NF && $1 ~ /^[0-9]+$/ {print $6, $7; exit}' "${filter_metrics_decoys}"
    )
    log "Decoy metrics: FINAL_PAIRED_READS=${decoys_final_paired} FINAL_UNPAIRED_READS=${decoys_final_unpaired}"
  fi


  # ---------- STEP 1C. Merge DECOYS + UNALIGNED ----------
  log "Merging decoys + unaligned outputs"
  

  # ----- Merge Paired -----
  if [[ -f "${bam_paired_unaligned_filt}" ]]; then
    if [[ -n "${bam_paired_decoys_filt}" && -f "${bam_paired_decoys_filt}" ]]; then
      samtools merge -@ "${threads}" "${bam_paired_merge}" \
        "${bam_paired_unaligned_filt}" "${bam_paired_decoys_filt}"
    else
      if [[ ${have_decoys} -eq 1 && ${decoys_final_paired} -gt 0 ]]; then
        die "Decoy metrics FINAL_PAIRED_READS=${decoys_final_paired} but paired decoy BAM missing: ${bam_paired_decoys_filt}"
      fi
      cp -f "${bam_paired_unaligned_filt}" "${bam_paired_merge}"
    fi
  else
    die "Unaligned metrics FINAL_PAIRED_READS=${unaln_final_paired} but paired unaligned BAM missing: ${bam_paired_unaligned_filt}"
  fi


  # ----- Merge Unpaired -----
  if [[ -f "${bam_unpaired_unaligned_filt}" ]]; then
    if [[ -n "${bam_unpaired_decoys_filt}" && -f "${bam_unpaired_decoys_filt}" ]]; then
      samtools merge -@ "${threads}" "${bam_unpaired_merge}" \
        "${bam_unpaired_unaligned_filt}" "${bam_unpaired_decoys_filt}"
    else
      if [[ ${have_decoys} -eq 1 && ${decoys_final_unpaired} -gt 0 ]]; then
        die "Decoy metrics FINAL_UNPAIRED_READS=${decoys_final_unpaired} but unpaired decoy BAM missing: ${bam_unpaired_decoys_filt}"
      fi
      cp -f "${bam_unpaired_unaligned_filt}" "${bam_unpaired_merge}"
    fi
  else
    die "Unaligned metrics FINAL_UNPAIRED_READS=${unaln_final_unpaired} but unpaired unaligned BAM missing: ${bam_unpaired_unaligned_filt}"
  fi


  # Validate merged outputs (merged BAMs won’t have .sbi; use generic check)
  ubam_check_or_die "${bam_paired_merge}"   "qcfilter: merged paired"
  ubam_check_or_die "${bam_unpaired_merge}" "qcfilter: merged unpaired"

  # ---------- Post-success cleanup ----------
  if [[ ${keep_intermediate} -eq 0 ]]; then
    log "Removing intermediate outputs (use --keep-intermediate to retain):"
    for f in \
      ${bam_input_decoys_rvt:+${bam_input_decoys_rvt}} \
      ${bam_paired_decoys_filt:+${bam_paired_decoys_filt}} \
      ${bam_unpaired_decoys_filt:+${bam_unpaired_decoys_filt}} \
      "${bam_paired_unaligned_filt}" \
      "${bam_paired_unaligned_filt}.sbi" \
      "${bam_unpaired_unaligned_filt}" \
      "${bam_unpaired_unaligned_filt}.sbi"
    do
      [[ -n "$f" && -e "$f" ]] && rm -f "$f"
    done
  fi

  log "qcfilter done → paired: ${bam_paired_merge} ; unpaired: ${bam_unpaired_merge}"
}
